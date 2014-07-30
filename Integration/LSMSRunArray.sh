#!/bin/bash
#------------------------------------------------------------------------------------
# RUNSCRIPT TO RUN lsm for FLUXNET sites
#------------------------------------------------------------------------------------
#BSUB -q io 
#BSUB -a mpich2pgi
#BSUB -n 1
if [ $LSB_JOBINDEX -eq 0 ] ; then
  echo script needs to be queued as array
  return
fi
#------------------------------------------------------------------------------------
username=luolh
#------------------------------------------------------------------------------------
# 1. EXPERIMENT SETTINGS
#------------------------------------------------------------------------------------
# 1.1 global settings
#
# name of experiment
exp_id=$(head -1 general_settings.txt | tail -1 | awk '{print $1}')
# code version to be used
codeversion=$(head -1 general_settings.txt | tail -1 | awk '{print $2}')
# whether to use existing restart (1) or not (0)
restart=$(head -1 general_settings.txt | tail -1 | awk '{print $3}')
# name of experiment to restart from  
restart_id=$(head -1 general_settings.txt | tail -1 | awk '{print $4}')
# forcing version used (TIMESTEP or DAILY)  
forcing_frequency=$(head -1 general_settings.txt | tail -1 | awk '{print $5}')
# whether time axis in forcing is local time or GMT
LOCALTIME=$(head -1 general_settings.txt | tail -1 | awk '{print $6}')
#
# 1.2 paths
#
# directory in which the code is located
code_dir=$(head -1 paths.txt | tail -1 | awk '{print $1}')
# compiler sub-directory for linking
compiler=$(head -1 paths.txt | tail -1 | awk '{print $2}')
# directory in which the model forcing is located
forcing_dir=$(head -1 paths.txt | tail -1 | awk '{print $3}')
# directory in which the observed fluxes are located
obs_dir=$(head -1 paths.txt | tail -1 | awk '{print $4}')
# directory for rerun and output storage
local_path=$(head -1 paths.txt | tail -1 | awk '{print $5}')
#
# 1.3 local settings
#
# site to be simulated
SITE=$(head -$LSB_JOBINDEX site_list.txt | tail -1 | awk '{print $1}') 
# start of experiment
STARTYEAR=$(head -$LSB_JOBINDEX site_list.txt | tail -1 | awk '{print $2}')
# end of experiment
ENDYEAR=$(head -$LSB_JOBINDEX site_list.txt | tail -1 | awk '{print $3}')
# number of tiles
NTILES=$(head -$LSB_JOBINDEX site_list.txt | tail -1 | awk '{print $4}')
# shift local time by X hours to compsenate for local-time time-zone time difference
TIMEOFFSET=$(head -$LSB_JOBINDEX site_list.txt | tail -1 | awk '{print $5}')

# data, work path etc.
out_dir=${local_path}/${SITE}
restart_dir=${local_path}/${SITE}
username=$(whoami)
scratch_dir=/scratch/${username}/TMP_lsm


#------------------------------------------------------------------------------------
# 1.1. RUN OPTIONS 
#------------------------------------------------------------------------------------
postprocess_flag=1         # whether to do postprocessing (1) or not (0)
NYEARINIT=15               # number of initial years for spin-up forcing generation
                           # (should be >10years)
nyear_cbalone=1000         # years of carbon spin-up (should be 1000)
writefreq_out=1800         # final output frequency in seconds 
out_filetype=GRIB          # format of final output (NETCDF or GRIB)
debug=FALSE                # write debug information 
dynveg=FALSE               # activate dynamic vegetation
read_fpc=FALSE             # take fpc from file 
withnitrogen=FALSE         # activate nitrogen cycling
read_ndepo=FALSE           # read Ndeposition data?

#------------------------------------------------------------------------------------
# 2. SETUP OF EXPERIMENT: nothing to be changed beyond this point 
#------------------------------------------------------------------------------------
startdate=${STARTYEAR}0101       # start of run
let date=${ENDYEAR}+1            
enddate=${date}0101              # end of run
nsdt=1800                        # time step in seconds
writefreq_spin=86400             # output during spin-up (has to be daily !!) in seconds

#-- link in script library
fpath=${code_dir}/contrib/site_scripts/
export PATH=$fpath:${PATH}
. lib_scripts.ksh
export cdo=/usr/local/bin/cdo-1.4.6

#-- define work directory and create if not existing
if [[ $LSB_JOBID = '' ]] ; then
 CWD=$PWD/work
else
 CWD=${scratch_dir}_${LSB_JOBID}_${LSB_JOBINDEX}
 #mrun=mpirun_pgi10.lsf
fi
if [ ! -d ${CWD} ] ; then
      mkdir -p ${CWD}
fi  
cd ${CWD}

#-- create output directory if not existing
if [ ! -d ${out_dir} ] ; then
      mkdir -p ${out_dir}/Restart
      mkdir -p ${out_dir}/Output
fi  

#-- cp input files to work directory
case "$forcing_frequency" in
    TIMESTEP)
	ftype=halfhourly
    ;;
    DAILY)
        ftype=daily
    ;;
    *)
    echo invalid forcing frequency selected, only TIMESTEP and DAILY allowed
    return
    ;;
esac
    
#-- copy input data and executables to woring directory
cp ${forcing_dir}/${ftype}/${SITE}.${STARTYEAR}-${ENDYEAR}.forcing.${ftype}.nc forcing.nc
cp ${forcing_dir}/surface/${SITE}.${NTILES}_tiles.surface.nc lsm.nc
cp ${code_dir}/util/running/adjunct_files/lsm/lctlib_nlct21.def lctlib.def
cp ${code_dir}/${compiler}/bin/lsm_${codeversion}.x lsm.x
cp ${code_dir}/${compiler}/bin/cbalance_${codeversion}.x cbalance.x
if [ $postprocess_flag -eq 1 ] ; then
 cp ${code_dir}/contrib/site_scripts/generate_FLUXNET_diagnostics.R .
fi
#------------------------------------------------------------------------------------
# 3. INITIALISATION OF EXPERIMENT from scratch if not reading restart files 
#------------------------------------------------------------------------------------

#--if this run starts from scratch...
if [ $restart -eq 0 ] ; then

#-- create basic namelist
setup_lsm_namelist ${forcing_frequency} ${NTILES} ${LOCALTIME} ${TIMEOFFSET} GRIB ${debug} \
                      ${dynveg} ${read_fpc} ${withnitrogen}

#-- create namelist for initial run 
read_cpool=FALSE
read_npool=FALSE
restart=FALSE
update_lsm_namelist ${code_dir} ${exp_id} ${startdate} ${enddate} ${nsdt} \
   ${writefreq_spin} ${restart} ${read_cpool} ${read_npool} ${read_ndepo}
set_stream_elements ${code_dir} LAI_yDayMean NPP_yDayMean topSoilTemp_yDayMean \
   alpha_yDayMean

#------------------------------------------------------------------------------------
# 3.1 run lsm to accumulate phenological triggers and equilibriate water cycle 
#     this provides the input for C cycle equilibration 
#------------------------------------------------------------------------------------
# 3.1.1 work out how many cycles are required to obtain NYEARINIT years
let nyear=$ENDYEAR-$STARTYEAR+1
NINIT=$((${NYEARINIT}/${nyear}))
[ $NINIT -lt 2 ] && NINIT=2
i=1
while [ $i -le $NINIT ] ; do

  #-- run the model
  time ${mrun} ${PWD}/lsm.x
  #-- modifiy date stamp of restart file 
  let setdate=${STARTYEAR}-1
  setdate=${setdate}1231
  touch_rerun_files ${exp_id} ${setdate}

  if [ $i -eq 1 ] ; then
    #-- in second and higher runs use restart file and update namelist 
    restart=TRUE
    startdate=00000000
    update_lsm_namelist ${code_dir} ${exp_id} ${startdate} ${enddate} ${nsdt} \
      ${writefreq_spin} ${restart} ${read_cpool} ${read_npool} ${read_ndepo}
    set_stream_elements ${code_dir} LAI_yDayMean NPP_yDayMean topSoilTemp_yDayMean \
      alpha_yDayMean
  fi

  let i=i+1 
done
#------------------------------------------------------------------------------------
# 3.2 prepare forcing for carbon spinup 
#------------------------------------------------------------------------------------
for i in `ls \${exp_id}*_veg`; do
 j=${i#${exp_id}_}
 $cdo -s selname,var165,var166,var168,var169 $i VEG.TMP
 $cdo -s -f nc -r -chvar,var165,LAI_yDayMean -chvar,var166,NPP_yDayMean \
                  -chvar,var168,topSoilTemp_yDayMean -chvar,var169,alpha_yDayMean \
                  VEG.TMP jsb_cb_${j:0:6}.lsm_yDay_Mean.nc
 rm VEG.TMP 
done

#------------------------------------------------------------------------------------
# 3.3 run cbalone to equilibrise carbon cycle 
#------------------------------------------------------------------------------------
create_cbalone_namelist $STARTYEAR ${ENDYEAR} 1 ${nyear_cbalone} ${dynveg} \
    ${read_fpc} ${withnitrogen} ${read_cpool} ${read_npool} ${read_ndepo}
time ${mrun} ${PWD}/cbalance.x

typeset -RZ6 i=${nyear_cbalone} 
cp Cbalone.${i}.????.nc Cpools.nc

#save restart and Cpools
save_rerun_files_simple ${exp_id} ${out_dir} ${SITE} ${NTILES} ${restart_id}

else 
 #-- if using old restart files, get them from repository
 get_rerun_files_simple ${exp_id} ${restart_dir} ${SITE} ${NTILES} ${restart_id} 
fi

#------------------------------------------------------------------------------------
# 4. run lsm for actual experiment 
#------------------------------------------------------------------------------------

#-- create namelist
setup_lsm_namelist ${forcing_frequency} ${NTILES} ${LOCALTIME} ${TIMEOFFSET} \
                      ${out_filetype} ${debug} ${dynveg} ${read_fpc} ${withnitrogen}
read_cpool=TRUE
if [[ ${withnitrogen} = 'TRUE' ]] ; then 
  read_npool=TRUE
else
  read_npool=FALSE
fi  
restart=TRUE
startdate=00000000
update_lsm_namelist ${code_dir} ${exp_id} ${startdate} ${enddate} ${nsdt} \
   ${writefreq_out} ${restart} ${read_cpool} ${read_npool} ${read_ndepo}
set_stream_elements ${code_dir} sensible_heat_flx net_radiation surface_temperature \
   latent_heat_flx transpiration evapotranspiration canopy_cond_limited canopy_conductance \
   soil_moisture runoff lai apar_acc par_acc CO2_conc_leaf net_assimilation \
   boxC_litter_wood boxC_green boxC_woods boxC_reserve boxC_litter_green_bg \
   boxC_slow boxC_litter_green_ag zCO2_flux_net box_soil_respiration \
   NPP_act_yDayMean air_temp precip_rain precip_snow spec_humidity frac_PAR_diffuse \
   rad_sw_down_pot wind_speed zchl_acc cdrag_acc
 
time ${mrun} ${PWD}/lsm.x

#------------------------------------------------------------------------------------
# 5. postprocessing 
#------------------------------------------------------------------------------------
set -vx
if [ $postprocess_flag -eq 0 ] ; then 
  mv ${exp_id}*nc $out_dir
else
  # join output files and generate additional output variables
  merge_jsb_output ${exp_id} ${out_filetype}
  generate_additional_vars 

  #aggregate in time
  if [ $writefreq_out -eq 1800 ] ; then
    mv lsm.all.nc lsm.halfhourly.nc
    cp lsm.halfhourly.nc ${out_dir}/Output/${SITE}.${exp_id}.${ftype}.lsm.halfhourly.nc
    $cdo -s -daymean lsm.halfhourly.nc lsm.daily.nc
  fi
  if [ $writefreq_out -eq 86400 ] ; then
    mv lsm.all.nc lsm.daily.nc
  fi
  $cdo -s -monmean lsm.daily.nc lsm.monthly.nc
  $cdo -s -yearmean lsm.daily.nc lsm.annual.nc

  cp lsm.daily.nc ${out_dir}/Output/${SITE}.${exp_id}.${ftype}.lsm.daily.nc
  cp lsm.monthly.nc ${out_dir}/Output/${SITE}.${exp_id}.${ftype}.lsm.monthly.nc
  cp lsm.annual.nc ${out_dir}/Output/${SITE}.${exp_id}.${ftype}.lsm.annual.nc
#-----------------------------------------------------------------------------------
# 5.1 statistics and visualisation  
#-----------------------------------------------------------------------------------
  cp ${obs_dir}/halfhourly/${SITE}.${STARTYEAR}-${ENDYEAR}.obs.halfhourly.nc obs.halfhourly.nc

  ./LSMSDiagnostics.R ${STARTYEAR} ${ENDYEAR} 

  cp Rplots.pdf ${out_dir}/Output/${SITE}.${exp_id}.${ftype}.plots.pdf
  cp summary.annual.txt ${out_dir}/Output/${SITE}.${exp_id}.${ftype}.summary.txt
  cp statistics.halfhourly.txt ${out_dir}/Output/${SITE}.${exp_id}.${ftype}.stats.halfhourly.txt
  cp statistics.daily.txt ${out_dir}/Output/${SITE}.${exp_id}.${ftype}.stats.daily.txt
  cp statistics.monthly.txt ${out_dir}/Output/${SITE}.${exp_id}.${ftype}.stats.monthly.txt
fi

#------------------------------------------------------------------------------------
# 6. clean up 
#------------------------------------------------------------------------------------
if [[ $LSB_JOBID = '' ]] ; then
 echo simulation terminated 
else
 rm -R $CWD
 echo simulation terminated 
fi
