#!/bin/bash
# ------------------------------------------------------------------------------------
# Run script for lsm at FLUXNET sites
# This script calls LSMSRunArray.sh for all sites given in the file "site_list"
# This script assumes that the SITE_LIST and LSMSRunArray.sh file exist
#  in the same directory
# ------------------------------------------------------------------------------------
# Settings:
NPROC=32                                      # number of CPUs to be used
SITE_LIST=all_sites_list.txt                  # site list file to be used
FORCING_TSTEP=TIMESTEP                        # time step of forcing (TIMESTEP or DAILY)
CODEVERSION=rev272                            # name extension of executable
EXP_ID=${CODEVERSION}                         # name of current experiment
RESTART=0                                     # whether to use existing restart file or not
RESTART_ID=${CODEVERSION}                     # name of experiment to restart from
LOCALTIME=TRUE                                # whether forcing timeaxis is local (TRUE) 
COMPILER=x86_64-pgi
                                              # or GMT (FALSE)
#PATHs
OUT_PATH=/public1/home/luolh/CLM4.0/Output
FORCING_PATH=/public1/home/luolh/FLUXNET/Forcing_data
OBS_PATH=/public1/home/luolh/FLUXNET/Obs_data

# set experiment id and code version to general_settings.txt
echo ${EXP_ID} ${CODEVERSION} ${RESTART} ${RESTART_ID} ${FORCING_TSTEP} ${LOCALTIME} > ${PWD}/general_settings.txt
echo ${CODE_BASE_PATH} ${COMPILER} ${FORCING_PATH} ${OBS_PATH} ${OUT_PATH} > ${PWD}/paths.txt

# cp SITE_LIST to site_list.txt file
if [[ ${SITE_LIST} != 'site_list.txt' ]] ; then
 cp ${SITE_LIST} site_list.txt
fi

NSITE=$(wc -l site_list.txt | cut -c1-8 | awk '{print $1}')
# array starts at 2 because first line in site_list.txt is a header
bsub -J lsm_s[2-${NSITE}]%${NPROC} -e run_lsm_site.%I.e -o run_lsm_site.%I.o < ./LSMSRunArray.sh
bsub -w "ended(lsm_s[2-${NSITE}])" -J cleanup -e run_lsm_cleanup.e -o run_lsm_cleanup.o < ./LSMSRunCLeanup.sh
