#!/bin/bash
# Date: 24/02/2012
# Written by Lihui Luo 
# E-mail: luolh@lzb.ac.cn
# Conver the output file of WRF model to the Forcing data set of Land Surface Models (Input file)

echo "requires the year of data as input"

# set the year of process 
export YEAR=$1

#delete some netcdf file after run the script
rm ${YEAR}-*.del.nc
rm ${YEAR}-*.ok.nc
rm ${YEAR}-*.all.nc
rm ${YEAR}-*.PRECTmms.nc

# process rainc
# month loop
for i in 01 02 03 04 05 06 07 08 09 10 11 12
do

    # get data from Pan Xiaoduo Directory, and process some variables
    # need to modify the data directory and the file name format of NetCDF data
    cdo selvar,Q2,T2,U10,V10,PSFC,RAINC,SWDOWN ~panxiaoduo/QZ/WRFV3/test/em_real/wrfout_d02_${YEAR}-${i}-01_00:00:00 ${YEAR}-${i}-wrfout.nc
    cdo expr,'WIND=sqrt(sqr(U10)+sqr(V10))' ${YEAR}-${i}-wrfout.nc ${YEAR}-${i}.wind.nc
    cdo selvar,RAINC ${YEAR}-${i}-wrfout.nc ${YEAR}-${i}.rainc.nc

    # forward to the first hour
    cdo shifttime,-1hours ${YEAR}-${i}.rainc.nc ${YEAR}-${i}.rainc.del.nc
    
    # delete the other month
    cdo selmon,${i} ${YEAR}-${i}.rainc.del.nc ${YEAR}-${i}.rainc.to.nc

    # process of WIND
    cdo shifttime,-1hours ${YEAR}-${i}.wind.nc ${YEAR}-${i}.wind.del.nc
    cdo selmon,${i} ${YEAR}-${i}.wind.del.nc ${YEAR}-${i}.wind.ok.nc

    # modify variable names
    cdo chname,Q2,QBOT,T2,TBOT,PSFC,PSRF,SWDOWN,FSDS ${YEAR}-${i}-wrfout.nc ${YEAR}-${i}.chname.nc
    cdo shifttime,-1hours ${YEAR}-${i}.chname.nc ${YEAR}-${i}.chname.del.nc
    cdo selmon,${i} ${YEAR}-${i}.chname.del.nc ${YEAR}-${i}.chname.ok.nc

    # hour
    for j in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23
    do
        #convert unit mm to mm/s
        cdo seltime,${j}:00:00 ${YEAR}-${i}.rainc.to.nc ${YEAR}-${i}.rainc.to.${j}.nc
        cdo expr,"RAINC=RAINC/(60*60*(${j}+1))" ${YEAR}-${i}.rainc.to.${j}.nc ${YEAR}-${i}.rainc.to.${j}.expr.nc
    done

    # merge all timestep
    cdo mergetime ${YEAR}-${i}.rainc.to.*.expr.nc ${YEAR}-${i}.rainc.ok.nc

    # delete redundancy netcdf file
    rm ${YEAR}-${i}.rainc.del.nc
    rm ${YEAR}-${i}.rainc.to*.nc
    rm ${YEAR}-${i}.wind.del.nc
    rm ${YEAR}-${i}.chname.del.nc

    # change RAINC to PRECTmms
    cdo chname,RAINC,PRECTmms ${YEAR}-${i}.rainc.ok.nc ${YEAR}-${i}.rainc.PRECTmms.nc
    rm ${YEAR}-${i}.rainc.ok.nc

    # merge PRECTmms, WIND, CHNAME to ALL
    cdo merge ${YEAR}-${i}.chname.ok.nc ${YEAR}-${i}.rainc.PRECTmms.nc ${YEAR}-${i}.wind.ok.nc ${YEAR}-${i}.all.nc
	
    # modify the attributions
    ncdump ${YEAR}-${i}-all.nc | sed -e "s/x/lon/g" | sed -e "s/y/lat/g" | sed -e "s/XLAT/LATIXY/g" | sed -e "s/XLONG/LONGXY/g" |ncgen -o ${YEAR}-${i}.nc

    # create the data directory of new LSMS Forcing data sets	
    mkdir LSMS_Forcing
    mv ${YEAR}-${i}.nc LSMS_Forcing/${YEAR}-${i}.nc
done
