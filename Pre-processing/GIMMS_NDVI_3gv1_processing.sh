# GIMMS NDVI 3gv1 data processing
#!/bin/bash

mkdir GIMMS_Summer_Output

for year in `seq 1982 2010`
do

cdo seldate,0000-00-05,0000-00-06 GIMMS_NDVI_3gv1/ndvi3g_geo_v1_${year}_0106.nc4 ndvi3g_geo_v1_${year}_MJ.nc4
cdo seldate,0000-00-07,0000-00-08 GIMMS_NDVI_3gv1/ndvi3g_geo_v1_${year}_0712.nc4 ndvi3g_geo_v1_${year}_JA.nc4
cdo seldate,0000-00-09 GIMMS_NDVI_3gv1/ndvi3g_geo_v1_${year}_0712.nc4 ndvi3g_geo_v1_${year}_S.nc4

cdo mergetime ndvi3g_geo_v1_${year}_MJ.nc4 ndvi3g_geo_v1_${year}_JA.nc4 ndvi3g_geo_v1_${year}_S.nc4 ndvi3g_geo_v1_${year}_MJJAS.nc4

rm ndvi3g_geo_v1_${year}_MJ.nc4
rm ndvi3g_geo_v1_${year}_JA.nc4
rm ndvi3g_geo_v1_${year}_S.nc4

cdo timmean ndvi3g_geo_v1_${year}_MJJAS.nc4 ndvi3g_geo_v1_${year}_MJJAS_mean.nc4

gdal_translate -of GTiff NetCDF:"ndvi3g_geo_v1_${year}_MJJAS_mean.nc4":ndvi ndvi3g_geo_v1_${year}_MJJAS_mean.TIF

rm ndvi3g_geo_v1_${year}_MJJAS.nc4

mv ndvi3g_geo_v1_${year}_MJJAS_mean.nc4 GIMMS_Summer_Output/ndvi3g_geo_v1_${year}_MJJAS_mean.nc4

done
