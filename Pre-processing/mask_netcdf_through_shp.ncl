;----------------------------------------------------------------------
;
; Concepts illustrated:
;   - Using a worldwide shapefile to create a land/ocean mask
;   - Using "mask" to set land or ocean values in your data to missing
;   - Masking a data array based on a geographical area
;   - Attaching shapefile polylines to a map plot
;   - Adding a _FillValue attribute to a variable

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;----------------------------------------------------------------------
; Function : create_mask_from_shapefile
;
; This function takes a shapefile of type "polygon" and a 2D data 
; array that has coordinate arrays, and creates a mask array 
; the same size as the data, that contains 0s and 1s, depending if
; the lat/lon values of the data are inside (1) or outside (0) of the 
; shapefile polygons.
;----------------------------------------------------------------------
undef("create_mask_from_shapefile")
function create_mask_from_shapefile(data[*][*]:numeric,fname[1]:string)
local f, segments, geometry, segsDims, geomDims, geom_segIndex, \
geom_numSegs, segs_xyzIndex, segs_numPnts, numFeatures, ilat, ilon, i, j, \
shp_lat, shp_lon, data_lat, data_lon, startSegment, numSegments, seg, \
startPT, endPT, mask_val
begin
;---Error checking
  if(.not.isfilepresent(fname)) then
    print("Error: create_mask_from_shapefile:")
    print("      '" + fname + "' doesn't exist.")
    print("       Mask array with all missing values will be returned.")
    return(new(dimsizes(data),integer,-999))
  end if

;---Check that "data" has 1D coordinate arrays.
  if(.not.isdimnamed(data,0).or..not.isdimnamed(data,1).or.\
     .not.iscoord(data,data!0).or..not.iscoord(data,data!1)) then
    print("Error: create_mask_from_shapefile:")
    print("       Input data doesn't have 1D coordinate arrays")
    print("       Mask array with all missing values will be returned.")
    return(new(dimsizes(data),integer,-999))
  end if

;---Open the shapefile
  f = addfile(fname,"r")

;---We can't use this routine to mask against point or line data
  if(f@geometry_type.ne."polygon") then
    print("Error: create_mask_from_shapefile: geometry_type attribute must be 'polygon'")
    print("       Mask array with all missing values will be returned.")
    return(new(dimsizes(data),integer,-999))
  end if

;---Read data off the shapefile
  segments = f->segments
  geometry = f->geometry
  segsDims = dimsizes(segments)
  geomDims = dimsizes(geometry)

;---Read global attributes  
  geom_segIndex = f@geom_segIndex
  geom_numSegs  = f@geom_numSegs
  segs_xyzIndex = f@segs_xyzIndex
  segs_numPnts  = f@segs_numPnts
  numFeatures   = geomDims(0)

;---Create mask array
  data_mask = new(dimsizes(data),integer,-999)
  mask_val  = 1     ; 1's represent being inside a polygon
  data_mask = 0     ; 0's represent being outside a polygon

;---Read lat/lon values
  shp_lon  = f->x
  shp_lat  = f->y
  data_lat = data&$data!0$
  data_lon = data&$data!1$
  nlat = dimsizes(data_lat)
  nlon = dimsizes(data_lon)

  do i=0, numFeatures-1  
    startSegment = geometry(i, geom_segIndex)
    numSegments  = geometry(i, geom_numSegs)
    do seg=startSegment, startSegment+numSegments-1
      startPT = segments(seg, segs_xyzIndex)
      endPT   = startPT + segments(seg, segs_numPnts) - 1
;
; Messy code to make the nested loop go faster.  Only check
; the area of data that is close to polygon that we are
; traversing.  If you check data's full lat/lon array against
; every shapefile polygon, this loop will be extremely slow.
;

;
; Get approx. index values where data's lat/lon is close
; to shapefile's current polygon.
;
      iilt_beg = ind(data_lat.le.min(shp_lat(startPT:endPT)))
      iilt_end = ind(data_lat.ge.max(shp_lat(startPT:endPT)))
      iiln_beg = ind(data_lon.le.min(shp_lon(startPT:endPT)))
      iiln_end = ind(data_lon.ge.max(shp_lon(startPT:endPT)))
      ilt_beg = 0
      iln_beg = 0
      ilt_end = nlat-1
      iln_end = nlon-1
      if(.not.any(ismissing(iilt_beg))) then
        ilt_beg = iilt_beg(dimsizes(iilt_beg)-1)
      end if
      if(.not.any(ismissing(iilt_end))) then
        ilt_end = iilt_end(0)
      end if
      if(.not.any(ismissing(iiln_beg))) then
        iln_beg = iiln_beg(dimsizes(iiln_beg)-1)
      end if
      if(.not.any(ismissing(iiln_end))) then
        iln_end = iiln_end(0)
      end if
;
; Loop across subset of data's lat/lon and check each point
; to see if it is inside or outside of the shapefile polygon.
;
      do ilat = ilt_beg,ilt_end
        do ilon = iln_beg,iln_end
          if(data_mask(ilat,ilon).ne.mask_val.and.\
            gc_inout(data_lat(ilat),data_lon(ilon),\
                     shp_lat(startPT:endPT),shp_lon(startPT:endPT))) then
            data_mask(ilat,ilon) = mask_val
          end if
        end do
      end do
      delete([/iilt_beg,iilt_end,iiln_beg,iiln_end/])
    end do
  end do
  return(data_mask)
end

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
  WRITE_MASK = True     ; Whether to write shapefile mask to file

;---Open file containing data to mask
  dir        = "./"
  cdf_prefix = "prec_ITPCAS-CMFD_V0106_B-01_201208"
  cdf_file   = dir + cdf_prefix + ".nc"
  fin        = addfile(cdf_file,"r")

;---Read "ts" and corresponding "ORO" mask.
  ts            = fin->prec(0,:,:)
  ts@_FillValue = default_fillvalue(typeof(ts))
  
;
; Create copy of "ts" with longitudes flipped. This
; is necessary for applying the shapefile mask, b/c
; the shapefile longitudes go from -180 to 180.
;
  ts_flip = lonFlip(ts)
;
; Create a mask array the same size as "ts", using
; lat/lon data read off a shapefile.
;
  shpfile  = "QTP_railway_highway_buffer30km.shp"
  shp_mask = create_mask_from_shapefile(ts_flip,shpfile)
;
; Mask "ts" using "ORO" mask on file and shapefile
; land mask, for comparison.
;
  ts_shp_mask = where(shp_mask.eq.1,ts_flip,ts_flip@_FillValue)
  
  copy_VarMeta(ts_flip,ts_shp_mask)
  
  time  =  fin->time 
  lat   =  fin->lat                   
  lon   =  fin->lon                   

ntim  =  dimsizes(time)
  nlat  =  dimsizes(shp_mask)               ;-- get dimension sizes of lat
  nlon  =  dimsizes(shp_mask)               ;-- get dimension sizes of lon


    new_cdf_file = cdf_prefix + "_with_mask.nc"
  
if (isfilepresent(new_cdf_file)) then    
     system("/bin/rm " + new_cdf_file)
end if 
   finout = addfile(new_cdf_file,"c")


setfileoption(finout,"DefineMode",True) 

  fAtt                  =  True        ;-- assign file attributes
  fAtt@title            = "NCL Efficient Approach to netCDF Creation"  
  fAtt@source_file      = "test.nc"
  fAtt@Conventions      = "CF"   
  fAtt@creation_date    =  systemfunc ("date")        
  fAtt@history          =  "NCL script: ex_netcdf_2.ncl" 
  fAtt@comment          = "Convert variable t: Kelvin to Celsius"       
  fileattdef(finout,fAtt)                ;-- copy file attributes    

  dimNames = (/"time", "lat", "lon"/)  
  dimSizes = (/ -1   ,  nlat,  nlon /) 
  dimUnlim = (/ True ,  False, False/)   
  filedimdef(finout,dimNames,dimSizes,dimUnlim)

  filevardef(finout, "time" ,typeof(time),getvardims(time))        
  filevardef(finout, "lat"  ,typeof(lat), getvardims(lat))                          
  filevardef(finout, "lon"  ,typeof(lon), getvardims(lon))                          
  filevardef(finout, "TS"   ,typeof(ts),  getvardims(ts))

  filevarattdef(finout,"time" ,time)       ;-- copy time attributes
  filevarattdef(finout,"lat"  ,lat)        ;-- copy lat attributes
  filevarattdef(finout,"lon"  ,lon)        ;-- copy lon attributes
  filevarattdef(finout,"TS",   ts)         ;-- copy tC attributes

  setfileoption(finout,"DefineMode",False)

  finout->time   =  (/time/)               ;-- write time to new netCDF file
  finout->lat    =  (/lat/)                ;-- write lat to new netCDF file
  finout->lon    =  (/lon/)                ;-- write lon to new netCDF file

    finout->TS = (/ts_shp_mask/)
  ;end if
end
