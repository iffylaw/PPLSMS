	program dataprocess
	
	implicit none

	character(len=30) sitename
	double precision lat, lon, ht
      integer firstyear
	integer  i, N_out, thisnpd, lst, filltype
	integer nobsmax
	parameter (nobsmax=17520)
	double precision Tair_out(nobsmax), RH_out(nobsmax)
      double precision wsp_out(nobsmax), precip_out(nobsmax)
      double precision press_out(nobsmax), Rg_in_out(nobsmax)
      double precision Rlong_in_out(nobsmax), dsince(nobsmax)
	integer year(nobsmax),month(nobsmax),day(nobsmax)
	integer hour(nobsmax),minute(nobsmax)
	character(len=100) fname
      character(len=4) yst

!       N_out: Number of data point
!       thisnpd: Number of data points in one day
!       lst: local standard time 

	sitename='HBG'
        lat=37.67
        lon=101.33
        ht=2.2
c	firstyear=2004
        print *,'Please input the first year'
        read(*,*)firstyear
        N_out=17520
        thisnpd=48
	lst=8
	filltype=0
      write(yst, '(I4)')    firstyear

	fname=trim(sitename)//'METE'//trim(yst)//'.txt'
        print *,fname
	open(1,file=fname)
		do i=1,nobsmax
			read(1,*)year(i),month(i),day(i),hour(i),minute(i),
     &			Rlong_in_out(i),Rg_in_out(i),precip_out(i),
     &			press_out(i),RH_out(i),Tair_out(i),wsp_out(i)
			press_out(i)=press_out(i)/10
		enddo
	close(1)
	
        call write_clm(sitename, lat, lon, ht, firstyear, N_out, 
     &     Tair_out, RH_out, wsp_out, precip_out, press_out, 
     &     Rg_in_out, Rlong_in_out, thisnpd, lst, filltype)

	stop
	end

      subroutine write_clm(sitename, lat, lon, ht, firstyear, N_out, 
     &     Tair_out, RH_out, wsp_out, precip_out, press_out, 
     &     Rg_in_out, Rlong_in_out, thisnpd, lst, filltype)

      implicit none
c      include 'fluxtowerdata.h'
      include 'netcdf.inc'

      integer  i, N_out, RCODE, NCID, thisnpd, lst, filltype
	integer nobsmax
	parameter (nobsmax=17520)
      integer thisy, thism, thisf, N_file, starti, dimid(4), dpm(12)
      integer firstyear, nmonths, varids(15), thisnyears
      double precision Tair_out(nobsmax), RH_out(nobsmax)
      double precision wsp_out(nobsmax), precip_out(nobsmax)
      double precision press_out(nobsmax), Rg_in_out(nobsmax)
      double precision Rlong_in_out(nobsmax), dsince(nobsmax)
      double precision temp3d(nobsmax,1,1)
      double precision lat, lon, ht, lat_in(1,1), lon_in(1,1)
      double precision ht_in(nobsmax,1,1)
      character(len=1) filltypest
      character(len=3) mst
      character(len=4) yst
      character(len=11) latst, lonst, htst
      character(len=30) sitename
      character(len=100) fname, cmd

      thisnyears = N_out/(8760*thisnpd/24)
      nmonths = thisnyears*12

      if (lon .lt. 0) lon=lon+360

      starti = 1+lst*thisnpd/24   !output files in GMT
      do thisf = 1,nmonths

      thisy = firstyear+(thisf-1)/12 
      thism = mod(thisf-1,12)+1
     
      DATA dpm/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

C      if (thism .eq. 2 .and. mod(thisy,4) .eq. 0) then 
C         dpm(2)=29
c      else
c         dpm(2)=28
c      end if
         
      write(latst, '(f9.5)')  lat
      write(lonst, '(f11.5)') lon
      write(htst,  '(f5.0)')  ht
      write(yst,   '(I4)')    thisy
      write(mst,   '(I3)')    100+thism

      N_file = dpm(thism)*thisnpd
      !don't include leap year
!      if (thism .eq. 2 .and. mod(thisy,4) .eq. 0) 
!     &     N_file=N_file+thisnpd
      !print*, thisy, thism, starti, N_file, lst, filltype

      cmd = 'mkdir -p input/' // trim(sitename)
      call system(cmd)
      fname= 'input/' // trim(sitename) // '/' //
     &     yst // '-' // mst(2:3) // '.nc'   
      print *, fname

      RCODE=NF_CREATE(trim(fname), NF_CLOBBER, NCID)
      RCODE = NF_DEF_DIM(NCID, 'scalar', 1, dimid(1))
      RCODE = NF_DEF_DIM(NCID, 'lon', 1, dimid(2))
      RCODE = NF_DEF_DIM(NCID, 'lat', 1, dimid(3))     
      RCODE = NF_DEF_DIM(NCID, 'time', N_file, dimid(4))

      !global attributes
      RCODE = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, 'institution', 
     &        29, 'Oak Ridge National Laboratory')
      RCODE = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, 'history', 
     &        99, 'File Origin - This file was created at '
     &     // 'Oak Ridge National Laboratory for the '
     &     // 'site-level MDC in 2010')
      RCODE = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, 'site_location', 
     &        68, 'Latitude: ' // latst // ' Longitude: '
     &        // lonst // ' Elevation (masl): ' // htst)        

      !time
      RCODE = NF_DEF_VAR(NCID, 'time', NF_DOUBLE, 1, dimid(4),
     &                   varids(1))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(1), 'long_name', 16, 
     &                        'observation time')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(1), 'units', 30, 
     &     'days since ' // yst // '-' // mst(2:3) // 
     &     '-01 00:00:00')
c      RCODE = NF_PUT_ATT_TEXT(NCID, varids(1), 'calendar', 9,
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(1), 'calendar', 6, 
c     &                        'gregorian')
     &                        'noleap')
      !longitude
      RCODE = NF_DEF_VAR(NCID, 'LONGXY', NF_DOUBLE, 2, dimid(2:3),
     &     varids(2))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(2), 'long_name', 9,
     &     'longitude')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(2), 'units', 9, 
     &     'degrees E')
      !latitude
      RCODE = NF_DEF_VAR(NCID, 'LATIXY', NF_DOUBLE, 2, dimid(2:3),
     &     varids(3))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(3), 'long_name', 8,
     &     'latitude')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(3), 'units', 9, 
     &     'degrees N')
      !level
      RCODE = NF_DEF_VAR(NCID, 'ZBOT', NF_DOUBLE, 3, dimid(2:4),
     &     varids(4))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(4), 'long_name', 20,
     &     'observational height')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(4), 'units', 1, 
     &     'm')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(4), 'Mode', 14, 
     &     'time-dependent')
      !EDGES
      RCODE=NF_DEF_VAR(NCID, 'EDGEW', NF_DOUBLE, 1, dimid(1), 
     &     varids(5))
      RCODE=NF_PUT_ATT_TEXT(NCID, varids(5), 'long_name', 32, 
     &     'western edge in atmospheric data')
      RCODE=NF_PUT_ATT_TEXT(NCID, varids(5), 'units', 9, 
     &     'degrees E')      
      RCODE=NF_DEF_VAR(NCID, 'EDGEE', NF_DOUBLE, 1, dimid(1), 
     &     varids(6))
      RCODE=NF_PUT_ATT_TEXT(NCID, varids(6), 'long_name', 32, 
     &     'eastern edge in atmospheric data')
      RCODE=NF_PUT_ATT_TEXT(NCID, varids(6), 'units', 9, 
     &     'degrees E')
      RCODE=NF_DEF_VAR(NCID, 'EDGES', NF_DOUBLE, 1, dimid(1), 
     &     varids(7))
      RCODE=NF_PUT_ATT_TEXT(NCID, varids(7), 'long_name', 33, 
     &     'southern edge in atmospheric data')
      RCODE=NF_PUT_ATT_TEXT(NCID, varids(7), 'units', 9, 
     &     'degrees N')
      RCODE=NF_DEF_VAR(NCID, 'EDGEN', NF_DOUBLE, 1, dimid(1), 
     &     varids(8))
      RCODE=NF_PUT_ATT_TEXT(NCID, varids(8), 'long_name', 33, 
     &     'northern edge in atmospheric data')
      RCODE=NF_PUT_ATT_TEXT(NCID, varids(8), 'units', 9, 
     &     'degrees N')

      !Meteorology
      !air temperature
      RCODE = NF_DEF_VAR(NCID, 'TBOT', NF_DOUBLE, 3, dimid(2:4), 
     &     varids(9))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(9), 'long_name', 42,
     &     'temperature at the lowest atm level (TBOT)')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(9), 'units', 1, 'K')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(9), 'mode', 14, 
     &     'time-dependent')
      !relative humidity
      RCODE = NF_DEF_VAR(NCID, 'RH', NF_DOUBLE, 3, dimid(2:4), 
     &     varids(10))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(10), 'long_name', 46,
     &     'relative humidity at the lowest atm level (RH)')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(10), 'units', 1, '%')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(10), 'mode', 14, 
     &     'time-dependent')
      !Wind
      RCODE = NF_DEF_VAR(NCID, 'WIND', NF_DOUBLE, 3, dimid(2:4), 
     &     varids(11))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(11), 'long_name', 35,
     &     'wind at the lowest atm level (WIND)')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(11), 'units', 3, 'm/s')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(11), 'mode', 14, 
     &     'time-dependent')
      !Shortwave radiation
      RCODE = NF_DEF_VAR(NCID, 'FSDS', NF_DOUBLE, 3, dimid(2:4), 
     &     varids(12))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(12), 'long_name', 21,
     &     'incident solar (FSDS)')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(12), 'units', 4, 'W/m2')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(12), 'mode', 14, 
     &     'time-dependent')
      !Longwave radiation
      RCODE = NF_DEF_VAR(NCID, 'FLDS', NF_DOUBLE, 3, dimid(2:4), 
     &     varids(13))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(13), 'long_name', 24,
     &     'incident longwave (FLDS)')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(13), 'units', 4, 'W/m2')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(13), 'mode', 14, 
     &     'time-dependent')
      !Surface pressure
      RCODE = NF_DEF_VAR(NCID, 'PSRF', NF_DOUBLE, 3, dimid(2:4), 
     &     varids(14))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(14), 'long_name', 39,
     &     'pressure at the lowest atm level (PSRF)')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(14), 'units', 2, 'Pa')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(14), 'mode', 14, 
     &     'time-dependent')
      !Precipitation
      RCODE = NF_DEF_VAR(NCID, 'PRECTmms', NF_DOUBLE, 3, dimid(2:4), 
     &     varids(15))
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(15), 'long_name', 24,
     &     'precipitation (PRECTmms)')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(15), 'units', 4, 'mm/s')
      RCODE = NF_PUT_ATT_TEXT(NCID, varids(15), 'mode', 14, 
     &     'time-dependent')

      RCODE=NF_ENDDEF(NCID)

      !create dsince
      do i=1,N_file
         dsince(i) = (i-1)*(1.0/thisnpd)
      end do
      !print*, dsince(1:1000)
      
      !put varibles
      RCODE = NF_PUT_VAR_DOUBLE(NCID, varids(1), dsince(1:N_file))
      lon_in(1,1)=lon
      lat_in(1,1)=lat
      ht_in(:,1,1)=ht
      RCODE = NF_PUT_VAR_DOUBLE(NCID, varids(2), lon_in)
      RCODE = NF_PUT_VAR_DOUBLE(NCID, varids(3), lat_in)
      RCODE = NF_PUT_VAR_DOUBLE(NCID, varids(4), ht_in(1:N_file,1,1))
      RCODE = NF_PUT_VAR_DOUBLE(NCID, varids(5), lon-0.1)
      RCODE = NF_PUT_VAR_DOUBLE(NCID, varids(6), lon+0.1)      
      RCODE = NF_PUT_VAR_DOUBLE(NCID, varids(7), lat-0.1)
      RCODE = NF_PUT_VAR_DOUBLE(NCID, varids(8), lat+0.1)
      
      !Tair (C to K)
      temp3d(1:N_file,1,1)=Tair_out(starti:starti+N_file-1)+273.15
      if (starti .lt. 1) temp3d(1:1-starti,1,1)=Tair_out(1)+273.15
      if (starti+N_file-1 .gt. nobsmax) then
        print *,thisf,starti+N_file-1,N_file-lst*thisnpd/24+1,N_file
        temp3d(N_file-lst*thisnpd/24+1:N_file,1,1)
     &  =Tair_out(nobsmax)+273.15
        endif
      RCODE=NF_PUT_VAR_DOUBLE(NCID, varids(9), temp3d(1:N_file,1,1))
      !RH
      temp3d(1:N_file,1,1)=RH_out(starti:starti+N_file-1)
      if (starti .lt. 1) temp3d(1:1-starti,1,1)=RH_out(1)
      if (starti+N_file-1 .gt. nobsmax)
     &   temp3d(N_file-lst*thisnpd/24+1:N_file,1,1)=
     &  RH_out(nobsmax)
      RCODE=NF_PUT_VAR_DOUBLE(NCID, varids(10), temp3d(1:N_file,1,1))
      !Wind
      temp3d(1:N_file,1,1)=wsp_out(starti:starti+N_file-1)
      if (starti .lt. 1) temp3d(1:1-starti,1,1)=wsp_out(1)
      if (starti+N_file-1 .gt. nobsmax)
     &   temp3d(N_file-lst*thisnpd/24+1:N_file,1,1)=wsp_out(nobsmax)
      RCODE=NF_PUT_VAR_DOUBLE(NCID, varids(11), temp3d(1:N_file,1,1))
      !SWdown
      temp3d(1:N_file,1,1)=Rg_in_out(starti:starti+N_file-1)
      if (starti .lt. 1) temp3d(1:1-starti,1,1)=Rg_in_out(1)
      if (starti+N_file-1 .gt. nobsmax)
     &   temp3d(N_file-lst*thisnpd/24+1:N_file,1,1)=
     &  Rg_in_out(nobsmax)
      RCODE=NF_PUT_VAR_DOUBLE(NCID, varids(12), temp3d(1:N_file,1,1))
      !LWdown
      temp3d(1:N_file,1,1)=Rlong_in_out(starti:starti+N_file-1)
      if (starti .lt. 1) temp3d(1:1-starti,1,1)=Rlong_in_out(1)
      if (starti+N_file-1 .gt. nobsmax)
     &   temp3d(N_file-lst*thisnpd/24+1:N_file,1,1)=
     &  Rlong_in_out(nobsmax)
      RCODE=NF_PUT_VAR_DOUBLE(NCID, varids(13), temp3d(1:N_file,1,1))
      !Psurf
      temp3d(1:N_file,1,1)=press_out(starti:starti+N_file-1)*1000
      if (starti .lt. 1) temp3d(1:1-starti,1,1)=press_out(1)*1000
      if (starti+N_file-1 .gt. nobsmax)
     &   temp3d(N_file-lst*thisnpd/24+1:N_file,1,1)
     &  =press_out(nobsmax)*1000
      RCODE=NF_PUT_VAR_DOUBLE(NCID, varids(14), temp3d(1:N_file,1,1))
      !Precip (mm to mm/s)
      temp3d(1:N_file,1,1)=precip_out(starti:starti+N_file-1)/
     &     (3600*24.0/thisnpd)
      !print*, precip_out(1:2000)
      if (starti .lt. 1) temp3d(1:1-starti,1,1)=precip_out(1)/
     &     (3600*24.0/thisnpd)
      if (starti+N_file-1 .gt. nobsmax)
     &   temp3d(N_file-lst*thisnpd/24+1:N_file,1,1)
     &  =precip_out(nobsmax)/
     &  (3600*24.0/thisnpd)
      RCODE=NF_PUT_VAR_DOUBLE(NCID, varids(15), temp3d(1:N_file,1,1))
      !print*, RCODE
      
  
      RCODE=NF_CLOSE(NCID)
      
      starti=starti+N_file

      !skip Dec. 31 of leap years
      !if (mod(thisy,4) .eq. 0 .and. thism .eq. 12) then 
      !   starti = starti+thisnpd
      !end if

      end do

      end subroutine write_clm
