;***********************************************
; color_6.ncl

;***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;************************************************
; Read the file
;************************************************
  ;do Year=2000,2001,1
  Year=2010
  ;files = systemfunc("ls /media/Data/CLM/wrfout-LH/wrfout_LH_"+Year+"*.nc")
  in = addfile("/media/Data/CLM/wrfout-LH/wrfout_LH_201001_monavg_final.nc","r")
  in2 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0.2010-01.nc","r")
  LatentHeat1 = in2->FSH
  printVarSummary(LatentHeat1)
  
  ;f = addfiles(files,"r")
  LatentHeat =LatentHeat1
  LatentHeat = in->LH
  ;XLAT = f[:]->XLAT
  ;XLONG = f[:]->XLONG

  ;LatentHeat!0="time"
  ;LatentHeat!1="XLAT"
  ;LatentHeat!2="XLONG"
  ;LatentHeat&Coordinates = "time lon lat"
  ;LatentHeat&lat@units="degrees_north"
  ;LatentHeat&lon@units="degrees_east"
  ;LatentHeat&latixy@units="degrees_north"
  ;LatentHeat&longxy@units="degrees_east"
  LatentHeat@long_name ="Latent Heat Flux of WRFOUT"
  LatentHeat@units ="W/m^2"
  ;LatentHeat@lat2d = f[:]->XLAT


  ;LatentHeat@lon2d = f[:]->XLONG

  printVarSummary(LatentHeat)

  ;in = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-01.nc","r")
  ;in2 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-02.nc","r")
  ;in3 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-03.nc","r")
  ;in4 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-04.nc","r")
  ;in5 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-05.nc","r")
  ;in6 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-06.nc","r")
  ;in7 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-07.nc","r")
  ;in8 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-08.nc","r")
  ;in9 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-09.nc","r")
  ;in10 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-10.nc","r")
  ;in11 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-11.nc","r")
  ;in12 = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-12.nc","r")

  ;LatentHeat = in->FSH
  ;LatentHeat2 = in2->FSH
  ;LatentHeat3 = in3->FSH
  ;LatentHeat4 = in4->FSH
  ;LatentHeat5 = in5->FSH
  ;LatentHeat6 = in6->FSH
  ;LatentHeat7 = in7->FSH
  ;LatentHeat8 = in8->FSH
  ;LatentHeat9 = in9->FSH
  ;LatentHeat10 = in10->FSH
  ;LatentHeat11 = in11->FSH
  ;LatentHeat12 = in12->FSH

;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("ps","CLM_LatentHeat_"+Year)                 ; open a png file

  gsn_define_colormap(wks, "BlAqGrYeOrReVi200")                ; choose colormap
  plot = new(13,graphic)
 
  res                      = True                  ; plot mods desired
 
  res@cnFillOn             = True                  ; turn on color fill
  res@gsnSpreadColors      = True                  ; use full range of colors
  res@cnLinesOn            = False

  res@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
  res@lbOrientation   = "Vertical" 
  res@lbLabelStride            = 2
  res@gsnAddCyclic = False

  ;res@mpMinLatF = 25.34 
  ;res@mpMaxLatF = 39.72 
  ;res@mpMinLonF = 72.73 
  ;res@mpMaxLonF = 107.27 

  res@cnLevelSelectionMode = "ExplicitLevels" 


  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-01)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(0) = gsn_csm_contour_map_ce(wks,LatentHeat(0,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-02)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(1) = gsn_csm_contour_map_ce(wks,LatentHeat(1,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-03)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(2) = gsn_csm_contour_map_ce(wks,LatentHeat(2,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-04)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(3) = gsn_csm_contour_map_ce(wks,LatentHeat(3,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-05)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(4) = gsn_csm_contour_map_ce(wks,LatentHeat(4,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-06)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(5) = gsn_csm_contour_map_ce(wks,LatentHeat(5,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-07)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(6) = gsn_csm_contour_map_ce(wks,LatentHeat(6,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-08)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(7) = gsn_csm_contour_map_ce(wks,LatentHeat(7,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-09)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(8) = gsn_csm_contour_map_ce(wks,LatentHeat(8,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-10)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(9) = gsn_csm_contour_map_ce(wks,LatentHeat(9,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-11)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(10) = gsn_csm_contour_map_ce(wks,LatentHeat(10,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux of WRFOUT in Qinghai-Tibet Plateau ("+Year+"-12)" 
  res@cnLevels            = (/0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(11) = gsn_csm_contour_map_ce(wks,LatentHeat(11,:,:), res) 
  delete (res@cnLevels)

  
  pres             = True                    ; panel mods desired
  pres@gsnMaximize = True                    ; fill the page
  gsn_panel(wks,plot,(/4,3/),pres)           ; draw: 2-across, 2-down
  ;end do
end
