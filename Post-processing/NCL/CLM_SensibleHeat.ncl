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
  do Year=2007,2010,1

  files = systemfunc("ls /media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"*.nc")
  f = addfiles(files,"r")
  SensibleHeat = f[:]->FSH
  ;print(SensibleHeat(0,:,:))

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

  ;SensibleHeat = in->FSH
  ;SensibleHeat2 = in2->FSH
  ;SensibleHeat3 = in3->FSH
  ;SensibleHeat4 = in4->FSH
  ;SensibleHeat5 = in5->FSH
  ;SensibleHeat6 = in6->FSH
  ;SensibleHeat7 = in7->FSH
  ;SensibleHeat8 = in8->FSH
  ;SensibleHeat9 = in9->FSH
  ;SensibleHeat10 = in10->FSH
  ;SensibleHeat11 = in11->FSH
  ;SensibleHeat12 = in12->FSH

;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("ps","CLM_SensibleHeat_"+Year)                 ; open a png file

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

  res@mpMinLatF = 25.34 
  res@mpMaxLatF = 39.72 
  res@mpMinLonF = 72.73 
  res@mpMaxLonF = 107.27 

  res@cnLevelSelectionMode = "ExplicitLevels" 


  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-01)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(0) = gsn_csm_contour_map_ce(wks,SensibleHeat(0,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-02)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(1) = gsn_csm_contour_map_ce(wks,SensibleHeat(1,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-03)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(2) = gsn_csm_contour_map_ce(wks,SensibleHeat(2,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-04)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(3) = gsn_csm_contour_map_ce(wks,SensibleHeat(3,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-05)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(4) = gsn_csm_contour_map_ce(wks,SensibleHeat(4,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-06)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(5) = gsn_csm_contour_map_ce(wks,SensibleHeat(5,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-07)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(6) = gsn_csm_contour_map_ce(wks,SensibleHeat(6,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-08)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(7) = gsn_csm_contour_map_ce(wks,SensibleHeat(7,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-09)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(8) = gsn_csm_contour_map_ce(wks,SensibleHeat(8,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-10)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(9) = gsn_csm_contour_map_ce(wks,SensibleHeat(9,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-11)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(10) = gsn_csm_contour_map_ce(wks,SensibleHeat(10,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-12)" 
  res@cnLevels            = (/-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150/)
  plot(11) = gsn_csm_contour_map_ce(wks,SensibleHeat(11,:,:), res) 
  delete (res@cnLevels)

  
  pres             = True                    ; panel mods desired
  pres@gsnMaximize = True                    ; fill the page
  gsn_panel(wks,plot,(/4,3/),pres)           ; draw: 2-across, 2-down
  end do
end
