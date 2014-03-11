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
  ;nlat = 150
  ;mlon = 300
 
  do Year=2000,2010,1
  
  in = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-01.nc","r")
  in2 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-02.nc","r")
  in3 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-03.nc","r")
  in4 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-04.nc","r")
  in5 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-05.nc","r")
  in6 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-06.nc","r")
  in7 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-07.nc","r")
  in8 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-08.nc","r")
  in9 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-09.nc","r")
  in10 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-10.nc","r")
  in11 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-11.nc","r")
  in12 = addfile("/media/Data/CLM/QZ_CLM_Forcing/"+Year+"-12.nc","r")
  
  FSDS1 = in->FSDS                                             
  ;S_Humi = clmMonTLL(FSDS)    ; Climatology
  S_Humi1 = dim_avg_n(FSDS1, 0)  
  delete(FSDS1)

  FSDS2 = in2->FSDS                                             
  S_Humi2 = dim_avg_n(FSDS2, 0) 
  delete(FSDS2) 
            
  FSDS3 = in3->FSDS                                             
  S_Humi3 = dim_avg_n(FSDS3, 0)  
  delete(FSDS3)

  FSDS4 = in4->FSDS                                             
  S_Humi4 = dim_avg_n(FSDS4, 0)  
  delete(FSDS4)

  FSDS5 = in5->FSDS                                             
  S_Humi5 = dim_avg_n(FSDS5, 0)  
  delete(FSDS5)

  FSDS6 = in6->FSDS                                             
  S_Humi6 = dim_avg_n(FSDS6, 0)  
  delete(FSDS6)

  FSDS7 = in7->FSDS                                             
  S_Humi7 = dim_avg_n(FSDS7, 0)  
  delete(FSDS7)

  FSDS8 = in8->FSDS                                             
  S_Humi8 = dim_avg_n(FSDS8, 0)  
  delete(FSDS8)

  FSDS9 = in9->FSDS                                             
  S_Humi9 = dim_avg_n(FSDS9, 0)  
  delete(FSDS9)

  FSDS10 = in10->FSDS                                             
  S_Humi10 = dim_avg_n(FSDS10, 0)  
  delete(FSDS10)

  FSDS11 = in11->FSDS                                             
  S_Humi11 = dim_avg_n(FSDS11, 0)  
  delete(FSDS11)

  FSDS12 = in12->FSDS                                             
  S_Humi12 = dim_avg_n(FSDS12, 0)  
  delete(FSDS12)

  ;S_Humi!0 ="time"
  ;S_Humi!1 ="lat"
  ;S_Humi!2 ="lon"
  ;S_Humi&lat =lat
  ;S_Humi&lon =lon
  ;S_Humi@long_name = " Specific Humidity at 2 M"
  ;S_Humi@units = "kg/kg"

  S_Humi1@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi1@units = "W/m^2"
  S_Humi1@lat2d = in->LATIXY
  S_Humi1@lon2d = in->LONGXY

  S_Humi2@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi2@units = "W/m^2"
  S_Humi2@lat2d = in->LATIXY
  S_Humi2@lon2d = in->LONGXY

  S_Humi3@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi3@units = "W/m^2"
  S_Humi3@lat2d = in->LATIXY
  S_Humi3@lon2d = in->LONGXY

  S_Humi4@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi4@units = "W/m^2"
  S_Humi4@lat2d = in->LATIXY
  S_Humi4@lon2d = in->LONGXY

  S_Humi5@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi5@units = "W/m^2"
  S_Humi5@lat2d = in->LATIXY
  S_Humi5@lon2d = in->LONGXY

  S_Humi6@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi6@units = "W/m^2"
  S_Humi6@lat2d = in->LATIXY
  S_Humi6@lon2d = in->LONGXY

  S_Humi7@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi7@units = "W/m^2"
  S_Humi7@lat2d = in->LATIXY
  S_Humi7@lon2d = in->LONGXY

  S_Humi8@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi8@units = "W/m^2"
  S_Humi8@lat2d = in->LATIXY
  S_Humi8@lon2d = in->LONGXY

  S_Humi9@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi9@units = "W/m^2"
  S_Humi9@lat2d = in->LATIXY
  S_Humi9@lon2d = in->LONGXY

  S_Humi10@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi10@units = "W/m^2"
  S_Humi10@lat2d = in->LATIXY
  S_Humi10@lon2d = in->LONGXY

  S_Humi11@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi11@units = "W/m^2"
  S_Humi11@lat2d = in->LATIXY
  S_Humi11@lon2d = in->LONGXY

  S_Humi12@long_name = "Downward Short Wave Flux at Ground Surface"
  S_Humi12@units = "W/m^2"
  S_Humi12@lat2d = in->LATIXY
  S_Humi12@lon2d = in->LONGXY

  ;FSDS&lat@units="LATIXY"
  ;FSDS&lon@units="LONGXY"

  ;printVarSummary(S_Humi) 
  ;print(min(S_Humi))
  ;print(max(S_Humi))

  ;print(min(FSDS))
  ;print(max(FSDS))
  ;print(FSDS@time)

 
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("ps","CLM_Forcing_FSDS_"+Year)
  ;colors = (/"white","black","royal blue","light sky blue","powder blue","light sea green","pale green","wheat","brown","pink"/)

  gsn_define_colormap(wks, "BlueRed")                ; choose colormap
  plot = new(12,graphic)
 
  res                      = True                  ; plot mods desired
 
  res@cnFillOn             = True                  ; turn on color fill
  res@gsnSpreadColors      = True                  ; use full range of colors
  res@cnLinesOn            = False


  res@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels

  res@lbOrientation   = "Vertical" 
  res@lbLabelStride   = 2
  res@gsnAddCyclic    = False

  res@mpMinLatF = 25 
  res@mpMaxLatF = 40 
  res@mpMinLonF = 72 
  res@mpMaxLonF = 108 

  res@cnLevelSelectionMode = "ExplicitLevels" 

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-01)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(0) = gsn_csm_contour_map_ce(wks,S_Humi1(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-02)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(1) = gsn_csm_contour_map_ce(wks,S_Humi2(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-03)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(2) = gsn_csm_contour_map_ce(wks,S_Humi3(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-04)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(3) = gsn_csm_contour_map_ce(wks,S_Humi4(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-05)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(4) = gsn_csm_contour_map_ce(wks,S_Humi5(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-06)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(5) = gsn_csm_contour_map_ce(wks,S_Humi6(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-07)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(6) = gsn_csm_contour_map_ce(wks,S_Humi7(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-08)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(7) = gsn_csm_contour_map_ce(wks,S_Humi8(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-09)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(8) = gsn_csm_contour_map_ce(wks,S_Humi9(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-10)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(9) = gsn_csm_contour_map_ce(wks,S_Humi10(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-11)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(10) = gsn_csm_contour_map_ce(wks,S_Humi11(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Downward Short Wave Flux at Ground Surface in Qinghai-Tibet Plateau ("+Year+"-12)"     ; title
  res@cnLevels            = (/0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450/)
  plot(11) = gsn_csm_contour_map_ce(wks,S_Humi12(:,:), res) 
  delete (res@cnLevels)
  
  pres             = True                    ; panel mods desired
  pres@gsnMaximize = True                    ; fill the page
  gsn_panel(wks,plot,(/4,3/),pres)           ; draw: 2-across, 2-down
  
  end do
end
