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
  
  WIND1 = in->WIND                                             
  ;S_Humi = clmMonTLL(WIND)    ; Climatology
  S_Humi1 = dim_avg_n(WIND1, 0)  
  delete(WIND1)

  WIND2 = in2->WIND                                             
  S_Humi2 = dim_avg_n(WIND2, 0) 
  delete(WIND2) 
            
  WIND3 = in3->WIND                                             
  S_Humi3 = dim_avg_n(WIND3, 0)  
  delete(WIND3)

  WIND4 = in4->WIND                                             
  S_Humi4 = dim_avg_n(WIND4, 0)  
  delete(WIND4)

  WIND5 = in5->WIND                                             
  S_Humi5 = dim_avg_n(WIND5, 0)  
  delete(WIND5)

  WIND6 = in6->WIND                                             
  S_Humi6 = dim_avg_n(WIND6, 0)  
  delete(WIND6)

  WIND7 = in7->WIND                                             
  S_Humi7 = dim_avg_n(WIND7, 0)  
  delete(WIND7)

  WIND8 = in8->WIND                                             
  S_Humi8 = dim_avg_n(WIND8, 0)  
  delete(WIND8)

  WIND9 = in9->WIND                                             
  S_Humi9 = dim_avg_n(WIND9, 0)  
  delete(WIND9)

  WIND10 = in10->WIND                                             
  S_Humi10 = dim_avg_n(WIND10, 0)  
  delete(WIND10)

  WIND11 = in11->WIND                                             
  S_Humi11 = dim_avg_n(WIND11, 0)  
  delete(WIND11)

  WIND12 = in12->WIND                                             
  S_Humi12 = dim_avg_n(WIND12, 0)  
  delete(WIND12)

  ;S_Humi!0 ="time"
  ;S_Humi!1 ="lat"
  ;S_Humi!2 ="lon"
  ;S_Humi&lat =lat
  ;S_Humi&lon =lon
  ;S_Humi@long_name = " Specific Humidity at 2 M"
  ;S_Humi@units = "kg/kg"

  S_Humi1@long_name = "Wind at the lowest atm level"
  S_Humi1@units = "m/s"
  S_Humi1@lat2d = in->LATIXY
  S_Humi1@lon2d = in->LONGXY

  S_Humi2@long_name = "Wind at the lowest atm level"
  S_Humi2@units = "m/s"
  S_Humi2@lat2d = in->LATIXY
  S_Humi2@lon2d = in->LONGXY

  S_Humi3@long_name = "Wind at the lowest atm level"
  S_Humi3@units = "m/s"
  S_Humi3@lat2d = in->LATIXY
  S_Humi3@lon2d = in->LONGXY

  S_Humi4@long_name = "Wind at the lowest atm level"
  S_Humi4@units = "m/s"
  S_Humi4@lat2d = in->LATIXY
  S_Humi4@lon2d = in->LONGXY

  S_Humi5@long_name = "Wind at the lowest atm level"
  S_Humi5@units = "m/s"
  S_Humi5@lat2d = in->LATIXY
  S_Humi5@lon2d = in->LONGXY

  S_Humi6@long_name = "Wind at the lowest atm level"
  S_Humi6@units = "m/s"
  S_Humi6@lat2d = in->LATIXY
  S_Humi6@lon2d = in->LONGXY

  S_Humi7@long_name = "Wind at the lowest atm level"
  S_Humi7@units = "m/s"
  S_Humi7@lat2d = in->LATIXY
  S_Humi7@lon2d = in->LONGXY

  S_Humi8@long_name = "Wind at the lowest atm level"
  S_Humi8@units = "m/s"
  S_Humi8@lat2d = in->LATIXY
  S_Humi8@lon2d = in->LONGXY

  S_Humi9@long_name = "Wind at the lowest atm level"
  S_Humi9@units = "m/s"
  S_Humi9@lat2d = in->LATIXY
  S_Humi9@lon2d = in->LONGXY

  S_Humi10@long_name = "Wind at the lowest atm level"
  S_Humi10@units = "m/s"
  S_Humi10@lat2d = in->LATIXY
  S_Humi10@lon2d = in->LONGXY

  S_Humi11@long_name = "Wind at the lowest atm level"
  S_Humi11@units = "m/s"
  S_Humi11@lat2d = in->LATIXY
  S_Humi11@lon2d = in->LONGXY

  S_Humi12@long_name = "Wind at the lowest atm level"
  S_Humi12@units = "m/s"
  S_Humi12@lat2d = in->LATIXY
  S_Humi12@lon2d = in->LONGXY

  ;WIND&lat@units="LATIXY"
  ;WIND&lon@units="LONGXY"

  ;printVarSummary(S_Humi) 
  ;print(min(S_Humi))
  ;print(max(S_Humi))

  ;print(min(WIND))
  ;print(max(WIND))
  ;print(WIND@time)

 
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("ps","CLM_Forcing_WIND_"+Year)
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

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-01)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(0) = gsn_csm_contour_map_ce(wks,S_Humi1(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-02)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(1) = gsn_csm_contour_map_ce(wks,S_Humi2(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-03)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(2) = gsn_csm_contour_map_ce(wks,S_Humi3(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-04)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(3) = gsn_csm_contour_map_ce(wks,S_Humi4(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-05)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(4) = gsn_csm_contour_map_ce(wks,S_Humi5(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-06)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(5) = gsn_csm_contour_map_ce(wks,S_Humi6(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-07)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(6) = gsn_csm_contour_map_ce(wks,S_Humi7(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-08)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(7) = gsn_csm_contour_map_ce(wks,S_Humi8(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-09)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(8) = gsn_csm_contour_map_ce(wks,S_Humi9(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-10)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(9) = gsn_csm_contour_map_ce(wks,S_Humi10(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-11)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(10) = gsn_csm_contour_map_ce(wks,S_Humi11(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Wind at the lowest atm level in Qinghai-Tibet Plateau ("+Year+"-12)"     ; title
  res@cnLevels            = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  plot(11) = gsn_csm_contour_map_ce(wks,S_Humi12(:,:), res) 
  delete (res@cnLevels)
  
  pres             = True                    ; panel mods desired
  pres@gsnMaximize = True                    ; fill the page
  gsn_panel(wks,plot,(/4,3/),pres)           ; draw: 2-across, 2-down
  
  end do
end
