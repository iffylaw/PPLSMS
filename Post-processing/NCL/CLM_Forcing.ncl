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
  nlat = 150
  mlon = 300
 
  
  in = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-01.nc","r")
  in2 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-02.nc","r")
  in3 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-03.nc","r")
  in4 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-04.nc","r")
  in5 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-05.nc","r")
  in6 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-06.nc","r")
  in7 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-07.nc","r")
  in8 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-08.nc","r")
  in9 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-09.nc","r")
  in10 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-10.nc","r")
  in11 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-11.nc","r")
  in12 = addfile("/media/Data/CLM/QZ_CLM_Forcing/2011-12.nc","r")
  
  TBOT1 = in->TBOT                                             
  ;S_Humi = clmMonTLL(TBOT)    ; Climatology
  S_Humi1 = dim_avg_n(TBOT1, 0)  

  TBOT2 = in2->TBOT                                             
  S_Humi2 = dim_avg_n(TBOT2, 0)  
            
  TBOT3 = in3->TBOT                                             
  S_Humi3 = dim_avg_n(TBOT3, 0)  

  TBOT4 = in4->TBOT                                             
  S_Humi4 = dim_avg_n(TBOT4, 0)  

  TBOT5 = in5->TBOT                                             
  S_Humi5 = dim_avg_n(TBOT5, 0)  

  TBOT6 = in6->TBOT                                             
  S_Humi6 = dim_avg_n(TBOT6, 0)  

  TBOT7 = in7->TBOT                                             
  S_Humi7 = dim_avg_n(TBOT7, 0)  

  TBOT8 = in8->TBOT                                             
  S_Humi8 = dim_avg_n(TBOT8, 0)  

  TBOT9 = in9->TBOT                                             
  S_Humi9 = dim_avg_n(TBOT9, 0)  

  TBOT10 = in10->TBOT                                             
  S_Humi10 = dim_avg_n(TBOT10, 0)  

  TBOT11 = in11->TBOT                                             
  S_Humi11 = dim_avg_n(TBOT11, 0)  

  TBOT12 = in12->TBOT                                             
  S_Humi12 = dim_avg_n(TBOT12, 0)  

  ;S_Humi!0 ="time"
  ;S_Humi!1 ="lat"
  ;S_Humi!2 ="lon"
  ;S_Humi&lat =lat
  ;S_Humi&lon =lon
  ;S_Humi@long_name = " Specific Humidity at 2 M"
  ;S_Humi@units = "kg/kg"

  S_Humi1@long_name = "Temperature at 2 M"
  S_Humi1@units = "K"
  S_Humi1@lat2d = in->LATIXY
  S_Humi1@lon2d = in->LONGXY

  S_Humi2@long_name = "Temperature at 2 M"
  S_Humi2@units = "K"
  S_Humi2@lat2d = in->LATIXY
  S_Humi2@lon2d = in->LONGXY

  S_Humi3@long_name = "Temperature at 2 M"
  S_Humi3@units = "K"
  S_Humi3@lat2d = in->LATIXY
  S_Humi3@lon2d = in->LONGXY

  S_Humi4@long_name = "Temperature at 2 M"
  S_Humi4@units = "K"
  S_Humi4@lat2d = in->LATIXY
  S_Humi4@lon2d = in->LONGXY

  S_Humi5@long_name = "Temperature at 2 M"
  S_Humi5@units = "K"
  S_Humi5@lat2d = in->LATIXY
  S_Humi5@lon2d = in->LONGXY

  S_Humi6@long_name = "Temperature at 2 M"
  S_Humi6@units = "K"
  S_Humi6@lat2d = in->LATIXY
  S_Humi6@lon2d = in->LONGXY

  S_Humi7@long_name = "Temperature at 2 M"
  S_Humi7@units = "K"
  S_Humi7@lat2d = in->LATIXY
  S_Humi7@lon2d = in->LONGXY

  S_Humi8@long_name = "Temperature at 2 M"
  S_Humi8@units = "K"
  S_Humi8@lat2d = in->LATIXY
  S_Humi8@lon2d = in->LONGXY

  S_Humi9@long_name = "Temperature at 2 M"
  S_Humi9@units = "K"
  S_Humi9@lat2d = in->LATIXY
  S_Humi9@lon2d = in->LONGXY

  S_Humi10@long_name = "Temperature at 2 M"
  S_Humi10@units = "K"
  S_Humi10@lat2d = in->LATIXY
  S_Humi10@lon2d = in->LONGXY

  S_Humi11@long_name = "Temperature at 2 M"
  S_Humi11@units = "K"
  S_Humi11@lat2d = in->LATIXY
  S_Humi11@lon2d = in->LONGXY

  S_Humi12@long_name = "Temperature at 2 M"
  S_Humi12@units = "K"
  S_Humi12@lat2d = in->LATIXY
  S_Humi12@lon2d = in->LONGXY

  ;TBOT&lat@units="LATIXY"
  ;TBOT&lon@units="LONGXY"

  ;printVarSummary(S_Humi) 
  ;print(min(S_Humi))
  ;print(max(S_Humi))

  ;print(min(TBOT))
  ;print(max(TBOT))
  ;print(TBOT@time)

 
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("ps","CLM_Forcing_TBOT_2011")
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

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-01)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(0) = gsn_csm_contour_map_ce(wks,S_Humi1(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-02)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(1) = gsn_csm_contour_map_ce(wks,S_Humi2(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-03)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(2) = gsn_csm_contour_map_ce(wks,S_Humi3(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-04)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(3) = gsn_csm_contour_map_ce(wks,S_Humi4(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-05)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(4) = gsn_csm_contour_map_ce(wks,S_Humi5(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-06)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(5) = gsn_csm_contour_map_ce(wks,S_Humi6(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-07)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(6) = gsn_csm_contour_map_ce(wks,S_Humi7(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-08)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(7) = gsn_csm_contour_map_ce(wks,S_Humi8(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-09)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(8) = gsn_csm_contour_map_ce(wks,S_Humi9(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-10)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(9) = gsn_csm_contour_map_ce(wks,S_Humi10(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-11)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(10) = gsn_csm_contour_map_ce(wks,S_Humi11(:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Temperature at 2 M in Qinghai-Tibet Plateau (2011-12)"     ; title
  res@cnLevels            = (/225,230,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310/)
  plot(11) = gsn_csm_contour_map_ce(wks,S_Humi12(:,:), res) 
  delete (res@cnLevels)
  
  pres             = True                    ; panel mods desired
  pres@gsnMaximize = True                    ; fill the page
  gsn_panel(wks,plot,(/4,3/),pres)           ; draw: 2-across, 2-down

end
