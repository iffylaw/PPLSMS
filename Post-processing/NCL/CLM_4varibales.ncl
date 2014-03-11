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
 
  Year="2000"
  Month="01"
  in = addfile("/media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"-"+Month+".nc","r")
  
  
  SoilWater = in->H2OSOI
  FCTR = in->FCTR
  FCEV = in->FCEV
  FGEV = in->FGEV

  LatentHeat = FCTR + FCEV + FGEV
  LatentHeat@long_name = "Latent Heat Flux"
  LatentHeat@units = "W/m^2"
  LatentHeat@lat3d = in->lat
  LatentHeat@lon3d = in->lon

  SensibleHeat = in->FSH
  GroundTemp  = in->TG
  
  ;tg_c=tg-273.15
  ;tg_c@long_name ="Ground Temperature"
  ;tg_c@units ="Centi Degree"
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("ps","CLM_Foring_4vars_"+Year+"-"+Month)                 ; open a png file
   colors = (/"white","black","white","royal blue","light sky blue",\
	     "powder blue","light sea green","pale green","wheat","brown",\
	     "pink"/)

  gsn_define_colormap(wks, "gui_default")                ; choose colormap
  plot = new(5,graphic)
 
  res                      = True                  ; plot mods desired
 
  res@cnFillOn             = True                  ; turn on color fill
  res@gsnSpreadColors      = True                  ; use full range of colors
  res@cnLinesOn            = False

; when using a blue/red color map, it is best to set the contour min/max
; to equal but opposite values so that the lightest colors are centered 
; on zero
  res@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
 ;res@cnMinLevelValF       =  min(TS)                ; set min contour level
; res@cnMaxLevelValF       =  max(TS)                  ; set max contour level
 ; res@cnLevelSpacingF      =   10                  ; set contour spacing
res@lbOrientation   = "Vertical" 
res@lbLabelStride            = 2
res@gsnAddCyclic = False
;res@mpDataSetName         = "Earth..4"
;res@mpDataBaseVersion     = "MediumRes"  ; Medium resolution database

  res@mpMinLatF = 25.34 
  res@mpMaxLatF = 39.72 
  res@mpMinLonF = 72.73 
  res@mpMaxLonF = 107.27 

;res@mpGeophysicalLineThicknessF= 2.      ; double the thickness of geophysical boundaries
 ; res@mpNationalLineThicknessF= 2.         ; double the thickness of national boundaries
 ; delete(res@mpOutlineSpecifiers)          ; about to change size of attribute, must delete
 ; res@mpOutlineSpecifiers = (/"China:states","Taiwan"/)
 
  res@cnLevelSelectionMode = "ExplicitLevels" 

;res@cnFillColors    = False
;print("pot:  min="+ min(TS) + "   max="+ max(TS) );

  res@tiMainString         = "Sensible Heat Flux in Qinghai-Tibet Plateau ("+Year+"-"+Month+")" 
  res@cnLevels            = (/-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170/)
  plot(0) = gsn_csm_contour_map_ce(wks,SensibleHeat(0,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Latent Heat Flux in Qinghai-Tibet Plateau ("+Year+"-"+Month+")"     ; title
  res@cnLevels            = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30/)
  plot(1) = gsn_csm_contour_map_ce(wks,LatentHeat(0,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "Ground Temperature in Qinghai-Tibet Plateau ("+Year+"-"+Month+")"     ; title
  res@cnLevels            = (/240,245,250,255,260,265,270,275,280,285,290/)
  plot(2) = gsn_csm_contour_map_ce(wks,GroundTemp(0,:,:), res) 
  delete (res@cnLevels)

  res@tiMainString         = "First Layer Soil Water Content for Qinghai-Tibet Plateau ("+Year+"-"+Month+")"     ; title
  res@cnLevels            = (/0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8/)
  plot(3) = gsn_csm_contour_map_ce(wks,SoilWater(0,0,:,:), res) 
  delete (res@cnLevels)
  
  pres             = True                    ; panel mods desired
  pres@gsnMaximize = True                    ; fill the page
  gsn_panel(wks,plot,(/2,2/),pres)           ; draw: 2-across, 2-down

end
