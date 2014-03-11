;***********************************************
; color_6.ncl

;***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin
;************************************************
; Read the file
;************************************************
 
  
  in    = addfile("/public1/home/root/QZ_2008/QZ_2008.clm2.h0.2008-03.nc","r")
    
  
  fsh    = in->FSH
  fpsn = in->FPSN
  qbot =in->QBOT
  tg = in->TG
  tg_c=(tg-32)/1.8
  tg_c@long_name ="ground temperature"
  tg_c@units ="Centi Degree"
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("png","china-2004")                 ; open a png file
   colors = (/"white","black","white","royal blue","light sky blue",\
	     "powder blue","light sea green","pale green","wheat","brown",\
	     "pink"/)

  gsn_define_colormap(wks, "gui_default")                ; choose colormap
  plot = new(4,graphic)
 
  res                      = True                  ; plot mods desired
  res@tiMainString         = "China (2004)"     ; title
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

  res@mpMinLatF = 25. 
  res@mpMaxLatF = 39. 
  res@mpMinLonF = 72. 
  res@mpMaxLonF = 107. 

;res@mpGeophysicalLineThicknessF= 2.      ; double the thickness of geophysical boundaries
 ; res@mpNationalLineThicknessF= 2.         ; double the thickness of national boundaries
 ; delete(res@mpOutlineSpecifiers)          ; about to change size of attribute, must delete
 ; res@mpOutlineSpecifiers = (/"China:states","Taiwan"/)
 
  res@cnLevelSelectionMode = "ExplicitLevels" 
 ; res@cnLevels            = (/-160,-120,-80,-40,0,40,80,120,160,200,240,280,320/)

;res@cnFillColors    = False
;print("pot:  min="+ min(TS) + "   max="+ max(TS) );

  res@cnLevels            = (/-240,-200,-160,-120,-80,-40,0,40,80,120,160,200,240,280,320,360,400,440,480,520/)
  plot(0) = gsn_csm_contour_map_ce(wks,fsh(0,:,:), res)   
  delete (res@cnLevels) 
  res@cnLevels            = (/1,3,5,7,9,11,13,15,17,19/)
  plot(1) = gsn_csm_contour_map_ce(wks,fpsn(0,:,:), res)  
  delete (res@cnLevels) 
  res@cnLevels            = (/0.002,0.006,0.01,0.014,0.018,0.022,0.026,0.03,0.034,0.038,0.042,0.046/)
  plot(2) = gsn_csm_contour_map_ce(wks,qbot(0,:,:), res)
  delete (res@cnLevels)  
  res@cnLevels            = (/200,240,280,320,360,400,440,480/)
  plot(3) = gsn_csm_contour_map_ce(wks,tg(0,:,:), res) 

  pres             = True                    ; panel mods desired
  pres@gsnMaximize = True                    ; fill the page
  gsn_panel(wks,plot,(/2,2/),pres)           ; draw: 2-across, 2-down

end
