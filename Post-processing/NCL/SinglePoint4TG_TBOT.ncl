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
do Year=2000,2004,1

  files = systemfunc("ls /media/Data/CLM/CLMRUN-QZ-MONTHLY/QZ_2008.clm2.h0."+Year+"*.nc")
  in = addfiles(files,"r")
       
  Group_Temp = in[:] ->TG
  Air_Temp = in[:] ->TBOT
 
      
print("Air Temp = "+Air_Temp(0,65,236) +" Group Temp = " +Group_Temp(0,65,236))
print("Air Temp = "+Air_Temp(1,65,236)+" Group Temp = " +Group_Temp(1,65,236))
print("Air Temp = "+Air_Temp(2,65,236) +" Group Temp = " +Group_Temp(2,65,236))
print("Air Temp = "+Air_Temp(3,65,236)+" Group Temp = " +Group_Temp(3,65,236))
print("Air Temp = "+Air_Temp(4,65,236)+" Group Temp = " +Group_Temp(4,65,236))
print("Air Temp = "+Air_Temp(5,65,236) +" Group Temp = " +Group_Temp(5,65,236))
print("Air Temp = "+Air_Temp(6,65,236) +" Group Temp = " +Group_Temp(6,65,236))
print("Air Temp = "+Air_Temp(7,65,236) +" Group Temp = " +Group_Temp(7,65,236))
print("Air Temp = "+Air_Temp(8,65,236) +" Group Temp = " +Group_Temp(8,65,236))
print("Air Temp = "+Air_Temp(9,65,236) +" Group Temp = " +Group_Temp(9,65,236))
print("Air Temp = "+Air_Temp(10,65,236) +" Group Temp = " +Group_Temp(10,65,236))
print("Air Temp = "+Air_Temp(11,65,236) +" Group Temp = " +Group_Temp(11,65,236))
 
end do
end

