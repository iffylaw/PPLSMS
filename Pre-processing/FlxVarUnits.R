#!/usr/bin/env Rscript
# Date: 15/02/2014
# Written by Lihui Luo 
# E-mail: luolh@lzb.ac.cn
# FLUXNET Variable Units website http://www.fluxdata.org/DataInfo/default.aspx
# In Flux/met variables in the 1/2 hourly files (download) -- NewHourly.xls
# Get the FLUXNET variable units

VarInfo <- read.csv("NewHourly.csv", header=TRUE)

VarName <- VarInfo[,2]
VarUnit <- VarInfo[,3]
VarDesc <- VarInfo[,4]
VarCubeDatum <- VarInfo[,5]
VarCubeExtDatum <- VarInfo[,6]
VarCubeOffset <- VarInfo[,7]
