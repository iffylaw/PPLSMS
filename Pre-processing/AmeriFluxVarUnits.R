#!/usr/bin/env Rscript
# Date: 15/02/2014
# Written by Lihui Luo 
# E-mail: luolh@lzb.ac.cn
# FLUXNET Variable Units website http://www.fluxdata.org/DataInfo/default.aspx
# LEVEL 4 VARIABLE DESCRIPTION (Half hourly dataset variables description)
# website ftp://cdiac.ornl.gov/pub/ameriflux/data/Level4/L4_README.txt ()
# Get the AmeriFlux variable units

AmeriVarInfo <- read.csv("AmeriFluxVar.csv", header=TRUE)

AmeriVarName <- AmeriVarInfo[,1]
AmeriVarUnit <- AmeriVarInfo[,3]
AmeriVarDesc <- AmeriVarInfo[,2]

# write site info table
AmeriVarUnitInfo <- data.frame(AmeriVarName, AmeriVarUnit, AmeriVarDesc)
