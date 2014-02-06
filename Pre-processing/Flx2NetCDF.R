#!/usr/bin/env Rscript
# Date: 01/02/2014
# Written by Lihui Luo luolh@lzb.ac.cn

rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows   

# Set parameters
Args <- commandArgs(trailingOnly=TRUE);
if(length(Args) != 3) {
  message("Flx2NetCDF.R requires data file name,  latitude and longitude as input. Terminating");
  quit()
}

# Input parameter
ncfile <- as.character(Args[1])
latitude <- as.numeric(Args[2])
longitude <- as.numeric(Args[3])
work_dir <- try(system("pwd",intern=TRUE))

setwd(work_dir)

library(ncdf)
library(RNetCDF)
  
# Read FLUXNET file
nc = read.csv(ncfile, header=TRUE, sep="")
attach(nc)

Variable_names <- names(nc)
Variable_unit <- nc[1,]
# Convert data.frame columns from factors to characters
Variable_units <- data.frame(lapply(Variable_unit, as.character), stringsAsFactors=FALSE) 
Variable_length <- length(nc[,1])

# Create NetCDF file
netcdf.from.fluxnet <- create.nc(paste(ncfile, ".nc", sep=""))

# dimensions
dim.def.nc(netcdf.from.fluxnet, "lon", 1)
dim.def.nc(netcdf.from.fluxnet, "lat", 1)
dim.def.nc(netcdf.from.fluxnet, "time", unlim=TRUE)

# variables
var.def.nc(netcdf.from.fluxnet, "lon", "NC_DOUBLE", "lon")
var.def.nc(netcdf.from.fluxnet, "lat", "NC_DOUBLE", "lat")
var.def.nc(netcdf.from.fluxnet, "time", "NC_DOUBLE", "time")
  
# loop for variables define (var & att)
# First & Sencond Variables are Day & Hour, time varibale will replace
att.put.nc(netcdf.from.fluxnet, "time", "long_name", "NC_FLOAT", "time")
att.put.nc(netcdf.from.fluxnet, "time", "units", "NC_FLOAT", "hour")

# define lat & lon
var.put.nc(netcdf.from.fluxnet, "lon", longitude)
var.put.nc(netcdf.from.fluxnet, "lat", latitude)
var.put.nc(netcdf.from.fluxnet, "time", as.numeric(as.character(Hour[2:Variable_length])), 1, Variable_length-1)

for (V in 3:length(Variable_names)){
    
    #dim.def.nc(netcdf.from.fluxnet, Variable_names[V], Variable_length-1)
    var.def.nc(netcdf.from.fluxnet, Variable_names[V], "NC_DOUBLE", c("lat","lon","time"))
    
    #att.put.nc(netcdf.from.fluxnet, "NEE", "long_name", "NC_CHAR", "gapfilled Net Ecosystem Exchange")
    att.put.nc(netcdf.from.fluxnet, Variable_names[V], "units", "NC_CHAR", Variable_units[1, V])
    
    att.put.nc(netcdf.from.fluxnet, Variable_names[V], "missing_value", "NC_DOUBLE", -9999.)
    
    # Write data out to NetCDF file
    var.put.nc(netcdf.from.fluxnet, Variable_names[V], as.numeric(as.character(nc[,V][2:Variable_length])), start=c(1,1,1), count=c(1,1,Variable_length-1))
    
}

# Global attribution
att.put.nc(netcdf.from.fluxnet, "NC_GLOBAL", "title", "NC_CHAR", "Data from FLUXNET")
att.put.nc(netcdf.from.fluxnet, "NC_GLOBAL", "author", "NC_CHAR", "Lihui Luo email:luolh@lzb.ac.cn")

# close the netcdf file
sync.nc(netcdf.from.fluxnet)
close.nc(netcdf.from.fluxnet)

