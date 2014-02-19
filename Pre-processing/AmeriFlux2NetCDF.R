#!/usr/bin/env Rscript
# Date: 09/02/2014
# Written by Lihui Luo 
# E-mail: luolh@lzb.ac.cn
# AmeriFlux website http://ameriflux.lbl.gov/Pages/default.aspx
# Data product come from Level 4 - Gap-filled & Adjusted Data Files with GEP & Re Estimates
# Convert AmeriFlux data to NetCDF format files

rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows   

# Set parameters
Args <- commandArgs(trailingOnly=TRUE);
if(length(Args) != 1) {
  message("AmeriFlux2NetCDF.R requires site name as input. Terminating");
  quit()
}

require(XML)
require(ncdf)
require(RNetCDF)

site_name <- as.character(Args[1])
work_dir <- try(system("pwd",intern=TRUE))
setwd(work_dir)

# include the script that get the AmeriFlux site infomation
source("AmeriFluxSiteInfo.R")
# include the script that get the variable units of AmeriFlux site
source("AmeriFluxVarUnits.R")

start_year <- as.numeric(siteInfo[[1]][4][siteInfo[[1]][1]==site_name])  
end_year <- as.numeric(siteInfo[[1]][5][siteInfo[[1]][1]==site_name])
longitude <- as.numeric(siteInfo[[1]][7][siteInfo[[1]][1]==site_name])
latitude <- as.numeric(siteInfo[[1]][8][siteInfo[[1]][1]==site_name])

# combine the multiple files of hourly with one header
listFiles <- list.files(pattern=".*_h_.*\\.txt$", recursive=TRUE)
listFilesData <- do.call("rbind", lapply(listFiles, read.csv, header = TRUE))

site_filename <- paste(site_name,"-", start_year,"-", end_year, ".csv", sep="")
write.csv(listFilesData, file=site_filename)

# Read FLUXNET file
nc = read.csv(site_filename, header=TRUE, sep=",")
attach(nc)

Variable_names <- names(nc)
Variable_length <- length(nc[,1])

# Create NetCDF file
netcdf.from.fluxnet <- create.nc(paste(site_filename, ".nc", sep=""))

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
var.put.nc(netcdf.from.fluxnet, "time", as.numeric(as.character(Hour[1:Variable_length])), 1, Variable_length)

for (V in 6:length(Variable_names)){
  
  #define dimensions
  var.def.nc(netcdf.from.fluxnet, Variable_names[V], "NC_DOUBLE", c("lat","lon","time"))
  
  #define long name of variables
  varDesc <- levels(droplevels(AmeriVarDesc[AmeriVarName==Variable_names[V]]))
  varDesc[is.na(varDesc)] <- as.character("--")
  att.put.nc(netcdf.from.fluxnet, Variable_names[V], "long_name", "NC_CHAR", varDesc)
  
  #define unit of variables
  varUnit <- levels(droplevels(AmeriVarUnit[AmeriVarName==Variable_names[V]]))
  varUnit[is.na(varUnit)] <- as.character("--")
  att.put.nc(netcdf.from.fluxnet, Variable_names[V], "units", "NC_CHAR", varUnit)
  
  # define missing value
  att.put.nc(netcdf.from.fluxnet, Variable_names[V], "missing_value", "NC_DOUBLE", -9999.)
  
  # Write data out to NetCDF file
  var.put.nc(netcdf.from.fluxnet, Variable_names[V], as.numeric(as.character(nc[,V][1:Variable_length])), start=c(1,1,1), count=c(1,1,Variable_length))
  
}

# Global attribution
att.put.nc(netcdf.from.fluxnet, "NC_GLOBAL", "title", "NC_CHAR", "Data from FLUXNET")
att.put.nc(netcdf.from.fluxnet, "NC_GLOBAL", "author", "NC_CHAR", "Lihui Luo email:luolh@lzb.ac.cn")

# close the netcdf file
sync.nc(netcdf.from.fluxnet)
close.nc(netcdf.from.fluxnet)
