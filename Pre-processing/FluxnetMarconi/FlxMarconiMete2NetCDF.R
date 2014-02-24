#!/usr/bin/env Rscript
# Date: 01/02/2014
# Written by Lihui Luo 
# E-mail: luolh@lzb.ac.cn
# download website FLUXNET data from ftp://daac.ornl.gov/data/fluxnet/gap_filled_marconi/data/Meteo/Halfhourly
# Convert FLUXNET Meteo Halfhourly data to NetCDF format files

rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows   

# Set parameters
Args <- commandArgs(trailingOnly=TRUE);
if(length(Args) != 3) {
  message("FlxMete2NetCDF.R requires site code (2 or 3 letter), lon & lat as input, e.g. AB97_hh.met is for AB. Terminating");
  quit()
}

# Input parameter
siteCode <- as.character(Args[1])
latitude <- as.numeric(Args[2])
longitude <- as.numeric(Args[3])
work_dir <- try(system("pwd",intern=TRUE))

setwd(work_dir)

require(ncdf)
require(RNetCDF)

# include the script that get the FLUXNET gap_filled_marconi data
source("FlxMarconiSiteInfo.R")

start_year <- as.numeric(substr(siteCodeInfo[6][siteCodeInfo[1]==siteCode], 1, 4))
end_year <- as.numeric(substr(siteCodeInfo[6][siteCodeInfo[1]==siteCode], 6, 9))
end_year[is.na(end_year)] <- start_year

#Time processing from start year %Y%m%d.%f, 20010101.0000000
# Days as time interval, one day has 48 halfhour
Obs_date <- strftime(seq(as.Date(paste(start_year,"-01-01",sep="")),  as.Date(paste(end_year,"-12-31",sep="")), by = "days"), "%Y%m%d")
Obs_time <- seq(0,1,length=49)[1:48]
Obs_date_time <- array(NA, dim=c(48, length(Obs_date)))
for (D in 1:length(Obs_date)){
  Obs_date_time[1:48, D] <- paste(Obs_date[D],".", substr(Obs_time, 3, 10)[1:length(Obs_time)], sep="")
}
Obs_date_time <- as.vector(Obs_date_time)

# combine the multiple files of hourly with one header
listFiles <- list.files(pattern=paste(siteCode, ".*\\.met$", sep=""), recursive=TRUE)
if (length(listFiles)>1){
  listFilesData <- do.call("rbind", lapply(listFiles, read.table, header = TRUE, skip = 1)) 
} else {
  listFilesData <- read.table(listFiles, header=TRUE, skip = 1)
}

site_filename <- paste(siteCode, "-", start_year,"-", end_year, ".csv", sep="")
write.table(listFilesData, file=site_filename, sep=",")

# Read FLUXNET file
nc = read.csv(site_filename, header=TRUE)
#attach(nc)

Variable_names <- MeteHalfhVar$Variable
# Variable_unit <- nc[1,]
# Convert data.frame columns from factors to characters
#Variable_units <- data.frame(lapply(Variable_unit, as.character), stringsAsFactors=FALSE) 
Variable_length <- length(nc[,1])

# Create NetCDF file
netcdf.from.fluxnet <- create.nc(paste(siteCode,"-", start_year,"-", end_year, "_obs_halfhourly.nc", sep=""))

# dimensions
dim.def.nc(netcdf.from.fluxnet, "lon", 1)
dim.def.nc(netcdf.from.fluxnet, "lat", 1)
dim.def.nc(netcdf.from.fluxnet, "time", unlim=TRUE)

# define lon & lat
var.def.nc(netcdf.from.fluxnet, "lon", "NC_DOUBLE", "lon")
var.def.nc(netcdf.from.fluxnet, "lat", "NC_DOUBLE", "lat")
att.put.nc(netcdf.from.fluxnet, "lon", "long_name", "NC_CHAR", "longitude")
att.put.nc(netcdf.from.fluxnet, "lon", "units", "NC_CHAR", "degrees_noth")
att.put.nc(netcdf.from.fluxnet, "lat", "long_name", "NC_CHAR", "latitude")
att.put.nc(netcdf.from.fluxnet, "lat", "units", "NC_CHAR", "degrees_east")

# loop for variables define (var & att)
# First & Sencond Variables are Day & Hour, time varibale will replace
var.def.nc(netcdf.from.fluxnet, "time", "NC_DOUBLE", "time")
att.put.nc(netcdf.from.fluxnet, "time", "long_name", "NC_FLOAT", "time")
att.put.nc(netcdf.from.fluxnet, "time", "units", "NC_CHAR", "day as %Y%m%d.%f")

# define lat & lon
var.put.nc(netcdf.from.fluxnet, "lon", longitude)
var.put.nc(netcdf.from.fluxnet, "lat", latitude)
var.put.nc(netcdf.from.fluxnet, "time", as.numeric(as.character(Obs_date_time)), 1, Variable_length)

for (V in 3:length(Variable_names)){
    
  #define dimensions
  var.def.nc(netcdf.from.fluxnet, levels(droplevels(Variable_names[V])), "NC_DOUBLE", c("lat","lon","time"))
  
  #define long name of variables
  varDesc <- levels(droplevels(MeteHalfhVar$Description[MeteHalfhVar$Variable==levels(droplevels(Variable_names[V]))]))
  varDesc[is.na(varDesc)] <- as.character("--")
  att.put.nc(netcdf.from.fluxnet, levels(droplevels(Variable_names[V])), "long_name", "NC_CHAR", varDesc)
  
  #define unit of variables
  varUnit <- levels(droplevels(MeteHalfhVar$Units[MeteHalfhVar$Variable==levels(droplevels(Variable_names[V]))]))
  varUnit[is.na(varUnit)] <- as.character("--")
  att.put.nc(netcdf.from.fluxnet, levels(droplevels(Variable_names[V])), "units", "NC_CHAR", varUnit)
    
  # define missing value
  att.put.nc(netcdf.from.fluxnet, levels(droplevels(Variable_names[V])), "missing_value", "NC_DOUBLE", -9999.)
    
  # Write data out to NetCDF file
  var.put.nc(netcdf.from.fluxnet, levels(droplevels(Variable_names[V])), as.numeric(as.character(nc[,V][1:Variable_length])), start=c(1,1,1), count=c(1,1,Variable_length))
    
}

# Global attribution
att.put.nc(netcdf.from.fluxnet, "NC_GLOBAL", "title", "NC_CHAR", "Data from FLUXNET gap_filled_marconi data")
att.put.nc(netcdf.from.fluxnet, "NC_GLOBAL", "author", "NC_CHAR", "Lihui Luo email:luolh@lzb.ac.cn")

# close the netcdf file
sync.nc(netcdf.from.fluxnet)
close.nc(netcdf.from.fluxnet)

