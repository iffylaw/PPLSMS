#!/usr/bin/env Rscript
# Date: 09/02/2014
# Written by Tong He (http://weibo.com/u/1635976784) and Lihui Luo
# E-mail: luolh@lzb.ac.cn
# This script generates some animation for NetCDF data files

rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows   

# Set parameters
Args <- commandArgs(trailingOnly=TRUE);
if(length(Args) != 2) {
  message("NetCDFAnimation.R requires file name of NetCDF data and Variable name as input. Terminating");
  quit()
}

require(ncdf)
require(RNetCDF)
require(ggplot2)
require(animation)

NetCDF_name <- as.character(Args[1])
Variable_name <- as.character(Args[2])
work_dir <- try(system("pwd",intern=TRUE))
setwd(work_dir)

# open & read NetCDF file
nc <- open.nc(NetCDF_name)
nc_data <- read.nc(nc)

# get the variable name of lon & lat & time
# sometimes these variable names have different name, such as LONGXY(y, x), LATIXY(y, x), lat(lon, lat), lon(lon, lat), etc
nc_variables <- names(nc_data)
lon <- nc_variables[1]
lat <- nc_variables[2]
time <- nc_variables[3]
time_length <- length(var.get.nc(nc, "time"))

# get the range of variable data, e.g. fivenum function
# min & max value may be singular values
# min_value <- min(unlist(nc_data[[Variable_name]]))
quartile_25 <- quantile(unlist(nc_data[[Variable_name]]))[2]
quartile_25[quartile_25>0] =  quartile_25*0.5
quartile_25[quartile_25<0] =  quartile_25*1.5
quartile_75 <- quantile(unlist(nc_data[[Variable_name]]))[4]
quartile_75[quartile_75>0] =  quartile_75*1.5
quartile_75[quartile_75<0] =  quartile_75*0.5
#max_value <- max(unlist(nc_data[[Variable_name]]))

# time process
date_int <- as.integer(nc_data[[time]])
date_hour <- round((nc_data[[time]]-date_int)*24)

# get the variable unit
Variable_unit <- att.get.nc(nc, Variable_name, "units")
# sometimes "long_name" & "description" are the same attribution 
Variable_long_name <- att.get.nc(nc, Variable_name, "description")

# Draw function, i is the time series
drawit = function(i, nc_data){
    dfr <- data.frame(Lon = as.vector(nc_data[[lon]]), Lat = as.vector(nc_data[[lat]]), vardata = as.vector(nc_data[[Variable_name]][,,i]))
    p <- ggplot(aes(x = Lon, y = Lat, color = vardata), data=dfr)
    p + geom_point() +
      scale_colour_gradientn(
        limits = c(quartile_25, quartile_75),
        colours = c("#0000FF","#00FF00","#FF0000"),
        guide = "colourbar") +
      labs(title = paste(Variable_long_name, "(", Variable_unit, ", ", date_int[i], " ", date_hour[i], "hr.)", sep=''))
}

# set current directory to output the gif animation
ani.options(outdir = getwd())

# Save GIF animation
saveGIF({
    for (i in 1:time_length) print(drawit(i, nc_data))
    }, 
    movie.name=paste(Variable_name,'_Animation.gif',sep=""), interval = 1, ani.width=1000, ani.height=600)

# close the netcdf file
close.nc(nc)
dev.off()