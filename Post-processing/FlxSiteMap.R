#!/usr/bin/env Rscript
# Date: 26/02/2014
# Written Lihui Luo
# E-mail: luolh@lzb.ac.cn
# FLUXNET site network data http://bwc.berkeley.edu/StaticReports/Fluxnet/SitesByNetwork.xls
# Making Maps for sites

rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows   

require(maps)
require(ggplot2)

work_dir <- try(system("pwd",intern=TRUE))
setwd(work_dir)

# Fluxnet sites infomation
sites <- read.csv("SitesByNetwork.csv", header=TRUE)
site_lon <- sites$Longitude
site_lat <- sites$Latitude
network <- sites$Network

# world map
world <- map_data("world")

# plot the fluxnet sites map
ggplot() +
  geom_polygon( data=world, aes(x=long, y=lat, group = group), colour="grey10", fill="white" ) +
  geom_point(aes(site_lon, site_lat, colour = network), size=2, pch = 20)

# save as png
ggsave(file="FlxSiteMap.png")

# close the file
dev.off()