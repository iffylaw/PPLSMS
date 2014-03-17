#!/usr/bin/env Rscript
library(ncdf)
library(gdata)
library(sp)
library(fields)
library(gplots)

site_surface_dir <- "/public1/home/luolh/FLUXNET/surface"
all_sites <- read.table("all_sites_list.txt", header=TRUE, as.is = TRUE)
attach(all_sites)
lonlat <- array(0,c(2,length(SITE_NAME)))
for (S in 1:length(SITE_NAME)) {
   site_nc <- open.ncdf(paste(site_surface_dir, "/", SITE_NAME[S], ".", NTILES[S], "_tiles.surface.nc", sep=""))
   print(SITE_NAME[S])
   lonlat[1,S] <- get.var.ncdf(site_nc, "lon")
   lonlat[2,S] <- get.var.ncdf(site_nc, "lat")[1]
}

LON <- lonlat[1,]

LON[LON>=180] = LON[LON>=180] - 360

print(LON)
LAT <- lonlat[2,]

all_sites_list <- data.frame(SITE_NAME, START_YEAR, END_YEAR, NTILES, TIMEOFFSET, LON, LAT)
print(all_sites_list)
write.table(all_sites_list, file= "all_sites.txt", col.names= names(all_sites_list))

# Plot lon&lat in world map
opar<-par( cex=0.6, pin=c(6.0, 3.0))
load("TM_WORLD_BORDERS_SIMPL-0.2.RData")
plot(wrld_simpl, axes = TRUE, xlim = c(-180,180), ylim = c(-120, 90), asp = 1)
grid()
for (L in 1:length(LON)){
points(LON[L],LAT[L], col=rainbow(length(LON))[L], pch=substr(SITE_NAME[L],1,1))
}
legend('bottom', SITE_NAME, pch=substr(SITE_NAME,1,1), col=1:length(SITE_NAME), bty="n", ncol=8, cex=0.5)
