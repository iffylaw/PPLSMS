#!/usr/bin/env Rscript
# Date: 15/02/2014
# Written by Lihui Luo 
# E-mail: luolh@lzb.ac.cn
# AmeriFlux website http://ameriflux.lbl.gov/Pages/default.aspx
# Get the AmeriFlux site infomation

require(XML)

siteInfoURL <- "http://ameriflux.lbl.gov/AmeriFluxSites/Pages/Site-Map.aspx"
tableID <- "ctl00_ctl33_g_a5f0c0fe_dfe6_435f_8d36_acb194e8c4b8_ctl00_sitesView"
siteInfo <- readHTMLTable(siteInfoURL, header=TRUE, as.data.frame = TRUE)
names(siteInfo)

SITE_ID <- siteInfo[[1]][1]  
SITE_NAME <- siteInfo[[1]][2]  
URL_AMERIFLUX <- siteInfo[[1]][3]
TOWER_BEGAN <- siteInfo[[1]][4]	
TOWER_END <- siteInfo[[1]][5]	
LOCATION_LAT <- siteInfo[[1]][6]	
LOCATION_LONG <- siteInfo[[1]][7]	
LOCATION_ELEV <- siteInfo[[1]][8]	
IGBP <- siteInfo[[1]][9]	
CLIMATE_KOEPPEN <- siteInfo[[1]][10]	
MAT <- siteInfo[[1]][11]	
MAP <- siteInfo[[1]][12]

# write site info table
allSiteInfo <- data.frame(SITE_ID,  SITE_NAME,	TOWER_BEGAN,	TOWER_END,	LOCATION_LAT,	LOCATION_LONG,	LOCATION_ELEV,	IGBP,	CLIMATE_KOEPPEN,	MAT,	MAP)
write.csv(allSiteInfo, file="AmeriFluxSiteInfo.csv")
