#!/usr/bin/env Rscript
# Date: 21/02/2014
# Written by Lihui Luo 
# E-mail: luolh@lzb.ac.cn
# FLUXNET gap_filled_marconi data website http://daac.ornl.gov//FLUXNET/guides/marconi_gap_filled.html
# Get the site infomation of FLUXNET gap_filled_marconi data 

require(XML)

siteInfoURL <- "http://daac.ornl.gov//FLUXNET/guides/marconi_gap_filled.html"
siteInfo <- readHTMLTable(siteInfoURL, header=TRUE, as.data.frame = TRUE)

siteCodeInfo <- siteInfo[[1]]
FlxHalfhVar <- siteInfo[[5]]
MeteHalfhVar <- siteInfo[[10]]