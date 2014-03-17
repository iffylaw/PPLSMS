#!/usr/bin/env Rscript
# Date 06/07/2011
# This script generates some txt file from statistics file 
# written by Lihui Luo
# contains:
# ConvertSite2Pars - function for convert each site's stats into each parameter's stats 

#---------------------------------------------------------------------------------------------------
rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows   

out_dir <- "/public1/home/luolh/CLM4.0/Output"
version <- "LSMS0.1"
FORCING_TSTEP <- "halfhourly"
work_dir <- try(system("pwd",intern=TRUE))

setwd(work_dir)

# define required compared variables, can be added later.
nvar=9
variable_names <- array(c("par_acc", "PPFD", "PPFD", "PPFD",
                          "net_radiation", "Rn", "Net Radiation", "Rn",
                          "sensible_heat_flx", "H", "Sensible Heat Flux", "Qh" , 
                          "latent_heat_flx", "LE", "Latent Heat Flux", "Qle",  
                          "canopy_cond_limited", "gsurf", "Canopy Conductance", "gscan" ,
                          "net_assimilation", "GPP", "GPP", "GPP",
                          "reco", "Reco", "Reco", "Reco",
                          "net_co2_flux", "NEE", "NEE", "NEE",
                          "fapar", "fapar", "FAPAR", "fapar",
                          "soil_moisture", "SWC2", "Soil Water Content", "SWC2")
                        , dim=c(4,nvar))
variable_units <- array(c("mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "W/m^2", "W/m^2", "W/m^2", 
                          "W/m^2", "W/m^2", "W/m^2", 
                          "W/m^2", "W/m^2", "W/m^2", 
                          "m/s","mmol m-2 s-1","m/s",
                          "mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "-","-","-",
                          "-","m","-")
                        , dim=c(3,nvar))
variable_units_conversion <- array(c(10^6,1,
                                     1,1,
                                     -1,1,
                                     -1,1,
                                      1,18*10^-6,
                                     10^6,1,
                                     10^6,1,
                                     10^6,1,
                                     1,1,
				     1,1)
                        , dim=c(2,nvar))

site_info <- read.table("all_sites_list.txt",header=TRUE,as.is = TRUE)
len_site <- length(site_info$SITE_NAME)
ab_var <- substr(variable_names[3,],1,1)
VAR_list <- variable_names[4,]
ab_site <- substr(site_info$SITE_NAME,1,1)

STAT_list <- c("intercept","slope", "r", "rr", "rmse", "nrmse", "sd_jsb", "sd_flx", "mef", "nae", "vr", "pbias", "nse", "rsr" )
summary <- c("JSBACH", "Observed", "Units")
 
# Convert statistical file of every site to statistical file of every parameter 
period <- c("halfhourly", "daily", "monthly")
for (I in 1:length(period)){
  for (VAR in 1:length(VAR_list)) {
    write.table(t(c("SITE", STAT_list)), file=paste(VAR_list[VAR], ".", version, ".",  FORCING_TSTEP, ".stats.", period[I], ".txt", sep=""), 
      col.names=FALSE, row.names=FALSE)
  } 
}

write.table(t(c("SITE", summary)), file=paste(VAR_list[VAR], ".", version, ".",  FORCING_TSTEP, ".summary.txt", sep=""),
      col.names=FALSE, row.names=FALSE)

# Function: Convert Site to Parameter
Convert_Site_to_Parameter <- function(){
  
  for (site in 1:len_site){
    site_data_halfhourly <- read.table(paste(out_dir, "/", site_info$SITE[site], "/Output/", site_info$SITE[site], ".", version, ".", 
      FORCING_TSTEP, ".stats.", period[1], ".txt", sep=""))
    site_data_daily <- read.table(paste(out_dir, "/", site_info$SITE[site], "/Output/", site_info$SITE[site], ".", version, ".",
      FORCING_TSTEP, ".stats.", period[2], ".txt", sep=""))
    site_data_monthly <- read.table(paste(out_dir, "/", site_info$SITE[site], "/Output/", site_info$SITE[site], ".", version, ".", 
      FORCING_TSTEP, ".stats.", period[3], ".txt", sep=""))
    site_data_summary <- read.table(paste(out_dir, "/", site_info$SITE[site], "/Output/", site_info$SITE[site], ".", version, ".", 
      FORCING_TSTEP, ".summary.txt", sep=""))
	
    for (VAR in 1:length(VAR_list)) {
      write.table(t(c(site_info$SITE[site], t(site_data_halfhourly)[((VAR-1)*length(STAT_list)+1):(VAR*length(STAT_list))])),
        file=paste(VAR_list[VAR], ".", version, ".",  FORCING_TSTEP, ".stats.halfhourly.txt", sep=""), append = TRUE, col.names=FALSE, row.names=FALSE)
      write.table(t(c(site_info$SITE[site], t(site_data_daily)[((VAR-1)*length(STAT_list)+1):(VAR*length(STAT_list))])), 
        file=paste(VAR_list[VAR], ".", version, ".",  FORCING_TSTEP, ".stats.daily.txt", sep=""), append = TRUE, col.names=FALSE, row.names=FALSE)
      write.table(t(c(site_info$SITE[site], t(site_data_monthly)[((VAR-1)*length(STAT_list)+1):(VAR*length(STAT_list))])), 
        file=paste(VAR_list[VAR], ".", version, ".",  FORCING_TSTEP, ".stats.monthly.txt", sep=""), append = TRUE, col.names=FALSE, row.names=FALSE)
      write.table(t(c(site_info$SITE[site], t(site_data_summary)[((VAR-1)*3+1):(VAR*3)])),
        file=paste(VAR_list[VAR], ".", version, ".",  FORCING_TSTEP, ".summary.txt", sep=""), append = TRUE, col.names=FALSE, row.names=FALSE)
    }
  }
}

Convert_Site_to_Parameter()
# The optimum value in each stats
optimum_value <- function(stat){
  if(stat=="pbias"|stat=="rsr"|stat=="intercept"|stat=="rmse"|stat=="nrmse"|stat=="mef"|stat=="nae")  h=0
  if(stat=="nse"|stat=="slope"|stat=="r"|stat=="rr"|stat=="vr")       h=1
  abline(h=h, col= "#0000ff22", cex=0.6) 
}

# For statistcs file in stat directory, if variables is not null, just need to plot the variable
plot_all_stat <- function (time, stat, variable=NA){
  opar<-par( cex=0.6, pin=c(6.0, 3.0))
  s_pos <- which(STAT_list==stat)
  p_pos <- which(period==time)
  stat_v <- NA 
  site <- site_info$SITE_NAME

  if (is.na(variable)){
    for (V in 1:nvar){
      stat_v[V] <- read.table(paste(VAR_list[V], ".", version, ".",  FORCING_TSTEP, ".stats.", time, ".txt", sep=""), header=TRUE, as.is = TRUE)[s_pos+1]
    }
    stat_u <- unlist(stat_v)
    stat_u[stat_u==-Inf|stat_u==Inf] <- NA
    
    ymin<-min(stat_u, na.rm=TRUE)
    ymax<-max(stat_u, na.rm=TRUE)
 
    for (V in 1:nvar){    
      if (V==1){
        plot(1:length(site), t(unlist(stat_v[V])),  ylim=c(ymin,ymax), col=V, pch=V, xaxt = "n",  xlab="", 
          ylab=toupper(stat), font.lab=2, type="b", las=1, main=toupper(paste(time, stat)), cex=0.5)
      }
      if (V>=2){
        lines(1:length(site), t(unlist(stat_v[V])), ylim=c(ymin,ymax), col=V, pch=V, type="b", cex=0.5)
      }
      if (V==nvar){    
        optimum_value(stat)
        grid(length(site), 5, lwd = 0.4)
        legend("bottomright", legend = VAR_list, col=1:nvar, pch=1:nvar, ncol=3, bty="n") 
        axis(1,1:length(site), labels = site, las=2)
      }
    }
  }
  else{
    v_pos <- which(VAR_list==variable)
    stat_v <- read.table(paste(VAR_list[v_pos], ".", version, ".",  FORCING_TSTEP, ".stats.", time, ".txt", sep=""), header=TRUE, as.is = TRUE)[s_pos+1]
    stat_v[stat_v==-Inf|stat_v==Inf]=NA
    ymin<-min(stat_v, na.rm=TRUE)
    ymax<-max(stat_v, na.rm=TRUE)
    plot(1:length(site), t(stat_v),  ylim=c(ymin,ymax), col="red", pch=1,  xaxt = "n", xlab="",
          ylab=toupper(stat), font.lab=2, type="b", las=1, main=toupper(paste(time, stat, "in", variable)), cex=0.5)
    optimum_value(stat)
    grid(length(site), 5, lwd = 0.4)
    axis(1,1:length(site), labels = site, las=2)
  }
    print(warnings())
}

plot_all_stat("halfhourly", "pbias")
plot_all_stat("halfhourly", "nse")
plot_all_stat("halfhourly", "rsr")
plot_all_stat("daily", "pbias")
plot_all_stat("daily", "nse")
plot_all_stat("daily", "rsr")
plot_all_stat("monthly", "pbias")
plot_all_stat("monthly", "nse")
plot_all_stat("monthly", "rsr")
plot_all_stat("daily", "pbias", variable="Qh")
plot_all_stat("daily", "pbias", variable="PPFD")
plot_all_stat("daily", "pbias", variable="Rn")
plot_all_stat("daily", "pbias", variable="GPP")
plot_all_stat("daily", "pbias", variable="NEE")
plot_all_stat("daily", "pbias", variable="gscan")
plot_all_stat("daily", "pbias", variable="Reco")
plot_all_stat("daily", "pbias", variable="Qle")
plot_all_stat("daily", "pbias", variable="fapar")
