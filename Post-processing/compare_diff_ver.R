#!/usr/bin/env Rscript
# Date 30/05/2011
# This script generates some simple diagnostics and plot for different version simulations
# written by Lihui Luo and Soenke Zaehle
# contains:
# earlier_current_GRAPHIC - function for ordinary current plots
# earlier_current_GRAPHIC_ERR - function for error bars 
# earlier_current_COMPUTE_GROUPING - function to aggregate to daily, monthly, annual and diurnal cycle
# earlier_current_EVALUATION - function to compute Root Mean Square Error, R sqaure, NSE, PBIAS, RSR
# 
# main 

# List of things to do
# account for quality flags in observations when aggregating / plotting
# strage behaviour of first plot (black squares

#---------------------------------------------------------------------------------------------------
rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows   

# Graphic Function (1) current STUFF-----------------------------------------------------------------
earlier_current_GRAPHIC <- function(earlier, current, FLUXNET_V, TIME, typestr, gtype, labs, variable_info, site, R.1, R.2, R.3, RMSE.1, RMSE.2, RMSE.3){

opar<-par( cex=0.6, pin=c(6.0, 2.0))
split.screen(c(2,1))
split.screen(c(1,3),screen=2)
screen(1)
  plot(TIME, earlier, type="n",  col="red", 
    ylim=c(min(earlier, current, FLUXNET_V, na.rm=TRUE), max(earlier, current, FLUXNET_V, na.rm=TRUE)), 
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
    main=paste(typestr, " ", variable_info[1], " at ", site, " (", variable_info[3], "-", variable_info[4], ")", sep=""),
    las=1, axes=F)
  axis(1, at=1:length(TIME), lab=labs)
  axis(2)
  box()  
  lines(TIME, earlier, type=gtype, col="green", pch=21)
  lines(TIME, current, type=gtype, col="red", pch=21)
  lines(TIME, FLUXNET_V, type=gtype, col="black", pch=21)
  legend("topleft", legend = c(earlier_diff,current_diff,"FLUXNET"),  col=c("green","red","black"), pch=c(21,21,21), bty="n", ncol=3)

screen(3)
  plot(FLUXNET_V, earlier, xlab='FLUXNET', ylab=earlier_diff, las=1,
      ylim=c(min(earlier, FLUXNET_V,na.rm=TRUE), max(earlier, FLUXNET_V,na.rm=TRUE)),
      xlim=c(min(earlier, FLUXNET_V,na.rm=TRUE), max(earlier, FLUXNET_V,na.rm=TRUE)),
      main= paste(variable_info[1], " at ", site, "(" ,  variable_info[2], ")", sep= "" ))
  line_FJ <- line(FLUXNET_V,earlier) 
  intercept=NA
  slope=NA
  if(!all(is.na(coef(line_FJ)))&&!all(is.nan(coef(line_FJ)))&&!all(is.infinite(coef(line_FJ)))) {
    intercept=sprintf("%.2f",coef(line_FJ,use.na=false)[1])
    slope=sprintf("%.2f",coef(line_FJ,use.na=false)[2])
	abline(line_FJ, col="red")
  }
  abline(0,1)
  legend("bottomright", legend = c(paste("slope = ", slope), paste("intercept = ", intercept), paste("R = ", R.1), paste("RMSE = ", RMSE.1)), bty="n")
  
screen(4)
  plot(FLUXNET_V, current, xlab='FLUXNET',ylab=current_diff,las=1,
      ylim=c(min(FLUXNET_V, current,na.rm=TRUE), max(FLUXNET_V, current,na.rm=TRUE)),
      xlim=c(min(FLUXNET_V, current,na.rm=TRUE), max(FLUXNET_V, current,na.rm=TRUE)),
      main= paste(variable_info[1], " at ", site, "(" ,  variable_info[2], ")", sep= "" ))
  line_FJ <- line(FLUXNET_V,current) 
  intercept=NA
  slope=NA
  if(!all(is.na(coef(line_FJ)))&&!all(is.nan(coef(line_FJ)))&&!all(is.infinite(coef(line_FJ)))) {
    intercept=sprintf("%.2f",coef(line_FJ,use.na=false)[1])
    slope=sprintf("%.2f",coef(line_FJ,use.na=false)[2])
    abline(line_FJ, col="red")
  }
  abline(0,1)
  legend("bottomright", legend = c(paste("slope = ", slope), paste("intercept = ", intercept), paste("R = ", R.2), paste("RMSE = ", RMSE.2)), bty="n")
  
screen(5)
  plot(current,earlier, xlab=current_diff,ylab=earlier_diff,las=1,
      ylim=c(min(earlier, current,na.rm=TRUE), max(earlier, current,na.rm=TRUE)),
      xlim=c(min(earlier, current,na.rm=TRUE), max(earlier, current,na.rm=TRUE)),
      main= paste(variable_info[1], " at ", site, "(" ,  variable_info[2], ")", sep= "" ))
  line_FJ <- line(current,earlier) 
  intercept=NA
  slope=NA
  if(!all(is.na(coef(line_FJ)))&&!all(is.nan(coef(line_FJ)))&&!all(is.infinite(coef(line_FJ)))) {
    intercept=sprintf("%.2f",coef(line_FJ,use.na=false)[1])
    slope=sprintf("%.2f",coef(line_FJ,use.na=false)[2])
	abline(line_FJ, col="red")
  }
  abline(0,1)
  legend("bottomright", legend = c(paste("slope = ", slope), paste("intercept = ", intercept), paste("R = ", R.3), paste("RMSE = ", RMSE.3)), bty="n")
  
close.screen(all = TRUE) 
} 
  
# Graphic Function (2) SATELLITE STUFF-----------------------------------------------------------------
earlier_current_GRAPHIC_FAPAR <- function(earlier, current, FLUXNET_V, FLUXNET_V2, FLUXNET_V3, TIME, typestr, gtype, variable_info, site, R_seawifs.1, R_seawifs.2, R_seawifs.3, RMSE_seawifs.1, RMSE_seawifs.2, RMSE_seawifs.3){

opar<-par( cex=0.6, pin=c(6.0, 2.0))
split.screen(c(2,1))
split.screen(c(1,3),screen=2)
screen(1)
  plot(TIME, earlier, type=gtype,  col="red", pch=19,
    ylim=c(min(earlier, current, c(FLUXNET_V,FLUXNET_V2,FLUXNET_V3),na.rm=TRUE), 
           max(earlier, current, c(FLUXNET_V,FLUXNET_V2,FLUXNET_V3),na.rm=TRUE)), 
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
    main=paste(typestr, " ", variable_info[1], " at ", site, " (", variable_info[3],"-",variable_info[4], ")", sep=""),
    las=1)
  lines(TIME, earlier, type=gtype, col="green", pch=21)
  lines(TIME, current, type=gtype, col="red", pch=21)
  lines(TIME, FLUXNET_V, type=gtype, col="black" )
  lines(TIME, FLUXNET_V2, type=gtype, col="blue" )
  lines(TIME, FLUXNET_V3, type=gtype, col="grey" )
  legend("bottomright", legend = c(earlier_diff,current_diff,"SeaWifs","Modis","Cyclopes"),  
    col=c("green","red","black","blue","grey"), pch=c(21,21,21,21,21), bty="n", ncol=5)

screen(3)
  plot(FLUXNET_V,earlier, xlab='SeaWifs',ylab=earlier_diff,las=1,xlim=c(0,1),ylim=c(0,1),
      main= paste(variable_info[1], " at ", site, "(" ,  variable_info[2], ")", sep= "" ))
  line_FJ <- line(FLUXNET_V,earlier) 
  intercept=NA
  slope=NA
  if(!all(is.na(coef(line_FJ)))&&!all(is.nan(coef(line_FJ)))&&!all(is.infinite(coef(line_FJ)))) {
    intercept=sprintf("%.2f",coef(line_FJ,use.na=false)[1])
    slope=sprintf("%.2f",coef(line_FJ,use.na=false)[2])
    abline(line_FJ, col="red")
  }
  abline(0,1)
  legend("bottomright", legend = c(paste("slope = ", slope), paste("intercept = ", intercept), paste("R = ", R_seawifs.1), paste("RMSE = ", RMSE_seawifs.1)), bty="n")
  
screen(4)
  plot(FLUXNET_V,current, xlab='SeaWifs',ylab=current_diff,las=1,xlim=c(0,1),ylim=c(0,1),
      main= paste(variable_info[1], " at ", site, "(" ,  variable_info[2], ")", sep= "" ))
  line_FJ <- line(FLUXNET_V,current) 
  intercept=NA
  slope=NA
  if(!all(is.na(coef(line_FJ)))&&!all(is.nan(coef(line_FJ)))&&!all(is.infinite(coef(line_FJ)))) {
    intercept=sprintf("%.2f",coef(line_FJ,use.na=false)[1])
    slope=sprintf("%.2f",coef(line_FJ,use.na=false)[2])
    abline(line_FJ, col="red")
  }
  abline(0,1)
  legend("bottomright", legend = c(paste("slope = ", slope), paste("intercept = ", intercept), paste("R = ", R_seawifs.2), paste("RMSE = ", RMSE_seawifs.2)), bty="n")
  
screen(5)
  plot(current,earlier, xlab=current_diff,ylab=earlier_diff,las=1,xlim=c(0,1),ylim=c(0,1),
      main= paste(variable_info[1], " at ", site, "(" ,  variable_info[2], ")", sep= "" ))
  line_FJ <- line(current,earlier) 
  intercept=NA
  slope=NA
  if(!all(is.na(coef(line_FJ)))&&!all(is.nan(coef(line_FJ)))&&!all(is.infinite(coef(line_FJ)))) {
    intercept=sprintf("%.2f",coef(line_FJ,use.na=false)[1])
    slope=sprintf("%.2f",coef(line_FJ,use.na=false)[2])
    abline(line_FJ, col="red")
  }
  abline(0,1)
  legend("bottomright", legend = c(paste("slope = ", slope), paste("intercept = ", intercept), paste("R = ", R_seawifs.3), paste("RMSE = ", RMSE_seawifs.3)), bty="n") 

close.screen(all = TRUE)
}

# Graphic Function (3) Error Bars----------------------------------------------------------------------
earlier_current_GRAPHIC_ERR <- function(earlier, current, TIME, earlier_SD, current_SD, variable_info, typestr, site){

opar<-par( cex=0.6, pin=c(6.0, 2.0))

split.screen(c(2,1))
split.screen(c(1,1))
screen(1)
  plotCI(1:length(TIME), earlier, earlier_SD, pch=21, col="red",sfrac=0.001,gap=0,
    ylim=c(min(earlier-earlier_SD, current-current_SD, na.rm=TRUE),
           max(earlier+earlier_SD, current+current_SD, na.rm=TRUE)),
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""),
    main=paste("Standard Deviation of ", typestr, " ", variable_info[1], " in ", site, " (", variable_info[3],"-",variable_info[4], ")", sep=""))
    #lines(1:length(TIME), earlier, type="l", col="pink")

screen(3)
  plotCI(1:length(TIME),current, current_SD, pch=21,col="black",sfrac=0.001,gap=0,
    ylim=c(min(earlier-earlier_SD, current-current_SD, na.rm=TRUE),
           max(earlier+earlier_SD, current+current_SD, na.rm=TRUE)),
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""),
    main=paste("Standard Deviation of ", typestr, " ", variable_info[1], " in ", site, " (", variable_info[3],"-",variable_info[4], ")", sep=""))
    #lines(1:length(TIME), current, type="l", col="light grey")
    legend("topright", legend = c("earlier","current"),
    col=c("red","black"), pch=c(21,21), bty="n", ncol=2)
	
close.screen(all = TRUE)
}

# Aggregation by time routine --------------------------------------------------------------------------
earlier_current_COMPUTE_GROUPING <- function(Value, TIME, period, method){

  if (period=="perhour"){
    time <- format(as.POSIXct(as.Date(as.numeric(TIME-1/24), origin=starttime)),"%m-%H-%M")
  }

  if (period=="perhour_15"){
    time <- format(as.POSIXct(as.Date(as.numeric(TIME-1/24), origin=starttime)),"%m-%d-%H-%M")
  }
  
  if (period=="daily"){
    time <- format(as.Date(as.numeric(TIME),origin=starttime), "%Y-%m-%d")
  }

  if (period=="halfmonth"){
    time <- format(as.Date(as.numeric(TIME),origin=starttime), "%Y-%m-%d")
  }
  
  if (period=="monthly"){
    time <- format(as.Date(as.numeric(TIME),origin=starttime), "%Y-%m")
  }
 
  if (period=="annual"){
    time <- format(as.Date(as.numeric(TIME),origin=starttime), "%Y")
  }

  JF_V <- data.frame(time, Value)
  JF_AGG <- aggregate(JF_V$Value, list(JF_V$time), method, na.rm=TRUE)

  if (period=="perhour_15"){
    time <- JF_AGG$Group.1
    time <- paste(substr(time,1,3), ifelse(as.integer(substr(time,4,5)) %in% c(1:15), 15, 30), substr(time,6,11),sep="")
    JF_V <- data.frame(time, Value=JF_AGG$x)
    JF_AGG <- aggregate(JF_V$Value, list(JF_V$time), method, na.rm=TRUE)
  }
  
  if (period=="halfmonth"){
    time <- JF_AGG$Group.1
    time <- paste(substr(time,1,8), ifelse(as.integer(substr(time,9,10)) %in% c(1:15), 15, 30),sep="")
    JF_V <- data.frame(time, Value=JF_AGG$x)
    JF_AGG <- aggregate(JF_V$Value, list(JF_V$time), method, na.rm=TRUE)
  }
  
  return (JF_AGG)

}

# Different versions to  be compared--------------------------------------------------------------------------------------------
earlier_current_EVALUATION <- function(method, earlier, current){
  
  d1=current[!is.na(current)&!is.na(earlier)]
  d2=earlier[!is.na(current)&!is.na(earlier)]

if(length(d1)>0) {
# Root Mean Square Error
   if (method=="rmse") r_value <- sqrt(mean((d1-d2)^2, na.rm=TRUE))

   # Normalized RMSE
   if (method=="nrmse") r_value <- (sqrt(mean((d1-d2)^2, na.rm=TRUE)))/mean(d1)

   # Coefficient of determination (R Squared)
   if (method=="RR") r_value <- summary(lm(d2 ~ d1, na.action=na.omit))[c("r.squared")]

   # Pearson's correlation coefficient
   if (method=="cor_test") r_value <- cor.test(d1, d2, na.action=na.omit )$p.value
     #cor.test(FLUXNET_V, JSBACH_V, method="spearman")
     #cor.test(FLUXNET_V, JSBACH_V, method="kendall")

   # Nash-Sutchliffe efficiency (NSE)
   if (method=="nse") r_value <- 1- sum((d1 - d2)^2)/sum((d1 - mean(d2))^2)

   # RMSE-observations standard deviation ratio (RSR)
   if (method=="rsr") r_value <- sqrt(sum((d1 - d2)^2))/sqrt(sum((d1 - mean(d2))^2))

   # Percent bias (PBIAS)
   if (method=="pbias") r_value <- sum(d2 - d1)*100/sum(d1)

   # Modelling efficiency
   if (method=="mef") r_value <- (sum((d1-mean(d1))^2)-sum((d2-d1)^2))/sum((d1-mean(d1))^2)

   # Normalized average error
   if (method=="nae") r_value <- (mean(d2)-mean(d1))/mean(d1)

   # Variance Ratio
   if (method=="vr") r_value <- var(d2)/var(d1)

  } else {
   r_value <- NA
  }

  return(sprintf("%.2f", r_value))

}

# define work directory and load reqiured library
work_dir <- try(system("pwd",intern=TRUE))
setwd(work_dir)
library(ncdf)
library(gdata)
library(gplots)
library(plotrix)
library(Hmisc)
#source('TAYLOR.R')

earlier_version <- "rev224" 
current_version <- "rev224"
earlier_forcing <- "halfhourly"
current_forcing <- "halfhourly"
parameter_forcing1 <- "halfhourly"
parameter_forcing2 <- "halfhourly_40dp"
earlier_dir <- "/fluxnet_runs/224/halfhourly"
current_dir <- "/fluxnet_runs/224/40dp"
obs_dir <- "/FLUXNET/halfhourly"

if (earlier_version == current_version && earlier_forcing == current_forcing){
    earlier_diff <- paste(earlier_version, " ", parameter_forcing1)
    current_diff <- paste(earlier_version, " ", parameter_forcing2)
}else{
    earlier_diff <- paste(earlier_version, " ", earlier_forcing, sep="")
    current_diff <- paste(current_version, " ", current_forcing, sep="")
}
# define required compared variables, can be added later.
nvar=9
aggregation_level = c( "perhour", "daily", "monthly" )
eval_arg <- c("RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr")
annual_values=array(0,c(3,nvar+3))
eval_stats=array(0, dim=c(12,nvar,3,3))
variable_names <- array(c("par_acc", "PPFD", "PPFD", "PPFD",
                          "net_radiation", "Rn", "Net Radiation", "Rn",
                          "sensible_heat_flx", "H", "Sensible Heat Flux", "Qh" , 
                          "latent_heat_flx", "LE", "Latent Heat Flux", "Qle",  
                          "canopy_cond_limited", "gsurf", "Canopy Conductance", "gscan" ,
                          "net_assimilation", "GPP", "GPP", "GPP",
                          "reco", "Reco", "Reco", "Reco",
                          "net_co2_flux", "NEE", "NEE", "NEE",
                          "fapar", "fapar", "FAPAR", "fapar" )
                        , dim=c(4,nvar))
variable_units <- array(c("mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "W/m^2", "W/m^2", "W/m^2", 
                          "W/m^2", "W/m^2", "W/m^2", 
                          "W/m^2", "W/m^2", "W/m^2", 
                          "m/s","mmol m-2 s-1","m/s",
                          "mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "-","-","-")
                        , dim=c(3,nvar))
variable_units_conversion <- array(c(10^6,1,
                                     1,1,
                                     -1,1,
                                     -1,1,
                                      1,18*10^-6,
                                     10^6,1,
                                     10^6,1,
                                     10^6,1,
				     1,1)
                        , dim=c(2,nvar))

site_info <- read.table("core_sites.txt",header=TRUE,as.is = TRUE)
len_site <- length(site_info$SITE)
ab_var <- substr(variable_names[3,],1,1)
ab_site <- substr(site_info$SITE,1,1)

for (s in 1:1){
  pdf(paste(site_info$SITE[s], ".compare.", earlier_diff, "_2_", current_diff, ".pdf", sep=""))
  
  for (V in 1:nvar){
    site <- site_info$SITE[s]
	startyear <- site_info$START_YEAR[s]
    starttime <- paste(startyear,"-01-01",sep='')
    endyear <- site_info$END_YEAR[s]
	
    earlier_nc <- open.ncdf(paste(earlier_dir, "/", site_info$SITE[s], "/Output/", site_info$SITE[s], ".", earlier_version, ".", earlier_forcing, ".lsms.halfhourly.nc", sep="")) 
    current_nc <- open.ncdf(paste(current_dir, "/", site_info$SITE[s], "/Output/", site_info$SITE[s], ".", current_version, ".", current_forcing, ".lsms.halfhourly.nc", sep=""))	 
	datafile <- open.ncdf(paste(obs_dir, "/", site_info$SITE[s], ".", site_info$START_YEAR[s], "-", site_info$END_YEAR[s], ".obs.halfhourly.nc", sep=""))	
    
	#earlier_nc <- open.ncdf("DE-Hai.rev168.halfhourly.lsms.halfhourly.nc")	
	#current_nc <- open.ncdf("DE-Hai.rev181.halfhourly.lsms.halfhourly.nc")
	#datafile <- open.ncdf("DE-Hai.2000-2006.obs.halfhourly.nc")
	
	earlier <- get.var.ncdf(earlier_nc, variable_names[1,V]) 
    earlier <- earlier*(variable_units_conversion[1,V])	
    earlier[earlier<=-9999.]=NA
    current <- get.var.ncdf(current_nc, variable_names[1,V])
	current <- current*(variable_units_conversion[1,V])
    current[current<=-9999.]=NA
	TIME <- get.var.ncdf(earlier_nc, "time")
    
	if(variable_names[1,V]=='fapar'){
	  
	  FLUXNET_V <- get.var.ncdf( datafile, 'fapar_seawifs' )
      FLUXNET_V[FLUXNET_V==-9999.]=NA 
      FLUXNET_V_FLAG <- get.var.ncdf( datafile, 'fapar_seawifs_flag' )
      FLUXNET_V_FILTERED=FLUXNET_V
      FLUXNET_V_FILTERED[FLUXNET_V_FLAG>0.5]=NA
	  
      FLUXNET_V2 <- get.var.ncdf( datafile, 'fapar_modis' )
      FLUXNET_V2[FLUXNET_V2==-9999.]=NA 
      FLUXNET_V2_FLAG <- get.var.ncdf( datafile, 'fapar_modis_flag' )
      FLUXNET_V2_FILTERED=FLUXNET_V2
      FLUXNET_V2_FILTERED[FLUXNET_V2_FLAG>0.5]=NA
	  
      FLUXNET_V3 <- get.var.ncdf( datafile, 'fapar_cyclopes' )
      FLUXNET_V3[FLUXNET_V3==-9999.]=NA 
      FLUXNET_V3_FLAG <- get.var.ncdf( datafile, 'fapar_cyclopes_flag' )
      FLUXNET_V3_FILTERED=FLUXNET_V3
      FLUXNET_V3_FILTERED[FLUXNET_V3_FLAG>0.5]=NA
	    
    } 
	else{
        
	  FLUXNET_V <- get.var.ncdf( datafile, variable_names[2,V] )
      FLUXNET_V[FLUXNET_V==-9999.]=NA
      FLUXNET_V <- FLUXNET_V*(variable_units_conversion[2,V]) 
      FLUXNET_V_FILTERED=FLUXNET_V
        
	  if(!variable_names[2,V]=='Rn') {
        FLUXNET_V_FLAG <- get.var.ncdf( datafile, paste(variable_names[2,V],'_flag',sep='') )
        FLUXNET_V_FILTERED[FLUXNET_V_FLAG>1]=NA
      }
	
      if(any(!is.na(FLUXNET_V_FILTERED))) {  
	    # eval_stats[first(differnt statistics), sencond(differnt variables), third(differnt periods), fourth(differnt combination)]
		eval_stats[1,V,1,1]=sprintf("%.2f",coef(line(FLUXNET_V_FILTERED[!is.na(FLUXNET_V_FILTERED)], earlier[!is.na(FLUXNET_V_FILTERED)]),
                 use.na=false)[1])
        eval_stats[2,V,1,1]=sprintf("%.2f",coef(line(FLUXNET_V_FILTERED[!is.na(FLUXNET_V_FILTERED)], earlier[!is.na(FLUXNET_V_FILTERED)]),
                 use.na=false)[2])
				 
        eval_stats[1,V,1,2]=sprintf("%.2f",coef(line(FLUXNET_V_FILTERED[!is.na(FLUXNET_V_FILTERED)], current[!is.na(FLUXNET_V_FILTERED)]),
                 use.na=false)[1])
        eval_stats[2,V,1,2]=sprintf("%.2f",coef(line(FLUXNET_V_FILTERED[!is.na(FLUXNET_V_FILTERED)], current[!is.na(FLUXNET_V_FILTERED)]),
                 use.na=false)[2])	
				 
        eval_stats[1,V,1,3]=sprintf("%.2f",coef(line(earlier, current), use.na=false)[1])
        eval_stats[2,V,1,3]=sprintf("%.2f",coef(line(earlier, current), use.na=false)[2])					 
      } 
	
	  eval_stats[3,V,1,1]=sprintf("%.2f", cor(earlier, FLUXNET_V_FILTERED))
      eval_stats[3,V,1,2]=sprintf("%.2f", cor(current, FLUXNET_V_FILTERED))
	  eval_stats[3,V,1,3]=sprintf("%.2f", cor(earlier, current))
	  
	  eval_stats[4,V,1,1]=sprintf("%.2f",sd(FLUXNET_V_FILTERED,na.rm=TRUE))
	  dummy_e=earlier
      dummy_e[is.na(FLUXNET_V_FILTERED)]==NA
      eval_stats[4,V,1,2]=sprintf("%.2f",sd(dummy_e,na.rm=TRUE))

	  dummy_c=earlier
      dummy_c[is.na(FLUXNET_V_FILTERED)]==NA
      eval_stats[4,V,1,3]=sprintf("%.2f",sd(dummy_c,na.rm=TRUE))

	  for (N in 1:length(eval_arg)){
        eval_stats[N+4,V,1,1]=earlier_current_EVALUATION(eval_arg, earlier, FLUXNET_V_FILTERED)
		eval_stats[N+4,V,1,2]=earlier_current_EVALUATION(eval_arg, current, FLUXNET_V_FILTERED)
		eval_stats[N+4,V,1,3]=earlier_current_EVALUATION(eval_arg, earlier, current)
      } 
	  
	} # end else
	
    for (time in aggregation_level ){ 
	   
      if ( time != "annual" ) {
         
        variable_info <- c(variable_names[2,V],variable_units[2,V],startyear[s],endyear[s])
		FLUXNET_AGG <- earlier_current_COMPUTE_GROUPING(FLUXNET_V, TIME, time, 'mean')
        earlier_AGG <- earlier_current_COMPUTE_GROUPING(earlier, TIME, time, 'mean')
        current_AGG <- earlier_current_COMPUTE_GROUPING(current, TIME, time, 'mean')

		
		R.1 <- sprintf("%.2f", cor(earlier_AGG$x, FLUXNET_AGG$x))
	    R.2 <- sprintf("%.2f", cor(current_AGG$x, FLUXNET_AGG$x))
		R.3 <- sprintf("%.2f", cor(earlier_AGG$x, current_AGG$x))
		
		fluxnet_SD <- (earlier_current_COMPUTE_GROUPING(FLUXNET_V, TIME, time, 'sd'))$x
        earlier_SD <- (earlier_current_COMPUTE_GROUPING(earlier, TIME, time, 'sd'))$x
        current_SD <- (earlier_current_COMPUTE_GROUPING(current, TIME, time, 'sd'))$x
		
		COR_TEST.1 <- earlier_current_EVALUATION("cor_test", earlier_AGG$x, FLUXNET_AGG$x)
		COR_TEST.2 <- earlier_current_EVALUATION("cor_test", current_AGG$x, FLUXNET_AGG$x)
        COR_TEST.3 <- earlier_current_EVALUATION("cor_test", earlier_AGG$x, current_AGG$x)
		
		RMSE.1 <- earlier_current_EVALUATION("rmse", earlier_AGG$x, FLUXNET_AGG$x)
		RMSE.2 <- earlier_current_EVALUATION("rmse", current_AGG$x, FLUXNET_AGG$x)
		RMSE.3 <- earlier_current_EVALUATION("rmse", earlier_AGG$x, current_AGG$x)
		
		RR.1 <- earlier_current_EVALUATION("RR", earlier_AGG$x, FLUXNET_AGG$x)
		RR.2 <- earlier_current_EVALUATION("RR", current_AGG$x, FLUXNET_AGG$x)
        RR.3 <- earlier_current_EVALUATION("RR", earlier_AGG$x, current_AGG$x)
		
		NSE.1 <- earlier_current_EVALUATION("nse", earlier_AGG$x, FLUXNET_AGG$x)
		NSE.2 <- earlier_current_EVALUATION("nse", current_AGG$x, FLUXNET_AGG$x)
        NSE.3 <- earlier_current_EVALUATION("nse", earlier_AGG$x, current_AGG$x)
		
		PBIAS.1 <- earlier_current_EVALUATION("pbias", earlier_AGG$x, FLUXNET_AGG$x)
		PBIAS.2 <- earlier_current_EVALUATION("pbias", current_AGG$x, FLUXNET_AGG$x)
        PBIAS.3 <- earlier_current_EVALUATION("pbias", earlier_AGG$x, current_AGG$x)
        
		RSR.1 <- earlier_current_EVALUATION("rsr", earlier_AGG$x, FLUXNET_AGG$x)
		RSR.2 <- earlier_current_EVALUATION("rsr", current_AGG$x, FLUXNET_AGG$x)
		RSR.3 <- earlier_current_EVALUATION("rsr", earlier_AGG$x, current_AGG$x)
		
		MEF.1 <- earlier_current_EVALUATION("mef", earlier_AGG$x, FLUXNET_AGG$x)
		MEF.2 <- earlier_current_EVALUATION("mef", current_AGG$x, FLUXNET_AGG$x)
		MEF.3 <- earlier_current_EVALUATION("mef", earlier_AGG$x, current_AGG$x)
        
		NAE.1 <- earlier_current_EVALUATION("nae", earlier_AGG$x, FLUXNET_AGG$x)
		NAE.2 <- earlier_current_EVALUATION("nae", current_AGG$x, FLUXNET_AGG$x)
		NAE.3 <- earlier_current_EVALUATION("nae", earlier_AGG$x, current_AGG$x)
        
		VR.1 <- earlier_current_EVALUATION("vr", earlier_AGG$x, FLUXNET_AGG$x)
		VR.2 <- earlier_current_EVALUATION("vr", current_AGG$x, FLUXNET_AGG$x)
		VR.3 <- earlier_current_EVALUATION("vr", earlier_AGG$x, current_AGG$x)         
		
        if (time == "daily" || time == "monthly" ) {
   
          if(time=="daily") si=2
          if(time=="monthly") si=3

		  eval_stats[1,V,si,1]=sprintf("%.2f",coef(line(earlier_AGG$x, FLUXNET_AGG$x),use.na=false)[1])
          eval_stats[2,V,si,1]=sprintf("%.2f",coef(line(earlier_AGG$x, FLUXNET_AGG$x),use.na=false)[2])
				 
          eval_stats[1,V,si,2]=sprintf("%.2f",coef(line(current_AGG$x, FLUXNET_AGG$x),use.na=false)[1])
          eval_stats[2,V,si,2]=sprintf("%.2f",coef(line(current_AGG$x, FLUXNET_AGG$x),use.na=false)[2])	
				 
          eval_stats[1,V,si,3]=sprintf("%.2f",coef(line(earlier_AGG$x, current_AGG$x), use.na=false)[1])
          eval_stats[2,V,si,3]=sprintf("%.2f",coef(line(earlier_AGG$x, current_AGG$x), use.na=false)[2])					 
     
	
	      eval_stats[3,V,si,1]=sprintf("%.2f", cor(earlier_AGG$x, FLUXNET_AGG$x))
          eval_stats[3,V,si,2]=sprintf("%.2f", cor(current_AGG$x, FLUXNET_AGG$x))
	      eval_stats[3,V,si,3]=sprintf("%.2f", cor(earlier_AGG$x, current_AGG$x))
	  
	      eval_stats[4,V,si,1]=sprintf("%.2f",sd(FLUXNET_AGG$x,na.rm=TRUE))
          eval_stats[4,V,si,2]=sprintf("%.2f",sd(earlier_AGG$x,na.rm=TRUE))
          eval_stats[4,V,si,3]=sprintf("%.2f",sd(current_AGG$x,na.rm=TRUE))
		  for (N in 1:length(eval_arg)){
            eval_stats[N+4,V,si,1]=earlier_current_EVALUATION(eval_arg, earlier_AGG$x, FLUXNET_AGG$x)
		    eval_stats[N+4,V,si,2]=earlier_current_EVALUATION(eval_arg, current_AGG$x, FLUXNET_AGG$x)
		    eval_stats[N+4,V,si,3]=earlier_current_EVALUATION(eval_arg, earlier_AGG$x, current_AGG$x)
          }
        }
		 
		if (time=="perhour") {
               typestr <- c("Per Hour Monthly")
               gtype <- c("b")
			   labs <- month.abb[as.numeric(gsub("0","",substr(earlier_AGG$Group.1, 1, 2)))]
        }
        else if (time=="perhour_15") {
               typestr <- c("Per Hour 15 days")
               gtype <- c("p")
			   labs <- month.abb[as.numeric(gsub("0","",substr(earlier_AGG$Group.1, 1, 2)))]
        }
        else if (time=="daily"){
               typestr <- c("Daily Mean")
               gtype <- c("p")
			   labs <- substr(earlier_AGG$Group.1, 1, 4) 
        }
	    else if (time=="halfmonth"){
               typestr <- c("Half Month")
               gtype <- c("b")
			   labs <- substr(earlier_AGG$Group.1, 1, 4)
        }
        else if (time=="monthly"){           
               typestr <- c("Monthly Mean")
               gtype <- c("b")
			   labs <- substr(earlier_AGG$Group.1, 1, 4)
        }
        else if (time=="annual"){           
               annual_values[1,V]=mean(earlier,na.rm=TRUE)
               annual_values[2,V]=mean(current,na.rm=TRUE)
			   labs <- format(as.Date(as.numeric(earlier_AGG$Group.1),origin=starttime), "%Y")
        }       
		
         if(variable_names[1,V]=="fapar") {
           if (time != "perhour" ) {
             FLUXNET_AGG2 <- earlier_current_COMPUTE_GROUPING(FLUXNET_V2,TIME, time,'mean')
             FLUXNET_AGG3 <- earlier_current_COMPUTE_GROUPING(FLUXNET_V3,TIME, time,'mean')

			 R_seawifs.1 <- sprintf("%.2f", cor(earlier_AGG$x, FLUXNET_AGG$x))
	         R_seawifs.2 <- sprintf("%.2f", cor(current_AGG$x, FLUXNET_AGG$x))
		     R_seawifs.3 <- sprintf("%.2f", cor(earlier_AGG$x, current_AGG$x))
		     
			 R_modis.1 <- sprintf("%.2f", cor(earlier_AGG$x, FLUXNET_AGG2$x))
	         R_modis.2 <- sprintf("%.2f", cor(current_AGG$x, FLUXNET_AGG2$x))
		     R_modis.3 <- sprintf("%.2f", cor(earlier_AGG$x, current_AGG$x))
			 
			 R_cyclopes.1 <- sprintf("%.2f", cor(earlier_AGG$x, FLUXNET_AGG3$x))
	         R_cyclopes.2 <- sprintf("%.2f", cor(current_AGG$x, FLUXNET_AGG3$x))
		     R_cyclopes.3 <- sprintf("%.2f", cor(earlier_AGG$x, current_AGG$x))
			 
			 RMSE_seawifs.1 <- earlier_current_EVALUATION("rmse", earlier_AGG$x, FLUXNET_AGG$x)
			 RMSE_seawifs.2 <- earlier_current_EVALUATION("rmse", current_AGG$x, FLUXNET_AGG$x)
             RMSE_seawifs.3 <- earlier_current_EVALUATION("rmse", earlier_AGG$x, current_AGG$x)
			 
			 RMSE_modis.1 <- earlier_current_EVALUATION("rmse", earlier_AGG$x, FLUXNET_AGG2$x)
			 RMSE_modis.2 <- earlier_current_EVALUATION("rmse", current_AGG$x, FLUXNET_AGG2$x)
             RMSE_modis.3 <- earlier_current_EVALUATION("rmse", earlier_AGG$x, current_AGG$x)
			 
			 RMSE_cyclopes.1 <- earlier_current_EVALUATION("rmse", earlier_AGG$x, FLUXNET_AGG3$x)
			 RMSE_cyclopes.2 <- earlier_current_EVALUATION("rmse", current_AGG$x, FLUXNET_AGG3$x)
             RMSE_cyclopes.3 <- earlier_current_EVALUATION("rmse", earlier_AGG$x, current_AGG$x)

             RR_seawifs.1 <- earlier_current_EVALUATION("RR", earlier_AGG$x, FLUXNET_AGG$x)
			 RR_seawifs.2 <- earlier_current_EVALUATION("RR", current_AGG$x, FLUXNET_AGG$x)
             RR_seawifs.3 <- earlier_current_EVALUATION("RR", earlier_AGG$x, current_AGG$x)
			 
			 RR_modis.1 <- earlier_current_EVALUATION("RR", earlier_AGG$x, FLUXNET_AGG2$x)
			 RR_modis.2 <- earlier_current_EVALUATION("RR", current_AGG$x, FLUXNET_AGG2$x)
             RR_modis.3 <- earlier_current_EVALUATION("RR", earlier_AGG$x, current_AGG$x)
			 
			 RR_cyclopes.1 <- earlier_current_EVALUATION("RR", earlier_AGG$x, FLUXNET_AGG3$x)
			 RR_cyclopes.2 <- earlier_current_EVALUATION("RR", current_AGG$x, FLUXNET_AGG3$x)
             RR_cyclopes.3 <- earlier_current_EVALUATION("RR", earlier_AGG$x, current_AGG$x)
			 
             earlier_current_GRAPHIC_FAPAR(earlier_AGG$x, current_AGG$x, FLUXNET_AGG$x, FLUXNET_AGG2$x, FLUXNET_AGG3$x, earlier_AGG$Group.1, typestr, gtype, variable_info, site, R_seawifs.1, R_seawifs.2, R_seawifs.3, RMSE_seawifs.1, RMSE_seawifs.2, RMSE_seawifs.3)
             
           }
          } else {
             earlier_current_GRAPHIC(earlier_AGG$x, current_AGG$x, FLUXNET_AGG$x, earlier_AGG$Group.1, typestr, gtype, labs, variable_info, site, R.1, R.2, R.3, RMSE.1, RMSE.2, RMSE.3)
             #earlier_current_GRAPHIC_ERR(earlier_AGG$x, current_AGG$x, FLUXNET_AGG$x, earlier_AGG$Group.1, earlier_SD, current_SD, fluxnet_SD, variable_info, typestr, site)
			 			 
            }			
          }
        }
    }
	
	# for halfhourly data
	write.table(t(eval_stats[,,1,1]),file=paste(site_info$SITE[s],".halfhourly.",earlier_diff, "_2_", current_diff,".txt",sep=""),col.names=c("Intercept", "Slope", "R", "SD", "RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr"), row.names=variable_names[2,])
    write.table(t(eval_stats[,,1,2]),file=paste(site_info$SITE[s],".halfhourly.",earlier_diff, "_2_", current_diff,".txt",sep=""),col.names=c("Intercept", "Slope", "R", "SD", "RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr"), row.names=variable_names[2,])
    write.table(t(eval_stats[,,1,3]),file=paste(site_info$SITE[s],".halfhourly.",earlier_diff, "_2_", current_diff,".txt",sep=""),col.names=c("Intercept", "Slope", "R", "SD", "RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr"), row.names=variable_names[2,])
	
	# for daily data
	write.table(t(eval_stats[,,2,1]),file=paste(site_info$SITE[s],".daily.",earlier_diff, "_2_", current_diff,".txt",sep=""),col.names=c("Intercept", "Slope", "R", "SD", "RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr"), row.names=variable_names[2,])
    write.table(t(eval_stats[,,2,2]),file=paste(site_info$SITE[s],".daily.",earlier_diff, "_2_", current_diff,".txt",sep=""),col.names=c("Intercept", "Slope", "R", "SD", "RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr"), row.names=variable_names[2,])
    write.table(t(eval_stats[,,2,3]),file=paste(site_info$SITE[s],".daily.",earlier_diff, "_2_", current_diff,".txt",sep=""),col.names=c("Intercept", "Slope", "R", "SD", "RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr"), row.names=variable_names[2,])
	
	# for monthly data
	write.table(t(eval_stats[,,3,1]),file=paste(site_info$SITE[s],".monthly.",earlier_diff, "_2_", current_diff,".txt",sep=""),col.names=c("Intercept", "Slope", "R", "SD", "RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr"), row.names=variable_names[2,])
    write.table(t(eval_stats[,,3,2]),file=paste(site_info$SITE[s],".monthly.",earlier_diff, "_2_", current_diff,".txt",sep=""),col.names=c("Intercept", "Slope", "R", "SD", "RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr"), row.names=variable_names[2,])
    write.table(t(eval_stats[,,3,3]),file=paste(site_info$SITE[s],".monthly.",earlier_diff, "_2_", current_diff,".txt",sep=""),col.names=c("Intercept", "Slope", "R", "SD", "RR", "nse", "rmse", "pbias", "rsr", "mef", "nae", "vr"), row.names=variable_names[2,])

	close.ncdf(datafile)
	close.ncdf(earlier_nc)
	close.ncdf(current_nc)
}

warnings()


