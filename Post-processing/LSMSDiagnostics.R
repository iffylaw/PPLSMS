#!/usr/bin/env Rscript
# Date: 25/02/2014
# This script generates some simple diagnostics and plot for LSMS fluxnet simulations
# contains:
# LSMS_FLUXNET_GRAPHIC - function for ordinary fluxnet plots
# LSMS_FLUXNET_GRAPHIC_FAPAR - function for satellite comparison plots
# LSMS_FLUXNET_GRAPHIC_ERR - function for error bars 
# LSMS_FLUXNET_COMPUTE_GROUPING - function to aggregate to daily, monthly, annual and diurnal cycle
# LSMS_FLUXNET_EVALUATION - function to compute Root Mean Square Error, R sqaure, NSE, PBIAS, RSR
# LSMS_FLUXNET_GRAPHIC_FORCING - function to compare forcing data with observation 
# 
# main 

# List of things to do
# strage behaviour of first plot (black squares

#---------------------------------------------------------------------------------------------------
rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows   

# Set parameters
Args <- commandArgs(trailingOnly=TRUE);
if(length(Args) != 2) {
  message("generate_FLUXNET_diagnostics.R requires startyear endyear as input. Terminating");quit()
}

# Graphic Function (1) FLUXNET STUFF-----------------------------------------------------------------
LSMS_FLUXNET_GRAPHIC <- function(LSMS_V, FLUXNET_V, TIME, mean, gtype, labs, variable_info, 
   RMSE, NRMSE, RR, NSE, PBIAS, RSR){

opar<-par( cex=0.6, pin=c(6.0, 2.0))

split.screen(c(2,1))
split.screen(c(1,2),screen=2)
screen(1)
  plot(TIME, LSMS_V, type='n',  col="red", 
    ylim=c(min(LSMS_V, FLUXNET_V,na.rm=TRUE), max(LSMS_V, FLUXNET_V,na.rm=TRUE)), 
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
    main=paste(typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")", sep=""),
    las=1, axes=F)
  axis(1, at=1:length(TIME), lab=labs,las=1)
  axis(2,las=1)
  box() 
  lines(TIME, FLUXNET_V, type=gtype, col="black")
  lines(TIME, LSMS_V, type=gtype, col="red" ,pch=19)
  legend("topleft", legend = c("LSMS","FLUXNET"),  col=c("red","black"), pch=c(19,21), bty="n", ncol=2)
  legend("topright", legend = c(paste("NSE = ", NSE), paste("PBIAS = ", PBIAS, "%", sep=""), paste("RSR = ", RSR)), bty="n")

screen(3)
  plot(FLUXNET_V,LSMS_V, xlab='Fluxnet',ylab='LSMS',las=1,
      ylim=c(min(LSMS_V, FLUXNET_V,na.rm=TRUE), max(LSMS_V, FLUXNET_V,na.rm=TRUE)),
      xlim=c(min(LSMS_V, FLUXNET_V,na.rm=TRUE), max(LSMS_V, FLUXNET_V,na.rm=TRUE)),
      main= paste(variable_info[1], "(" ,  variable_info[2], ")", sep= "" ))
  if(!all(is.na(FLUXNET_V))) {
    line_FJ <- line(FLUXNET_V,LSMS_V) 
    a=NA;b=NA
    if(!any(is.na(coef(line_FJ)))) {
      a=sprintf("%.2f",coef(line_FJ,use.na=false)[1])
      b=sprintf("%.2f",coef(line_FJ,use.na=false)[2])
      abline(line_FJ, col="red")
    }
    abline(0,1)
    legend("bottomright", legend = c(paste("Int = ", a), paste("slope = ", b), paste("RMSE = ", RMSE), 
      paste("NRMSE = ", NRMSE), as.expression(bquote(R^2 == .(RR)))), bty="n") 
  }
screen(4)
  #JF_T <- c(rep("LSMS",length(TIME)),rep("FLUXNET",length(TIME)))
  #JF_V <- c(LSMS_V, FLUXNET_V)
  #bwplot(JF_T ~ JF_V, panel=panel.bpplot, nout=.05, scat1d.opts=list(frac=.01))
  #bpplot(LSMS_V, FLUXNET_V, name=c("LSMS", "FLUXNET"))
  boxplot(LSMS_V, FLUXNET_V, names=c("LSMS", "FLUXNET"),las=1,
    ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
    main = paste("Boxplot of ", variable_info[1]))

close.screen(1:4)
}

# Graphic Function (2) SATELLITE STUFF-----------------------------------------------------------------
LSMS_FLUXNET_GRAPHIC_FAPAR <- function(LSMS_V, FLUXNET_V, FLUXNET_V2, FLUXNET_V3, TIME, typestr, gtype, variable_info, 
  RMSE_seawifs, RMSE_modis, RMSE_cyclopes, RR_seawifs, RR_modis, RR_cyclopes, NSE, PBIAS, RSR){

opar<-par( cex=0.6, pin=c(6.0, 2.0))

split.screen(c(2,1))
split.screen(c(1,2),screen=2)
screen(1)
  plot(TIME, LSMS_V, type=gtype,  col="red", pch=19,
    ylim=c(min(LSMS_V, c(FLUXNET_V,FLUXNET_V2,FLUXNET_V3),na.rm=TRUE), 
           max(LSMS_V, c(FLUXNET_V,FLUXNET_V2,FLUXNET_V3),na.rm=TRUE)), 
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
    main=paste(typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")", sep=""),
    las=1)
  lines(TIME, FLUXNET_V, type=gtype, col="black" )
  lines(TIME, FLUXNET_V2, type=gtype, col="blue" )
  lines(TIME, FLUXNET_V3, type=gtype, col="grey" )
  lines(TIME, LSMS_V, type=gtype, col="red" ,pch=19)
  legend("bottomright", legend = c("LSMS","SeaWifs","Modis","Cyclopes"),  
    col=c("red","black","blue","grey"), pch=c(19,21,21,21), bty="n", ncol=4)
  legend("topright", legend = c(paste("NSE = ", NSE), paste("PBIAS = ", PBIAS, "%", sep=""), paste("RSR = ", RSR)), bty="n")

screen(3)
  plot(FLUXNET_V,LSMS_V, xlab='SeaWifs',ylab='LSMS',las=1,xlim=c(0,1),ylim=c(0,1),
      main= paste(variable_info[1], "(" ,  variable_info[2], ")", sep= "" ))
  legend("topleft", legend = c(paste("RMSE (SeaWifs) = ", RMSE_seawifs), paste("RMSE (Modis) = ", RMSE_modis), paste("RMSE (cyclopes) = ", RMSE_cyclopes),
    as.expression(bquote(R^2 (SeaWifs) == .(RR_seawifs))), as.expression(bquote(R^2 (Modis) == .(RR_modis))), as.expression(bquote(R^2 (Cyclopes) == .(RR_cyclopes)))), bty="n")
  points(FLUXNET_V2, LSMS_V,col='blue',pch=1)
  points(FLUXNET_V3, LSMS_V,col='grey',pch=15)

screen(4)
  boxplot(LSMS_V, FLUXNET_V, names=c("LSMS", "FLUXNET"), col = "light grey", las=1,
    ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
    main = paste("Boxplot of ", variable_info[1]))

close.screen(1:4)
}

# Graphic Function (3) compare Forcing data with simulation and evaluation------------------------
LSMS_FLUXNET_GRAPHIC_FORCING <- function(LSMS_V, FLUXNET_V, FORCING_V, Parameter_info, TIME, typestr, gtype, variable_info){

  opar<-par( cex=0.6, pin=c(6.0, 2.0))
  plot(TIME, LSMS_V, type="n",  col="green", 
    ylim=c(min(LSMS_V, FLUXNET_V,na.rm=TRUE), max(LSMS_V, FLUXNET_V,na.rm=TRUE)), 
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
    main=paste(typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")", sep=""),
    las=1)
  lines(TIME, LSMS_V, lwd=2, xlab='Time', ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), ylim=c(min(LSMS_V, FLUXNET_V,na.rm=TRUE), max(LSMS_V, FLUXNET_V,na.rm=TRUE)), col="green")
  lines(TIME, FLUXNET_V, type=gtype, col="blue")
 
  par(new=T)
  plot(TIME, FORCING_V, col='black', lwd=2, ann=F, axes=F)
  mtext(Parameter_info, side=4, line=3, col='black')
  axis(4, col.axis='black', col='black')
  legend("topleft", legend = c("LSMS","FLUXNET","FORCING"),  col=c("green","blue","black"), pch=19, bty="n", ncol=2)
 
}

# Graphic Function (4) Error Bars----------------------------------------------------------------------
LSMS_FLUXNET_GRAPHIC_ERR <- function(LSMS_V, FLUXNET_V, TIME, LSMS_SD, FLUXNET_SD, variable_info, typestr){

opar<-par( cex=0.6, pin=c(6.0, 2.0))

split.screen(c(2,1))

screen(1)
  plotCI(1:length(TIME), LSMS_V, LSMS_SD, pch=21, col="red",sfrac=0.001,gap=0,
    ylim=c(min(LSMS_V-LSMS_SD, FLUXNET_V-FLUXNET_SD, na.rm=TRUE),
           max(LSMS_V+LSMS_SD, FLUXNET_V+FLUXNET_SD, na.rm=TRUE)),
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""),
    main=paste("Standard Deviation of ", typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")", sep=""))


  plotCI(1:length(TIME),FLUXNET_V, FLUXNET_SD, pch=21,col="black",sfrac=0.001,gap=0, add=TRUE,
    ylim=c(min(LSMS_V-LSMS_SD, FLUXNET_V-FLUXNET_SD, na.rm=TRUE),
           max(LSMS_V+LSMS_SD, FLUXNET_V+FLUXNET_SD, na.rm=TRUE)),
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""),
    main=paste("Standard Deviation of ", typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")", sep=""))
    #lines(1:length(TIME), FLUXNET_V, type="l", col="light grey")
    legend("topleft", legend = c("LSMS","FLUXNET"),
    col=c("red","black"), pch=c(21,21), bty="n", ncol=2)

screen(2)
  if (typestr=="Per Hour Monthly"){
     TIME <- substr(TIME, 4, 5)
     at <- 1:24
     labs <- 1:24
  }
  if (typestr=="Half Month"){
     TIME <- substr(TIME, 6, 10)
     at <- 1:24
     labs <-NA
     for (V in 1:12) labs<-c(labs,rep(month.abb[V],2))
     labs <- labs[2:25] 
  }
  if (typestr=="Monthly Mean"){
     TIME <- substr(TIME, 6, 7)
     at <- 1:12
     labs <- month.abb
  }
  if (typestr=="Per Hour Monthly" || typestr=="Half Month" || typestr=="Monthly Mean"){
    LSMS_V <- data.frame(TIME,LSMS_V)
    FLUXNET_V <- data.frame(TIME,FLUXNET_V)
    boxplot(LSMS_V ~ TIME, LSMS_V, col = "tomato", boxwex = 0.25, at = at - 0.2, 
      ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
      main = paste(typestr, variable_info[1]),axes=F)
    axis(1, at=at, lab=labs)
    axis(2)
    box()   
    boxplot(FLUXNET_V ~ TIME, FLUXNET_V, add= TRUE, col = "orange", boxwex = 0.25, at = at + 0.2, axes=F,
      ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
      main = paste(typestr, variable_info[1])) 
    smartlegend(x="left",y="top", inset = 0, c("LSMS", "FLUXNET"), fill = c("tomato", "orange"), bty="n")
  }
close.screen(all = TRUE)
}

# Graphic Function (5) Variables Statistics----------------------------------------------------------------------
LSMS_FLUXNET_STATS_GRAPHIC <- function(stat_n, period){
  def.par <- par(no.readonly = TRUE)
  nf <- layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE), heights=c(1,1,1), widths=c(3,1))
  times <- length(period)
    
  # JUST FOR 3 STATISTICS
  for (V in 1:3){
    s_pos <- which(stat_name==stat_n[V])
    ymin <- min(as.numeric(eval_stats[s_pos,,]), na.rm=TRUE)
    ymax <- max(as.numeric(eval_stats[s_pos,,]), na.rm=TRUE)
    plot(1:nvar, eval_stats[s_pos,,which(aggregation_level==period[1])], ylim=c(ymin, ymax), xlab="", ylab=toupper(stat_n[V]), main="Variables Statistics", type="o", col=rainbow(times)[1], xaxt = "n", yaxt="n")
    axis(1, at=1:nvar, lab=variable_names[4,1:nvar])
    axis(2)
    box()
    if (times >2){
      for (p in 2:times) {
           lines(1:nvar, eval_stats[s_pos,,which(aggregation_level==period[p])], type="o", col=rainbow(times)[p])
      }
    }
    legend("topright", legend = period, col=rainbow(times), pch=21, bty="n")
    boxplot(as.numeric(as.vector(eval_stats[s_pos,,])))	
  }
  par(def.par)
}

# Graphic Function (5) Taylor Diagram----------------------------------------------------------------------
LSMS_FLUXNET_TAYLOR <- function(modelfile, datafile, time){
  ab_var <- substr(variable_names[3,],1,1)
  
  for (V in 1:nvar){     
    
    LSMS_V <- LSMS_DATA(modelfile, V)
    FLUXNET_V <- FLUXNET_DATA(datafile, V)$FLUXNET_V

    LSMS_AGG <- LSMS_FLUXNET_COMPUTE_GROUPING(LSMS_V, TIME, time,'mean')$x
    FLUXNET_AGG <- LSMS_FLUXNET_COMPUTE_GROUPING(FLUXNET_V, TIME, time,'mean')$x
    
    if (time == "halfhourly"){
       LSMS_AGG <- LSMS_V
       FLUXNET_AGG <- FLUXNET_V
    }

    if (V==1){ 
      ad=FALSE 
    }else{
      ad=TRUE 
    }
    
    graph <- taylor.diagram(FLUXNET_AGG, LSMS_AGG, add=ad, pch=ab_var[V], col=rainbow(nvar)[V], normalize=TRUE, main=time)
   
    if (V==1){
       legend('bottomleft', variable_names[3, ], pch=ab_var, col=rainbow(nvar), cex=0.8, bty="n")
    }
  }
}

# Data Function (1) Get LSMS Data-----------------------------------------------------------------------
LSMS_DATA <- function(modelfile, V){
  LSMS_V <- get.var.ncdf( modelfile, variable_names[1,V] )
  LSMS_V <- LSMS_V*(variable_units_conversion[1,V])       
  LSMS_V[LSMS_V<=-9999.]=NA
  return(LSMS_V)
}

# Data Function (2) Get FLUXNET Data-----------------------------------------------------------------------
FLUXNET_DATA <- function(datafile, V){
  if(variable_names[1,V]=='fapar') {
     FLUXNET_V <- get.var.ncdf( datafile, 'fapar_seawifs' )
     FLUXNET_V[FLUXNET_V==-9999.]=NA 
     FLUXNET_V_FLAG <- get.var.ncdf( datafile, 'fapar_seawifs_flag' )
     FLUXNET_V[FLUXNET_V_FLAG>0.5]=NA
     FLUXNET_V2 <- get.var.ncdf( datafile, 'fapar_modis' )
     FLUXNET_V2[FLUXNET_V2==-9999.]=NA 
     FLUXNET_V2_FLAG <- get.var.ncdf( datafile, 'fapar_modis_flag' )
     FLUXNET_V2[FLUXNET_V2_FLAG>0.5]=NA
     FLUXNET_V3 <- get.var.ncdf( datafile, 'fapar_cyclopes' )
     FLUXNET_V3[FLUXNET_V3==-9999.]=NA 
     FLUXNET_V3_FLAG <- get.var.ncdf( datafile, 'fapar_cyclopes_flag' )
     FLUXNET_V3[FLUXNET_V3_FLAG>0.5]=NA
   } else {
     FLUXNET_V <- get.var.ncdf( datafile, variable_names[2,V] )
     FLUXNET_V[FLUXNET_V==-9999.]=NA
     FLUXNET_V <- FLUXNET_V*(variable_units_conversion[2,V]) 
     FLUXNET_V_FILTERED=FLUXNET_V
     if(!variable_names[2,V]=='Rn') {
        FLUXNET_V_FLAG <- get.var.ncdf( datafile, paste(variable_names[2,V],'_flag',sep='') )
        FLUXNET_V[FLUXNET_V_FLAG>1]=NA
     }
     FLUXNET_V2=NA
     FLUXNET_V3=NA
  }
  list(FLUXNET_V=FLUXNET_V, FLUXNET_V2=FLUXNET_V2, FLUXNET_V3=FLUXNET_V3)   
}

# Aggregation by time routine --------------------------------------------------------------------------
LSMS_FLUXNET_COMPUTE_GROUPING <- function(Value, TIME, period, method){

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

# Model Evalauation--------------------------------------------------------------------------------------------
LSMS_FLUXNET_EVALUATION <- function(method, LSMS_V, FLUXNET_V){
  
  d1=FLUXNET_V[!is.na(FLUXNET_V)&!is.na(LSMS_V)]
  d2=LSMS_V[!is.na(FLUXNET_V)&!is.na(LSMS_V)]

  if(length(d1)>0) {

   # Intercept
   if (method=="intercept")  r_value <- coef(line(d1, d2),use.na=false)[1]

   # Slope
   if (method=="slope")  r_value <- coef(line(d1, d2),use.na=false)[2]
   
   # R
   if (method=="r")  r_value <- cor(d2, d1,use="na.or.complete")
   
   # Standard Deviation of LSMS
   if (method=="sd_jsb")  r_value <- sd(d2, na.rm=TRUE)
   
   # Standard Deviation of FLUXNET
   if (method=="sd_flx")  r_value <- sd(d1, na.rm=TRUE)

   # Root Mean Square Error
   if (method=="rmse") r_value <- sqrt(mean((d1-d2)^2, na.rm=TRUE))

   # Normalized RMSE
   if (method=="nrmse") r_value <- (sqrt(mean((d1-d2)^2, na.rm=TRUE)))/mean(d1)

   # Coefficient of determination (R Squared)
   if (method=="rr") r_value <- summary(lm(d2 ~ d1, na.action=na.omit))[c("r.squared")]
  
   # Pearson's correlation coefficient
   if (method=="cor_test") r_value <- cor.test(d1, d2, na.action=na.omit )$p.value
     #cor.test(FLUXNET_V, LSMS_V, method="spearman")
     #cor.test(FLUXNET_V, LSMS_V, method="kendall")

   # Nash-Sutchliffe efficiency (NSE)
   if (method=="nse") r_value <- 1- sum((d1 - d2)^2)/sum((d1 - mean(d2))^2)

   # RMSE-observations standard deviation ratio (RSR)
   if (method=="rsr") r_value <- sqrt(sum((d1 - d2)^2))/sqrt(sum((d1 - mean(d2))^2)) 

   # Percent bias (PBIAS) 
   if (method=="pbias") r_value <- sum(d1 - d2)*100/sum(d1)

   # Modelling efficiency
   if (method=="mef") r_value <- (sum((d1-mean(d1))^2)-sum((d2-d1)^2))/sum((d1-mean(d1))^2)

   # Normalized average error
   if (method=="nae") r_value <- (mean(d2)-mean(d1))/mean(d1)
  
   # Variance Ratio
   if (method=="vr") r_value <- var(d2)/var(d1) 
   
   # Index of Agreement 
   if (method=="d") r_value <- d(d2, d1, na.rm=TRUE) 
    
  } else {
   r_value <- NA 
  }

  return(sprintf("%.2f", r_value))
}

#--- Here comes the main code -----------------------------------------------------------------------------------
startyear <- as.numeric(Args[1])
starttime <- paste(Args[1],"-01-01",sep='')
endyear <- as.numeric(Args[2])
work_dir <- try(system("pwd",intern=TRUE))

setwd(work_dir)
require(ncdf)
require(gdata)
require(gplots)
require(plotrix)
require(Hmisc)
require(lattice)
require(hydroGOF)
#source('TAYLOR.R')

# define required compared variables, can be added later.
nvar=9
aggregation_level = c("halfhourly", "perhour", "daily", "halfmonth", "monthly", "annual" )
stat_name = c("intercept","slope", "r", "rr", "rmse", "nrmse", "sd_jsb", "sd_flx", "mef", "nae", "vr", "pbias", "nse", "rsr", "d" )
annual_values=array(NA,c(3,nvar+3))
eval_stats=array(NA, dim=c(15,nvar,5))
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

datafile <- open.ncdf(paste("obs.halfhourly.nc", sep=""))
modelfile <- open.ncdf(paste("LSMS.halfhourly.nc", sep=""))
TIME <- get.var.ncdf( modelfile, "time")
TIME_F <- get.var.ncdf( datafile, "time" )

if(length(TIME) == length(TIME_F)) {
  for (V in 1:nvar){   
      print(variable_names[1,V])
      LSMS_V <- LSMS_DATA(modelfile, V)
      FLUXNET_V <- FLUXNET_DATA(datafile, V)$FLUXNET_V
      FLUXNET_V2 <- FLUXNET_DATA(datafile, V)$FLUXNET_V2
      FLUXNET_V3 <- FLUXNET_DATA(datafile, V)$FLUXNET_V3
      
      if(any(!is.na(FLUXNET_V))&& !variable_names[1,V]=='fapar') {
         for(s in 1:length(stat_name)){
            eval_stats[s,V,1]=LSMS_FLUXNET_EVALUATION(stat_name[s], LSMS_V, FLUXNET_V)
         }
      }

      si=1
      for (time in aggregation_level[2:length(aggregation_level)] ){ 

        if ( time != "annual" ) {
         
          variable_info <- c(variable_names[3,V],variable_units[3,V],startyear,endyear)
          LSMS_AGG <- LSMS_FLUXNET_COMPUTE_GROUPING(LSMS_V, TIME, time,'mean')
          FLUXNET_AGG <- LSMS_FLUXNET_COMPUTE_GROUPING(FLUXNET_V, TIME, time,'mean')

          RMSE <- LSMS_FLUXNET_EVALUATION("rmse", LSMS_AGG$x, FLUXNET_AGG$x)
          NRMSE <- LSMS_FLUXNET_EVALUATION("nrmse", LSMS_AGG$x, FLUXNET_AGG$x)
          RR <- LSMS_FLUXNET_EVALUATION("rr", LSMS_AGG$x, FLUXNET_AGG$x)
          LSMS_SD <- (LSMS_FLUXNET_COMPUTE_GROUPING(LSMS_V,TIME, time,'sd'))$x
          FLUXNET_SD <- (LSMS_FLUXNET_COMPUTE_GROUPING(FLUXNET_V,TIME, time,'sd'))$x
          NSE <- LSMS_FLUXNET_EVALUATION("nse", LSMS_AGG$x, FLUXNET_AGG$x)
          PBIAS <- LSMS_FLUXNET_EVALUATION("pbias", LSMS_AGG$x, FLUXNET_AGG$x)
          RSR <- LSMS_FLUXNET_EVALUATION("rsr", LSMS_AGG$x, FLUXNET_AGG$x)
          COR_TEST <- LSMS_FLUXNET_EVALUATION("cor_test", LSMS_AGG$x, FLUXNET_AGG$x)

          si=si+1
          for(s in 1:length(stat_name)){
             eval_stats[s,V,si]=LSMS_FLUXNET_EVALUATION(stat_name[s], LSMS_AGG$x, FLUXNET_AGG$x)
          }

          if (time=="perhour") {
               typestr <- c("Per Hour Monthly")
               gtype <- c("b")
	           labs <- month.abb[as.numeric(gsub("0","",substr(LSMS_AGG$Group.1, 1, 2)))]
          }
          else if (time=="perhour_15") {
               typestr <- c("Per Hour 15 days")
               gtype <- c("p")
	           labs <- month.abb[as.numeric(gsub("0","",substr(LSMS_AGG$Group.1, 1, 2)))]
          }
          else if (time=="daily"){
               typestr <- c("Daily Mean")
               gtype <- c("p")
	           labs <- substr(LSMS_AGG$Group.1, 1, 4) 
          }
          else if (time=="halfmonth"){
               typestr <- c("Half Month")
               gtype <- c("b")
	           labs <- substr(LSMS_AGG$Group.1, 1, 4)
          }
          else if (time=="monthly"){           
               typestr <- c("Monthly Mean")
               gtype <- c("b")
               labs <- substr(LSMS_AGG$Group.1, 1, 4)
          }

          if(variable_names[1,V]=="fapar") {
            if (time != "perhour" ) {
              FLUXNET_AGG2 <- LSMS_FLUXNET_COMPUTE_GROUPING(FLUXNET_V2,TIME, time,'mean')
              FLUXNET_AGG3 <- LSMS_FLUXNET_COMPUTE_GROUPING(FLUXNET_V3,TIME, time,'mean')

              RMSE_seawifs <- LSMS_FLUXNET_EVALUATION("rmse", LSMS_AGG$x, FLUXNET_AGG$x)
              RMSE_modis <- LSMS_FLUXNET_EVALUATION("rmse", LSMS_AGG$x, FLUXNET_AGG2$x)
              RMSE_cyclopes <- LSMS_FLUXNET_EVALUATION("rmse", LSMS_AGG$x, FLUXNET_AGG3$x) 
              RR_seawifs <- LSMS_FLUXNET_EVALUATION("rr", LSMS_AGG$x, FLUXNET_AGG$x)
              RR_modis <- LSMS_FLUXNET_EVALUATION("rr", LSMS_AGG$x, FLUXNET_AGG2$x)
              RR_cyclopes <- LSMS_FLUXNET_EVALUATION("rr", LSMS_AGG$x, FLUXNET_AGG3$x)

              LSMS_FLUXNET_GRAPHIC_FAPAR(LSMS_AGG$x, FLUXNET_AGG$x, FLUXNET_AGG2$x, FLUXNET_AGG3$x, 
	         LSMS_AGG$Group.1, typestr, gtype,variable_info, RMSE_seawifs, RMSE_modis, RMSE_cyclopes, 
                 RR_seawifs, RR_modis, RR_cyclopes, NSE, PBIAS, RSR)
             
            }
          } else {
              LSMS_FLUXNET_GRAPHIC(LSMS_AGG$x, FLUXNET_AGG$x, LSMS_AGG$Group.1, 
                typestr, gtype, labs, variable_info, RMSE, NRMSE, RR, NSE, PBIAS, RSR)
              LSMS_FLUXNET_GRAPHIC_ERR(LSMS_AGG$x,FLUXNET_AGG$x,LSMS_AGG$Group.1, LSMS_SD, FLUXNET_SD, 
                variable_info, typestr)
          }     
      } else {
          annual_values[1,V]=mean(LSMS_V,na.rm=TRUE)
          annual_values[2,V]=mean(FLUXNET_V,na.rm=TRUE)
      }
    }
  }

  LSMS_V <- get.var.ncdf( modelfile, 'vegC')
  annual_values[1,nvar+1]=mean(LSMS_V,na.rm=TRUE)
  LSMS_V <- get.var.ncdf( modelfile, 'litterC')
  annual_values[1,nvar+2]=mean(LSMS_V,na.rm=TRUE)
  LSMS_V <- get.var.ncdf( modelfile, 'soilC')
  annual_values[1,nvar+3]=mean(LSMS_V,na.rm=TRUE)

} else {
  print('length of output and evaluation file do not match, verify simulation settings')
}

# Write Mean Values Into File
rownames=c(variable_names[4,],'vegC','litterC','soilC')
annual_values[3,1:nvar]=variable_units[3,1:nvar]
annual_values[3,nvar+(1:3)]='gC/m2'
write.table(t(annual_values),file="summary.annual.txt",col.names=c('LSMS','Observed','Units'),
  row.names=rownames) 

# Write Statistic Values Into File
for (si in 1:(length(aggregation_level)-1))
 write.table(t(eval_stats[,,si]),file=paste("statistics.", aggregation_level[si], ".txt", sep=''),
   col.names=toupper(stat_name), row.names=variable_names[4,]) 

# Statistic Graphic 
LSMS_FLUXNET_STATS_GRAPHIC(c("pbias","nse","rsr"), c("perhour","daily","halfmonth","monthly"))

# Taylor Diagram
LSMS_FLUXNET_TAYLOR(modelfile, datafile, "halfhourly")
LSMS_FLUXNET_TAYLOR(modelfile, datafile, "perhour")
LSMS_FLUXNET_TAYLOR(modelfile, datafile, "daily")
LSMS_FLUXNET_TAYLOR(modelfile, datafile, "halfmonth")
LSMS_FLUXNET_TAYLOR(modelfile, datafile, "monthly")
LSMS_FLUXNET_TAYLOR(modelfile, datafile, "annual")
