#!/usr/bin/env Rscript
library(ncdf)

# Data Function (1) Get JSBACH Data-----------------------------------------------------------------------
JSBACH_DATA <- function(modelfile, V){
  JSBACH_V <- get.var.ncdf( modelfile, V )
  JSBACH_V[JSBACH_V<=-9999.]=NA
  return(JSBACH_V)
}

# Data Function (2) Get FLUXNET Data-----------------------------------------------------------------------
FLUXNET_DATA <- function(datafile, V){
  FLUXNET_V <- get.var.ncdf( datafile, V )
  FLUXNET_V[FLUXNET_V==-9999.]=NA
  if(!V=='Rn') {
     FLUXNET_V_FLAG <- get.var.ncdf( datafile, paste(V,'_flag',sep='') )
     FLUXNET_V[FLUXNET_V_FLAG>1]=NA
  }
  return(FLUXNET_V)   
}

datafile <- open.ncdf(paste("US-Ho1.1996-2004.obs.daily.nc", sep=""))
modelfile <- open.ncdf(paste("US-Ho1.rev270.daily.lsms.daily.nc", sep=""))
forcing <- open.ncdf(paste("US-Ho1.1996-2004.forcing.daily.nc", sep=""))
TIME <- get.var.ncdf( modelfile, "time")
TIME_F <- get.var.ncdf( datafile, "time" )

pdf("rel_hum.pdf")
timeformat <- c("%Y-%m-%d", "%Y-%m", "%m")
for (T in 1:3){
#timeformat <- "%Y-%m-%d"
time <- format(as.Date(as.numeric(TIME),origin="1996-01-01"), timeformat[T])
#which(time==)
if(length(TIME) == length(TIME_F)) {
  J_LE <- JSBACH_DATA(modelfile, "latent_heat_flx")
  J_LE <- J_LE*(-1)
  F_LE <- FLUXNET_DATA(datafile, "LE")

  J_EVA <- JSBACH_DATA(modelfile, "evaporation")*(-1000)
  F_RHUM <- FLUXNET_DATA(forcing, "rel_humidity")
  F_SHUM <- JSBACH_DATA(modelfile, "soil_moisture")

  J_LET <- data.frame(time, J_LE)
  J_LE <- aggregate(J_LET$J_LE, list(J_LET$time), "mean", na.rm=TRUE)

  F_LET <- data.frame(time, F_LE)
  F_LE <- aggregate(F_LET$F_LE, list(F_LET$time), "mean", na.rm=TRUE)

  J_EVA_T <- data.frame(time, J_EVA)
  J_EVA <- aggregate(J_EVA_T$J_EVA, list(J_EVA_T$time), "mean", na.rm=TRUE)

  F_RHUM_T <- data.frame(time, F_RHUM)
  F_RHUM <- aggregate(F_RHUM_T$F_RHUM, list(F_RHUM_T$time), "mean", na.rm=TRUE) 

  F_SHUM_T <- data.frame(time, F_SHUM)
  F_SHUM <- aggregate(F_SHUM_T$F_SHUM, list(F_SHUM_T$time), "mean", na.rm=TRUE)

  time <- unique(time)
  opar<-par(pin=c(6.0, 2.0), cex=0.6)
  #plot(1:length(time), J_LE$x, ylim=c(min(J_LE$x, F_LE$x, na.rm=TRUE), max(J_LE$x, F_LE$x, na.rm=TRUE)),col="red", type="n", xlab="Month", ylab="Latent Heat (w/m^2)", main="")
  #labs <- month.abb
  #axis(1, at=1:12, lab=labs,las=1, cex=0.6)
  #axis(2,las=1)
  #box()
  #lines(1:length(time), F_LE$x, col="black", type="p", cex=0.6)
  #lines(1:length(time), J_LE$x, col="red", type="p", cex=0.6)
  #grid()
  #polygon(c(1:length(time),length(time):1),c(F_LE$x, J_LE$x[length(time):1]), col="grey", lty=0)

  split.screen(c(2,1))
  screen(1)  
  plot(1:length(time), J_EVA$x, type='o', col='green', lwd=0.5, cex=0.6, ylab="Evaporation(g/m^2s)")
  par(new=T)
  plot(1:length(time), F_RHUM$x, type='o', col='blue1', lwd=0.5, ann=F, axes=F, cex=0.6, ylab="Relative Humidity (%)")
  #lines(1:length(time), F_SHUM$x, type='o', col='grey', lwd=0.5, ann=F, axes=F, cex=0.6, ylab="Soil Moisture (m)")
  axis(4, col="blue1")
  legend("bottom", legend=c("Evaporation","Relative Humidity"), col=c("green","blue1"), pch=1)

  screen(2)
  plot(1:length(time), F_RHUM$x, type='o', col='blue1', ylim=c(0,100), ann=F, axes=F, lwd=0.5,  cex=0.6, ylab="Relative Humidity (%)")
  par(new=T)
  plot(1:length(time), F_SHUM$x, type='o', col='green', lwd=0.5,  cex=0.6, ylab="Soil Moisture (m)")
  axis(4, col="green") 
  legend("bottom", legend=c("Relative Humidity", "Soil Moisture"), col=c("blue1", "green"), pch=1)

  #par(new=T)
  #plot(1:length(time), F_SHUM$x, type='o', col='green', lwd=0.5, ann=F, axes=F, cex=0.6, ylab="Relative Humidity")
  close.screen(1:2)

  int <- sprintf("%.4f", coef(line(F_RHUM$x, J_EVA$x),use.na=false)[1])
  slope <- sprintf("%.4f", coef(line(F_RHUM$x, J_EVA$x),use.na=false)[2])
  r  <- sprintf("%.4f",cor(F_RHUM$x, J_EVA$x,use="na.or.complete"))
  rr <- sprintf("%.4f",summary(lm(J_EVA$x ~ F_RHUM$x, na.action=na.omit))[c("r.squared")])

  int1 <- sprintf("%.4f", coef(line(F_RHUM$x, F_SHUM$x),use.na=false)[1])
  slope1 <- sprintf("%.4f", coef(line(F_RHUM$x, F_SHUM$x),use.na=false)[2])
  r1 <- sprintf("%.4f",cor(F_RHUM$x, F_SHUM$x,use="na.or.complete"))
  rr1 <- sprintf("%.4f",summary(lm(F_SHUM$x ~ F_RHUM$x, na.action=na.omit))[c("r.squared")])

  opar<-par(pin=c(3.0, 3.0))
  #split.screen(c(2,1))
  #screen(1)
  plot(F_RHUM$x, J_EVA$x, xlab="Relative Humidity (%)", ylab="Evaporation(g/m^2s)", cex=1)
  abline(lm(J_EVA$x ~ F_RHUM$x), col="red")
  #abline(0,1)
  p1 <- t.test(F_RHUM$x, J_EVA$x)$p.value
  legend("bottomleft", legend=c(paste("Intercept=",int),paste("Slope=",slope),paste("R=",r), as.expression(bquote(R^2 == .(rr)))), col=c("red"), cex=0.6, bty="n")
  #print(t.test(F_RHUM$x, J_EVA$x))

  #screen(2)
  opar<-par(pin=c(3.0, 3.0))
  plot(F_RHUM$x, F_SHUM$x, xlab="Relative Humidity (%)", ylab="Soil Moisture(m)", cex=1)
  abline(lm(F_SHUM$x ~ F_RHUM$x), col="red")
  #abline(0,1)
  p2 <- t.test(F_RHUM$x, F_SHUM$x)$p.value
  legend("bottomleft", legend=c(paste("Intercept=",int1),paste("Slope=",slope1),paste("R=",r1),as.expression(bquote(R^2 == .(rr1)))), col=c("red"), cex=0.6, bty="n")
  #print(t.test(F_RHUM$x, F_SHUM$x))
  close.screen(1:2) 
}
}
