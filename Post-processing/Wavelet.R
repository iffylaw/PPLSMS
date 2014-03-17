library(ncdf)
library(dplR)
library(wmtsa)

# Compute the Global Wavelet Spectrum
# Power Numeric. The squared power.
power <- NA
GWS <- function(cwt){  
  for (i in 1:ncol(cwt)) {
    power[i] <- sum(cwt[,i])/nrow(cwt)
  }
  return(power)
}  

js <- open.ncdf("lsms.daily.nc")
fl <- open.ncdf("obs.daily.nc")

js_time <- get.var.ncdf(js,"time")
fl_time <- get.var.ncdf(fl,"time")

#opar<- par(mfrow = c(2, 1), mar = c(0,4.2,0,1))

js_data <- get.var.ncdf(js,"sensible_heat_flx")
js_data <- -js_data

fl_data <- get.var.ncdf(fl,"H_f")

js_morlet <- morlet(y1=js_data, x1=seq_along(js_data), dj=0.25, siglvl=0.95)
fl_morlet <- morlet(y1=fl_data, x1=seq_along(fl_data), dj=0.25, siglvl=0.95)

js_gws <- GWS(js_morlet$Power)
fl_gws <- GWS(fl_morlet$Power)

# Compute the Global Wavelet Spectrum from JSBACH divided by the corresponding spectrum from FLUXNET
jf_gws <- js_gws/fl_gws

# Use log function, JSBACH > FLUXNET, log(JSBACH/FLUXNET) > 0 ; ELSE < 0
plot(1:length(jf_gws), log2(jf_gws),xlab="Scale (days)", ylim=c(-2,2.5), ylab="log(Global Wavelet Spectrum)", type="l", main="Wavelet Analysis for DK-Sor", lwd=4, axes=F)
abline(h=0,col="grey", lty=1)
axis(1, at=1:length(js_morlet$Scale), lab=js_morlet$Scale,las=1)
axis(2,las=1)
box() 

#
js_data <- get.var.ncdf(js,"latent_heat_flx")
js_data <- -js_data

fl_data <- get.var.ncdf(fl,"LE_f")

js_morlet <- morlet(y1=js_data, x1=seq_along(js_data), dj=0.25, siglvl=0.95)
fl_morlet <- morlet(y1=fl_data, x1=seq_along(fl_data), dj=0.25, siglvl=0.95)

js_gws <- GWS(js_morlet$Power)
fl_gws <- GWS(fl_morlet$Power)

# Compute the Global Wavelet Spectrum from JSBACH divided by the corresponding spectrum from FLUXNET
jf_gws <- js_gws/fl_gws


# Use log function, JSBACH > FLUXNET, log(JSBACH/FLUXNET) > 0 ; ELSE < 0
lines(1:length(jf_gws), log2(jf_gws),xlab="Scale (days)", ylab="log(Global Wavelet Spectrum)", type="l", lwd=4, col="red")
#abline(h=0,col="red")
#axis(1, at=1:length(js_morlet$Scale), lab=js_morlet$Scale,las=1)
#axis(2,las=1)
#box() 
legend("topleft", legend=c("Sensible Heat", "Latent Heat"), col=c("black","red"), lty=1, lwd=2, bty="n")
