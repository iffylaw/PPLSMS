require(RNetCDF)
require(ggplot2)
require(animation)

ncdat = read.nc(open.nc('wrfout_LH_201001-12_monavg_final.nc'))

drawit = function(i,ncdat)
{
    dfr = data.frame(Lon = as.vector(ncdat$XLONG),Lat = as.vector(ncdat$XLAT),LH = as.vector(ncdat$LH[,,i]))
    p = ggplot(aes(x=Lon,y=Lat,color=LH),data=dfr)
    p+geom_point()+
        scale_colour_gradientn(
            limits = c(0,150),
            colours = c("#0000FF","#00FF00","#FF0000"),
            guide = "colourbar"
        )+
        labs(title = paste('Latent Heat Flux in Qinghai-Tibet Plateau (2010-',i,', W/m^2)',sep=''))
}

saveGIF({
    for (i in 1:12) print(drawit(i,ncdat))
    },movie.name='NetCDFwithR.gif',interval = 1,ani.width=1000,ani.height=600)

