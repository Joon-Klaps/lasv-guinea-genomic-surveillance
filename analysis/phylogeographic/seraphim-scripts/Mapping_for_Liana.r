library(raster)
library(rgdal)

e_MRU = extent(-15, -4, 4, 13.04)

background = raster("Environmental_files/LASV_rasters/Elevation_MRU_0.008.asc")
background = crop(background, e_MRU); background[!is.na(background[])] = 1; r = background
background_cols = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
borders = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="International_borders"), e_MRU)
lakes = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="Natural_Earth_lakes"), e_MRU)
rivers = crop(readOGR(dsn="Environmental_files/LASV_rivers/", layer="FAO_GeoNetwork"), e_MRU)
rivers@data[which(rivers@data[,"Strahler"]<3),"Strahler"] = 0
riverCols = rep("#AFBED5",dim(rivers@data)[1]); riverCols[which(rivers@data[,"Strahler"]<3)] = NA

pdf(paste0("NEW_MAP_FIGURE.pdf"), width=12, height=10.2)
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(1.5,2.7,1.0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
plot(background, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col="white", colNA="gray90")
plot(lakes, add=T, lwd=0.7, col=unique(riverCols)[!is.na(unique(riverCols))], border=NA)
plot(rivers, add=T, lwd=rivers@data[,"Strahler"]/3.5, col=riverCols)
plot(borders, add=T, lwd=1, col="gray60", lty=1)
rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.2)
axis(1,c(ceiling(xmin(background)), floor(xmax(background))), pos=ymin(background), mgp=c(0,0.7,0), 
	 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
axis(2, c(ceiling(ymin(background)), floor(ymax(background))), pos=xmin(background), mgp=c(0,0.8,0),
	 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
dev.off()
