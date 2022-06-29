library(raster)
library(sp)
library(ggplot2)
library(scales)




r1<-raster('ndvi_2021.tif')
plot(r1)
r2<-raster('NDVI_2019-0000023296-0000000000.tif')
r3<-raster('NDVI_2019-0000000000-0000023296.tif')
r4<-raster('NDVI_2019-0000023296-0000023296.tif')




rc <- reclassify(r1, c(-Inf,0.1523,1, 0.1523,0.6524,2, 0.6524,0.7219,3, 
                       0.7219,0.9381,4, 0.9381,Inf,5))
plot(rc)
tifoptions <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")#
writeRaster(rc, "SC_reclass_r.tif",
            options = tifoptions, overwrite = TRUE)

r1.min = cellStats(r1, "min")
r1.max = cellStats(r1, "max")

r.scale <- ((r1 - r1.max) / (r1.min - r1.max))

plot(r.scale)
writeRaster(r.scale, "ndvi_2021_fuz_r_B.tif",
            options = tifoptions, overwrite = TRUE)




rmoz<- mosaic(r1,r2,r3,r4
              ,fun = mean,tolerance=0.05, overwrite = TRUE) 

tifoptions <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")#

writeRaster(rmoz, "ndvi_2019_compressed.tif",
            options = tifoptions, overwrite = TRUE)
plot(rmoz)



writeRaster(rmoz
            ,fun = mean, filename = 'ndvi_2020mos.tif',tolerance=0.05, overwrite = TRUE)













r5[r5==0] <- NA
R6 <- cover(r5)
ras<-raster('CE_2019_NEW.tif')
y <- reclassify(ras, cbind(-9999, 0),filename = 'R_CE_2019.tif')
plot(y)
q <- writeraster(y
                 ,fun = mean, filename = 'R_CE_2019.tif',tolerance=0.05, overwrite = TRUE)
