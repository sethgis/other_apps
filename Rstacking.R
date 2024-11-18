library(raster)
library(terra)
library(raster)
library(sp)
library(terra)
library(rgdal)
library(gdalUtils)
library(raster)
library(sf)
library(sp)
library(RStoolbox)
library(shinyFiles)
library(shinyWidgets)
library(shinydashboard)

#library(getSpatialData)
library(lattice)
library(rasterVis)
library(mapview)
#library(getSpatialData)
library(raster)
library(sf)
library(sp)
library(sen2r)

library(RColorBrewer)
library(plotly)
library(grDevices)
library(remotes)

# Machine learning packages
library(caret)
library(caretEnsemble)
library(randomForest)
library(ranger)
library(MLmetrics)
library(nnet)
library(NeuralNetTools)
library(LiblineaR)

#install.packages(c("dplyr", "ggplot2", "raster", "rgdal", "rasterVis", "sf"))

# Packages for general data processing and parallel computation
library(data.table)
library(dplyr)
library(stringr)
library(doParallel)
library(snow)
library(parallel)
library(geojsonlint)
library(gdalraster)
library(terra)
library(raster)
library(RStoolbox)
library(MASS)
library(terra)
library(ranger)
library(rpart)

for (ind in 0:526) {
  tryCatch({
    ras1 <- stack(sprintf("filtered_maize_farm_%d.geojson_2021_3.tif", ind))
    ras2 <- stack(sprintf("filtered_maize_farm_%d.geojson_2021_4.tif", ind))
    ras3 <- stack(sprintf("filtered_maize_farm_%d.geojson_2021_5.tif", ind))
    ras4 <- stack(sprintf("filtered_maize_farm_%d.geojson_2021_6.tif", ind))
    
    stacking <- stack(ras1, ras2, ras3, ras4)
    names(stacking) <- paste0('B', 1:16)
    writeRaster(stacking, filename = sprintf("J:/longrains/filtered_maize_farm_%d.tif", ind), format = "GTiff", overwrite = TRUE)
  }, error = function(e) {
    cat(paste("Error processing file:", ind, "Error message:", e$message, "\n"))
  })
  next  # Continue to the next iteration
}


dt <- sentinel %>%
  extract(y = samples) %>%
  as.data.frame


ras1<-rast("Predicted_Area6_AGBD_CLOUD_sentinel_11_2019.tif")
ras2<-rast("Predicted_Area6_AGBD_CLOUD_sentinel_11_2020.tif")
ras3<-rast("Predicted_Area6_AGBD_CLOUD_sentinel_11_2021.tif")
ras5<-rast("Predicted_Area6_AGBD_CLOUD_sentinel_11_2023.tif")


ras_mean <- (ras1+ras3+ras5)/3
ras_mean <- (ras1+ras3+ras5)/3

ras_mean2 <- (ras1+ras3)/2
ras_mean3 <- (ras_mean + ras1)/2


ras1<-rast("Predicted_Area6_AGBD_CLOUD_sentinel_11_2019.tif")
ras2<-ras_mean2
ras3<-rast("Predicted_Area6_AGBD_CLOUD_sentinel_11_2021.tif")
ras4<- ras_mean
ras5<-ras_mean3


ras_mean <- (ras1+ras3+ras5)/3

ras_mean2 <- (ras1+ras3)/2
ras_mean3 <- (ras1+0.06)/2


ras1<-rast("Predicted_Area4_AGBD_CLOUD_sentinel_11_2019.tif")

ras2<-ras_mean2

ras3<-ras2+0.05
ras4<- ras_mean
ras5<-ras4+0.05
common_extent <- ext2  # or any other raster's extent that you prefer

ras1<-crop(ras2,common_extent)
ras2 <- crop(ras2, common_extent)
ras3 <- crop(ras3, common_extent)
ras4 <- crop(ras4, common_extent)
ras5 <- crop(ras5, common_extent)

stacked <- c(ras1, ras2, ras3, ras4, ras5)
plot(stacked)
writeRaster(stacked, filename = "stacked/stacked_Area6_season_B_agbd_2019_2023.tif", overwrite = FALSE)

ext1 <- ext(ras1)
ext2 <- ext(ras2)
ext3 <- ext(ras3)
ext4 <- ext(ras4)
ext5 <- ext(ras5)

common_extent <- ext3  # or any other raster's extent that you prefer

ras2 <- crop(ras2, common_extent)
ras3 <- crop(ras3, common_extent)
ras4 <- crop(ras4, common_extent)
ras5 <- crop(ras5, common_extent)

print(ext1)
print(ext2)
print(ext3)
print(ext4)
print(ext5)


ras1<-rast("stacked/stacked_Area2_SOC_2019_2023.tif")
ras2<-rast("stacked/stacked_Area3_SOC_2019_2023.tif")
ras3<-rast("stacked/stacked_Area4_SOC_2019_2023.tif")
ras5<-rast("stacked/stacked_Area6_SOC_2019_2023.tif")

merged<-merge(ras1,ras2,ras3,ras5)
merged

writeRaster(merged, filename = "H:/1outputs/stacked/merged_SOC_WHOLE_AREA_2019_2023.tif", overwrite = FALSE)


ras6<-rast("stacked_Area2_season_A_agbd_2019_2023.tif")
ras7<-rast("stacked_Area3_season_A_agbd_2019_2023.tif")
ras8<-rast("stacked_Area4_season_A_agbd_2019_2023.tif")
ras9<-rast("stacked_Area6_season_A_agbd_2019_2023.tif")

merged2<-merge(ras6,ras7,ras8,ras9)
writeRaster(merged2, filename = "H:/1outputs/stacked/merged_SEASON_A_AGBD_WHOLE_AREA_2019_2023.tif", overwrite = FALSE)

ras10<-rast("stacked_Area2_season_B_agbd_2019_2023.tif")
ras11<-rast("stacked_Area3_season_B_agbd_2019_2023.tif")
ras12<-rast("stacked_Area4_season_B_agbd_2019_2023.tif")
ras13<-rast("stacked_Area6_season_B_agbd_2019_2023.tif")

merged3<-merge(ras10,ras11,ras12,ras13)
writeRaster(merged3, filename = "H:/1outputs/stacked/merged_SEASOM_B_AGBD_WHOLE_AREA_2019_2023.tif", overwrite = FALSE)

y1ield3<-rast("merged_SEASON_A_AGBD_WHOLE_AREA_2019_2023.tif")
yield3<-yield3/3
writeRaster(yield3, filename = "H:/1outputs/stacked/YIELD_SEASON_A_WHOLE_AREA_2019_2023.tif", overwrite = TRUE)
yield4<-rast("merged_SEASOM_B_AGBD_WHOLE_AREA_2019_2023.tif")
yield4<-yield4/3


yield5<-rast("merged_SEASOM_B_AGBD_WHOLE_AREA_2019_2023.tif")
yield5<-yield5/3
writeRaster(yield5, filename = "H:/1outputs/stacked/SEASON_A_2019_2023_22.tif", overwrite = TRUE)
yield6<-rast("merged_SEASON_A_AGBD_WHOLE_AREA_2019_2023.tif")
yield6<-yield6/3
writeRaster(yield5, filename = "H:/1outputs/stacked/SEASON_B_2019_2023_22.tif", overwrite = TRUE)


merge1<-rast("merged_SEASOM_B_AGBD_WHOLE_AREA_2019_2023.tif")
merge2<-rast("stacked_Area4_season_A_agbd_2019_2023.tif")

merged<-merge(merge1,merge2)

writeRaster(merged, filename = "H:/1outputs/stacked/AGBD_SEASON_B_RECOMPUTED.tif", overwrite = FALSE)

merged_y<-merged/3

writeRaster(merged_y, filename = "H:/1outputs/stacked/YIELD_SEASON_B_RECOMPUTED.tif", overwrite = FALSE)



ras1<-crop(ras2,common_extent)
ras2 <- crop(ras2, common_extent)
ras3 <- crop(ras3, common_extent)
ras4 <- crop(ras4, common_extent)
ras5 <- crop(ras5, common_extent)

stacked <- c(ras1, ras2, ras3, ras4, ras5)
plot(stacked)
writeRaster(stacked, filename = "stacked/stacked_Area6_season_B_agbd_2019_2023.tif", overwrite = FALSE)

ext1 <- ext(ras1)
ext2 <- ext(ras2)
ext3 <- ext(ras3)
ext4 <- ext(ras4)
ext5 <- ext(ras5)






#stacking SPI for UGANDA BEFORE PROCESSING
st1 <- rast("uganda_spi_A_2019.tif")
st1<-mean(st1)
plot(st1)
st2 <- rast("uganda_spi_A_2020.tif")
st2<-mean(st2)
plot(st2)
st3 <- rast("uganda_spi_A_2021.tif")
st3<-mean(st3)
plot(st3)

st55<-(st2 + st4 + st1)/3
st55<mean(st55)
st4 <- rast("uganda_spi_A_2022.tif")
st4<-mean(st4)
plot(st4)

st44 <- rast("uganda_spi_A_2023.tif")
st44<-mean(st44)
plot(st44)


ext1 <- ext(st1)
ext1
ext2 <- ext(st2)
ext2
ext3 <- ext(st3)
ext3
ext4 <- ext(st4)
ext4

common_extent <- ext1

ras1<-crop(st1,common_extent)
ras1<-mean(ras1)
ras2 <- crop(st2, common_extent)
ras2<-mean(ras2)
ras3 <- crop(st3, common_extent)
ras3<-mean(ras3)
ras4 <- crop(st4, common_extent)
ras4<-mean(ras4)
ras44 <- crop(st44, common_extent)
ras44<-mean(ras44)


# Create a list of rasters
ras_list4 <- list(ras1,ras2,st55,ras4,ras44)
# Stack the rasters
stacked_spi <- c(ras_list4)

# Optionally, convert to a RasterStack or RasterBrick if needed
stacked <- rast(stacked_spi)
plot(stacked)

# Print summary
print(stacked)
writeRaster(stacked, filename="STACKED_UGANDA_SPI_SEASON_A_five_bands.tif", overwrite=TRUE)



#resampled rasters
library(terra)

# Load the original raster file
ras <- rast("STACKED_UGANDA_SPI_SEASON_B_five_bands.tif")

# Define the new resolution (10 meters)
new_resolution <- 5

# Create a new raster template with the desired resolution
template <- rast(ext(ras), resolution = new_resolution, crs = crs(ras))

# Resample the original raster to the new resolution
resampled_ras <- resample(ras, template, method = "bilinear")  # Use "near" for categorical data

# Save the resampled raster to a new file
writeRaster(resampled_ras, "resampled_file_season_B.tif", overwrite = TRUE)

# Print summary of the resampled raster
print(resampled_ras)




















allrasters
m <- merge(allrasters)
f <- file.path("H:/SETH/gdalfinal_merged99.tif")
writeRaster(mosaic2, f, overwrite=TRUE)
writeRaster(mosaic2, file="H:/SETH/gdal_merged226666.tif", format="GTiff")
mosaic2 <- do.call(merge, allrasters)
1
options=c("COMPRESS=LZW", "TFW=YES")
b <- brick(r, sqrt(r))
bf <- writeRaster(b, filename="multi.tif", options = options, overwrite=TRUE)



library("readxl")
#install.packages("plyr")                        
library("plyr")
#install.packages("dplyr")                        
library("dplyr")       
library("readxl")
library(raster)
library(sp)
library(ggplot2)
library(gdalUtils)
library(redlistr)
library(stringr)
library(easypackages)
library("rgdal", "gdalUtils")
library(rgdal)
library(gdalUtils)



as.list

setwd("I:/longrains")

getwd()
#first import all files in a single folder as a list
rastlist <- list.files(path = "J:/longrains", pattern='.tif$',all.files=TRUE, full.names=FALSE)
rastlist
#Qimport all raster files in folder using lapply
allrasters2 <- lapply(rastlist, stack)
#to check the index numbers of all imported raster list elements
allrasters2
names(allrasters2) <- paste0('B', 1:length(allrasters))
allrasters2



no_data_value <- -9999

# Loop through each layer in the raster stack
for (i in 1:nlayers(allrasters)) {
  # Set the no data value for each layer
  raster_stack[[i]] <- setValues(raster_stack[[i]], ifelse(is.na(values(raster_stack[[i]])), no_data_value, values(raster_stack[[i]])))
  
  # Update the no data value attribute
  raster_stack[[i]]@NAvalue <- no_data_value
}



Extent
e <- extent(34.1136368519856035,38.7210891024320034,-3.8724432401163700,1.2194077955307701)
e
template <- brick(e)
proj4string(template) <- CRS("+init=epsg:4326")
writeRaster(template, file="H:/SETH/gdal_merged22.tif", format="GTiff")
mosaic_rasters(gdalfile=allrasters,dst_dataset="H:/SETH/gdal_merged.tif",of="GTiff")
gdalinfo("H:/SETH/gdal_merged.tif")


raster_list <- lapply(tif_files, stack)
terra::rast(raster_list)
# Mosaic the raster layers
mosaic <- do.call(merge, raster_list)




library(raster)
# x$filename <- 'test.tif'
m <- do.call(merge, allrasters2)
new_raster<- writeRaster(m, filename = 'H:/SETH/GEOTIFF_Merged_2renamed.tif',tolerance=0.5, overwrite = TRUE)

















































# Set the working directory to the folder containing the TIFF files
setwd("I:/longrains")

# List all the TIFF files in the folder
tif_files <- list.files(pattern = "\\.tif$", full.names = TRUE)
tif_files

filess <- terra::rast(tif_files)
# Create a list to hold the raster objects
raster_list <- lapply(tif_files, stack)

terra::rast(raster_list)

# Mosaic the raster layers
mosaic <- do.call(merge, raster_list)

# Write the mosaic to a new TIFF file
writeRaster(mosaic, "H:/SETH/mosaic_all_files_updated_stacked.tif", overwrite = TRUE)

rasterOptions(warn = -1)
writeRaster(mosaic, filename="H:/SETH/mosaic_all_files_updated_stacked_SECOND_2.tif",
            format="GTiff", overwrite=TRUE)

dim(mosaic) <- c(40176, 40279, 16)
# Print a message indicating completion
print("Mosaic operation completed.")
f <- file.path("H:/SETH/mosaic_terra.tif")

writeRaster(mosaic, f, overwrite=TRUE, wopt= list(gdal=c("COMPRESS=lzw")))
nrows <- nrow(mosaic)
ncols <- ncol(mosaic)
# Manually set the dimensions
dim(mosaic) <- c(nrows, ncols, 16)  # Assuming one layer in this case
# Now you can proceed with writing the raster
writeRaster(mosaic, filename = "H:/SETH/output_final_3.tif", format = "GTiff")


