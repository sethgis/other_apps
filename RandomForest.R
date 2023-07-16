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
getwd()


sentinel <- brick("AUS_BANDS.tif")
extent(sentinel)
s1<-brick("S1_AUS_2_2.tif")

extent(sentinel)<-extent(s1)
extent(s1)
extent(elevation)<-extent(sentinel)
raster_list <- list(s1, sentinel)  # Add more rasters to this list if needed

# Find the minimum number of columns
min_cols <- min(sapply(raster_list, ncol))
min_cols









bands<-stack(s1, sentinel)
bands






samples <- st_read("AUTRALIA_MERGED_GEDI.shp")
samples
samples_2 <- na.omit(samples)

dt <- sentinel %>% 
  extract(y = samples_2) %>% 
  as.data.frame %>% 
  mutate(agdb = samples_2$agbd)
dt
as.data.frame(dt)
lc_samples.df2 <- na.omit(dt)
lc_samples.df2
#column.names(lc_samples.df2)


(predictors <- colnames(lc_samples.df2)[1:144])


(fo <- as.formula(paste("agdb ~", paste(predictors, collapse = "+"))))
fo

fit <- lda(fo, data = lc_samples.df2)

fit

pred <- predict(fit, newdata = lc_samples.df2)$agdb

summary(pred)

library("ranger")
fit <- ranger(fo, data = lc_samples.df2, importance="impurity")

fit
fit$variable.importance


pred <- predict(fit, data =lc_samples.df2 , type = "response")
pred



###predict for the whole raster using fitted model

classified <- predict(bands, fit, 
                      fun = function(model, ...) predict(model,   ...)$predictions)

classified

plot(classified)
##function that is called for error computation

ovAcc <- function(conmat)
{
  # number of total cases/samples
  n = sum(conmat)
  # number of correctly classified cases per class
  diag = diag(conmat)
  # Overall Accuracy
  OA = sum(diag) / n
  #
  # observed (true) cases per class
  rowsums = apply(conmat, 1, sum)
  p = rowsums / n
  # predicted cases per class
  colsums = apply(conmat, 2, sum)
  q = colsums / n
  expAccuracy = sum(p*q)
  kappa = (OA - expAccuracy) / (1 - expAccuracy)
  #
  # Producer accuracy
  PA <- diag / colsums
  # User accuracy
  UA <- diag / rowsums
  outAcc <- data.frame(producerAccuracy = PA, userAccuracy = UA)
  print(outAcc)
  #
  global_acc = data.frame(overallAccuracy=OA, overallKappa=kappa)
  print(global_acc)
}
#function to call the error computing matrix

conmat <- table(pred = classified, obs = lc_samples.df2$agbd)
cart_acc <- ovAcc(conmat)




writeRaster(classified, filename="RF_AUS_TIF2_s1+s2.tif", 
            format="GTiff", overwrite=TRUE)

classified$layer

plot(classified)

###CART for random forest regressor - orintingg the classification tree
library("rpart")
fit <- rpart(fo, data = lc_samples.df2)

#print the trained model
library(rpart.plot) 
rpart.plot(fit)

#References
##https://cran.r-project.org/web/packages/sperrorest/vignettes/spatial-modeling-use-case.html

##https://doi.org/10.5194/nhess-5-853-2005 -publication


##https://cran.r-project.org/web/packages/sperrorest/vignettes/spatial-modeling-use-case.html
##produce fruit tree crop in central chile














































set.seed(321)
# A stratified random split of the data



set.seed(321)
# A stratified random split of the data
idx_train <- createDataPartition(dt$agdb,
                                 p = 0.8, # percentage of data as training
                                 list = FALSE)
dt_train1 <- dt[idx_train]
dt_test1 <- dt[-idx_train]
dt_test
dt_train
table(dt_test)
table(dt_train)
as.data.frame(dt_train)


n_folds <- 10
set.seed(321)
folds <- createFolds(1:nrow(dt_train), k = n_folds)
# Set the seed at each resampling iteration. Useful when running CV in parallel.
seeds <- vector(mode = "list", length = n_folds + 1) # +1 for the final model
for(i in 1:n_folds) seeds[[i]] <- sample.int(1000, n_folds)
seeds[n_folds + 1] <- sample.int(1000, 1) 


ctrl <- trainControl(summaryFunction = multiClassSummary,
                     method = "cv",
                     number = n_folds,
                     search = "grid",
                     classProbs = TRUE, # not implemented for SVM; will just get a warning
                     savePredictions = TRUE,
                     index = folds,
                     seeds = seeds)

ctrl

cl <- makeCluster(3/4 * detectCores())
registerDoParallel(cl)
mat <- as.matrix(dt_train)
mat
model_rf <- caret::train(agdb ~ . , method = "rf", data = mat,
                         importance = TRUE, # passed to randomForest()
                         # run CV process in parallel;
                         # see https://stackoverflow.com/a/44774591/5193830
                         allowParallel = TRUE,
                         tuneGrid = data.frame(mtry = c(2, 3, 4, 5, 8)),
                         trControl = ctrl)
stopCluster(cl); remove(cl)
# Unregister the doParallel cluster so that we can use sequential operations
# if needed; details at https://stackoverflow.com/a/25110203/5193830
registerDoSEQ()
saveRDS(model_rf, file = "./cache/model_rf.rds")






















class(dt_train)
print(dt_train)
as.data.frame(dt_train)
nrow(dt_train)
colnames(dt_train)
view(dt_train)
table(dt_train$agdb)

table(dt_test$agdb)
n_folds <- 10
set.seed(321)
folds <- createFolds(1:nrow(dt_train), k = n_folds)
folds
# Set the seed at each resampling iteration. Useful when running CV in parallel.
seeds <- vector(mode = "list", length = n_folds + 1) # +1 for the final model
for(i in 1:n_folds) seeds[[i]] <- sample.int(1000, n_folds)
seeds[n_folds + 1] <- sample.int(1000, 1) # seed for the final model

ctrl <- trainControl(summaryFunction = multiClassSummary,
                     method = "cv",
                     number = n_folds,
                     search = "grid",
                     classProbs = TRUE, # not implemented for SVM; will just get a warning
                     savePredictions = TRUE,
                     index = folds,
                     seeds = seeds)


cl <- makeCluster(3/4 * detectCores())
registerDoParallel(cl)
ncol(dt_train)

predictors <- colnames(dt_train)[1:192]
predictors
(fo <- as.formula(paste("id_cls ~", paste(predictors, collapse = "+"))))
library(MASS)

class(dt_train)

column_numbers_to_remove <- c(14,15,16,30,31,32,46,47,48,62,63,64,78,79,80,94,95,96,110, 111, 112, 126, 127, 128, 142, 143, 144, 158, 159, 160, 174, 175, 176, 190, 191, 192)
dt_train <- dt_train[, -column_numbers_to_remove]
class(dt_train)
print(dt_train)
dt_train <- dt_train[,-c(14,15,16,30,31,32,46,47,48,62,63,64,78,79,80,94,95,96,110, 111, 112, 126, 127, 128, 142, 143, 144, 158, 159, 160, 174, 175, 176, 190, 191, 192)]

fit <- lda(fo, data = dt_train)

class(dt_train)

column_number <- 2  # Specify the column number you want to remove
#dt_train <- dt_train[, -14]
fit <- lda(fo, data = dt_train)
dt_train

#id_cls

mat <- as.matrix(dt_train)
mat
fo2 <- as.formula(paste("class ~", paste(mat, collapse = "+")))

fit <- lda(fo2, data = mat)





























model_rf <- caret::train(id_cls ~ . , method = "rf", data = mat,
                         importance = TRUE, # passed to randomForest()
                         # run CV process in parallel;
                         # see https://stackoverflow.com/a/44774591/5193830
                         allowParallel = TRUE,
                         tuneGrid = data.frame(mtry = c(2, 3, 4, 5, 8)),
                         trControl = ctrl)















coord <- st_coordinates(samples)
coord
samples$x <- coord[,1]
samples
samples$y <- coord[,2]
samples

samples = st_drop_geometry(samples)
class(samples)
samples
summary(samples)
samples_2 <- na.omit(samples)
samples_2
column.names(samples_2)

lc_samples.df = st_drop_geometry(samples)
class(lc_samples.df)
lc_samples.df2 <- na.omit(lc_samples.df)

(predictors <- colnames(samples_2)[15])
predictors
(fo <- as.formula(paste("class ~", paste(predictors, collapse = "+"))))
fit <- lda(fo, data = lc_samples.df2)
















