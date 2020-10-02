##Script for running fuzzy iAA using bootstrapping

library(sf)
library(raster)
library(data.table)
library(foreach)
library(tidyverse)
library(fasterize)
library(raster)
library(clhs)
library(ranger)
library(scales)
Rcpp::sourceCpp('FuzzySpread.cpp')

##load covariates
covDir <- "./InputData/Covariates/"
covs <- list.files(covDir)
ancDat <- stack(paste0(covDir,covs))

##training data
trDat <- st_read(dsn = "./InputData/", layer = "transect1_all_pts")
trDat <- trDat[grep("SBSmc2",trDat$mapunit1),]
trDat$id2 <- gsub("_[[:alpha:]].*","",trDat$id)

# ##for buffer
# trDat <- as.data.table(trDat)
# legend <- data.table(Unit = as.character(unique(trDat$mapunit1)))
# legend[,pID := seq_along(Unit)]
# trDat[legend, pID := i.pID, on = c(mapunit1 = "Unit")]
# trDat <- st_as_sf(trDat)
# 
# trRast <- rasterize(trDat,ancDat[[1]],field = "pID")
# 
# test <- foreach(unit = legend$pID, .combine = rbind) %do% {
#   tempR <- trRast
#   tempR[tempR != unit] <- NA
#   mat1 <- as.matrix(tempR)
#   
#   ## spread to adjacent cells
#   rSpread(mat1, diag = (unit-0.5), near = (unit-0.5)) ##updates mat1 in place
#   values(tempR) <- mat1
#   pts <- rasterToPoints(tempR,spatial = T) %>% st_as_sf()
#   pts
# }

test <- raster::extract(ancDat, trDat)
atts <- st_drop_geometry(trDat) %>% select ("mapunit1","mapunit2","id")
trAll <- cbind(atts, test) %>% as.data.table()
trAll[,id2 := gsub("_[[:alpha:]].*","",id)]
trID <- unique(trAll$id2) %>% as.character()

### read in matrix of fuzzy values
fuzzMat <- fread("SBSmc2_FuzzyClass.csv")
fMat <- data.table::melt(fuzzMat)
setnames(fMat, c("mapunit1","Pred","fVal"))

bsRes <- foreach(it = 1:100, .combine = rbind) %do% {
  testID <- sample(trID,1)
  testDat_sf <- trDat[trDat$id2 == testID,]
  testDat <- trAll[id2 == testID,]
  trainDat <- trAll[id2 != testID,]
  trainDat <- trainDat[is.na(mapunit2),]
  trainDat[,`:=`(mapunit2 = NULL,id = NULL,id2 = NULL)]
  numTr <- trainDat[,.(Num = .N), by = .(mapunit1)]
  numTr[,New := as.integer(rescale(log(Num), to = c(50,800)))]
  
  ###clhs to reduce tps within unit
  trainClean <- foreach(unit = numTr$mapunit1, .combine = rbind) %do% {
    if(numTr[mapunit1 == unit,(Num)] > numTr[mapunit1 == unit,(New)]){
      dat <- trainDat[mapunit1 == unit,]
      lhs <- clhs(dat[,!"mapunit1"], size = numTr[mapunit1 == unit,(New)], 
                  use.cpp = T, simple = F,iter = 5000)
      res <- as.data.table(lhs$sampled_data)
      res[,mapunit1 := unit]
    }else{
      dat <- trainDat[mapunit1 == unit,]
      dat
    }
  }
  trainClean <- unique(trainClean)
  trainClean[,mapunit1 := as.factor(as.character(mapunit1))]
  
  ##create model
  mod <- ranger(mapunit1 ~ ., data = trainClean, 
                num.trees = 501,importance = "impurity",
                splitrule = "extratrees",
                classification = T)
  
  temp <- predict(mod,data = testDat,predict.all = T)
  predMat <- as.data.table(temp$predictions)
  predMat[,pID := seq_along(1:nrow(predMat))]
  predMat[,mapunit1 := testDat$mapunit1]
  predMat <- data.table::melt(predMat, id.vars = c("pID","mapunit1"))
  predMat <- predMat[,.(Num = .N), by = .(mapunit1,pID,value)]
  predMat[,Prop := Num/501]
  leg2 <- data.table(UnitName = mod$forest$levels, ID = seq_along(mod$forest$levels))
  predMat[leg2, PredUnit := i.UnitName, on = c(value = "ID")]
  predMat[fMat, fVal := i.fVal, on = c("mapunit1",PredUnit = "Pred")]
  predMat[,TotProp := Prop * fVal]
  predMat <- predMat[complete.cases(predMat),]
  
  testDat$Pred <- predict(mod,data = testDat)$predictions
  testDat <- testDat[,.(mapunit1,mapunit2,Pred)]
  testDat[fMat, fVal1 := i.fVal, on = c("mapunit1","Pred")]
  testDat[fMat, fVal2 := i.fVal, on = c(mapunit2 = "mapunit1","Pred")]
  testDat[is.na(fVal2), fVal2 := 0]
  testDat[,fValMax := max(fVal1,fVal2), by = seq(1:nrow(testDat))]
  
  res <- foreach(unit = unique(testDat$mapunit1), .combine = rbind) %do% {
    sub <- testDat[mapunit1 == unit,]
    sub2 <- predMat[mapunit1 == unit,]
    out <- data.table(Unit = unit, Acc1 = sum(sub$fVal1)/nrow(sub),
                      Acc2 = sum(sub$fValMax)/nrow(sub),
                      Acc3 = sum(sub2$TotProp)/length(unique(sub2$pID)))
    out
  }
  res[,It := it]
  res
}





###old####
bsRes[,Unit := as.factor(as.character(Unit))]
boxplot(Acc2 ~ Unit, data = bsRes)
bsRes <- data.table::melt(bsRes, id.vars = c("Unit","It"))
bsRes[,variable := as.factor(variable)]

ggplot(bsRes, aes(x = Unit, y = value, fill = variable))+
  geom_boxplot()

map <- raster("./InputData/response.tif")
tps <- st_read("./InputData/transect1_all_pts_att.gpkg")
tps <- tps[,c("mapunit1")]
tps <- tps[grep("SBSmc2",tps$mapunit1),]
tps <- as.data.table(tps) %>% st_as_sf()
legend <- fread("./InputData/response_names.csv")
legend <- legend[complete.cases(legend),]
setnames(legend, c("Code","mapunit1"))

### read in matrix of fuzzy values
fuzzMat <- fread("SBSmc2_FuzzyClass.csv")
fMat <- data.table::melt(fuzzMat)
setnames(fMat, c("Unit1","Unit2","fVal"))

fMat[legend, U1Val := i.Code, on = c(Unit1 = "mapunit1")]
fMat[legend, U2Val := i.Code, on = c(Unit2 = "mapunit1")]
fMat[,`:=`(Unit1 = NULL, Unit2 = NULL)]

tps <- legend[tps, on = "mapunit1"]
tps <- st_as_sf(tps)
realRast <- rasterize(tps,map,field = "Code")
mapAll <- mask(map,realRast)
realRast <- crop(realRast, tps)

##
tpRast <- realRast
mapMask <- crop(mapAll, realRast)
tempStack <- stack(tpRast, mapMask)
temp <- values(tempStack)
temp <- temp[!is.na(temp[,1]),]
dat <- as.data.table(temp)
setnames(dat,c("Real","Predict"))
dat[fMat, FuzzProp := i.fVal, on = c(Real = "U1Val", Predict = "U2Val")]
dat <- dat[complete.cases(dat),]
##overall fuzzy accuracy
sum(dat$FuzzProp)/nrow(dat)

fuzzAA <- foreach(unit = unique(fMat$U1Val), .combine = rbind) %do% {
  cat("Processing", unit, "...\n")
  tpRast <- realRast
  tpRast[tpRast != unit] <- NA
  tpRast[!is.na(tpRast)] <- 1
  mat1 <- as.matrix(tpRast)
  
  ## spread to adjacent cells
  rSpread(mat1, diag = 0.3, near = 0.5) ##updates mat1 in place
  values(tpRast) <- mat1
  
  mapMask <- crop(mapAll, tpRast)
  tempStack <- stack(tpRast, mapMask)
  temp <- values(tempStack)
  temp <- temp[!is.na(temp[,1]),]
  dat <- as.data.table(temp)
  setnames(dat,c("Real","Predict"))
  
  ##fuzzy similarity
  fMat2 <- fMat[U1Val == unit,]
  dat[fMat2, FuzzProp := i.fVal, on = c(Predict = "U2Val")]
  dat <- dat[complete.cases(dat),]
  
  ### Real contains fuzzy values based on area: 1 if actual call
  ### FuzzProp contains fuzzy values from matrix based on similarity
  ### Predict is code name for predicted values of cell
  dat[,FuzzSum := Real*FuzzProp]
  dat2 <- dat[Real == 1,]
  dat3 <- dat[(Predict == unit | Real == 1),]
  data.table(Unit = unit, Acc2Fuzz = sum(dat$FuzzSum)/sum(dat$Real), 
             Acc1Fuzz = sum(dat2$FuzzSum)/sum(dat2$Real), Acc2Fuzzv2 = sum(dat3$FuzzSum)/sum(dat3$Real))
}

