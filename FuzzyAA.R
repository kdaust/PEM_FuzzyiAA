library(sf)
library(raster)
library(data.table)
library(foreach)
library(tidyverse)
library(fasterize)

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

