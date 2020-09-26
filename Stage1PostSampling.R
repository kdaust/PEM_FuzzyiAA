###Post-Processing PEM Transects
###Kiri Daust, Colin Chisholm, Aug 2019

library(sf)
library(lwgeom)
library(raster)
library(sp)
library(gdistance)
library(stars)
library(plyr)
library(dplyr)
library(foreach)
library(doSNOW)
library(stars)
library(iterators)
library(LearnGeom)
library(tmap)
library(viridis)
library(shinyjs)
library(formatR)
library(tidyr)
library(fasterize)


###read data
setwd("C:/Users/Kiri Daust/Desktop/Deception_SamplingPlan/Stage1Analysis")
trans <- "SBSmc2_5.1_21_cLHS" ###name of output zip file
setwd("TransectData")
unzip(paste(trans,".zip", sep = ""), exdir = trans)###unzip
t1 <- st_read(dsn = trans, layer = "SBSmc2new") %>% ##read transect data
  st_transform(3005) %>%
  st_zm()

###import original TransectPairs cLHS file to get ID
LHSPoints <- st_read(dsn = "../LHS Points", layer = "SBSmc2_TransectPairs4_2019-07-25") %>%
  st_transform(3005)
LHSPoints <- st_buffer(LHSPoints, dist = 10)
t1 <- st_join(t1, LHSPoints, join = st_intersects)
t1 <- t1[!is.na(t1$ID),]
allTrans <- t1


###Loop through each polygon, convert to raster and merge
setwd("../")
modKey <- read.csv("SSMap/response_key.csv")[,-1]###SS Unit - model key crosswalk
rTemp <- raster("2.5m/MRVBF.tif")###for template

rastAll <- foreach(trs = unique(allTrans$ID), .combine = function(x,y){raster::merge(x,y)}) %do% {
  t1 <- allTrans[allTrans$ID == trs,]
  t1$Order <- gsub("\\D+","", t1$name)
  t1 <- distinct(t1, .keep_all = T)
  t1$Order <- as.numeric(t1$Order)
  t1 <- t1[order(t1$Order),]
  rownames(t1) <- NULL
  t1 <- t1[!is.na(t1$X2MapUnit),]
  
  calls <- t1[,c("name","X2MapUnit")] %>% st_drop_geometry()
  
  lines <- cbind(t1, st_coordinates(t1)) %>% as_tibble() %>% ###connect points
    mutate(Xend = lead(X),
           Yend = lead(Y)) %>%  # collect the coordinates of the next point
    filter(!is.na(Yend)) %>%  # drops the last row (start point with no end)
    unite(Start, X, Y, sep = " ") %>%
    unite(End, Xend, Yend, sep = " ") %>%
    gather(key = "Start_End", value = coords, Start, End) %>% # converts to long table format
    dplyr::select(-geometry) %>%  # old point geometry is dropped
    separate(coords, into = c("X", "Y"), sep = " ") %>%       # coordinate pairs are seperated
    st_as_sf(coords = c("X", "Y"))  %>% # geometry is created for each point (start and end)
    group_by(name) %>% # with the point ID created above it is now possible to cast multiple points to a single geometry
    summarise() %>%
    st_cast("LINESTRING") %>%
    st_set_crs(st_crs(t1))
  
  lines <- merge(lines, calls, by = "name", all.x = TRUE)
  ##st_write(lines, dsn = "TestTransect",layer = "lines", driver = "ESRI Shapefile")
  
  lines <- lines %>%
    group_by(X2MapUnit) %>%
    summarise()
  
  lBuff <- st_buffer(lines, dist = 2.5, endCapStyle = "FLAT", joinStyle = "MITRE") ##buffer
  lBuff$X2MapUnit <- gsub("(.*[[:digit:]]).*","\\1",lBuff$X2MapUnit)
  plot(lBuff)
  
  lBuff$X2MapUnit <- gsub("/","_",lBuff$X2MapUnit)
  lBuff <- merge(lBuff, modKey, by.x = "X2MapUnit", by.y = "category", all.x = T)
  lBuff <- st_cast(lBuff, "MULTIPOLYGON")
  rast <- fasterize(lBuff, rTemp, field = "ID") ###rasterize
  rast <- crop(rast, lBuff)
}

#add other plot
# rastAll <- merge(rastAll, rast)
writeRaster(rastAll, "AllTransects.tif", format = "GTiff")


###################################################################
####Accuracy####################
#################################################################

setwd("C:/Users/Kiri Daust/Desktop/Deception_SamplingPlan/Stage1Analysis")
##Overal Accuracy 

mod <- raster("SSMap/map.tif") ##predicted map
rast <- rastAll
test <- crop(mod, rast)
plot(test)
plot(rast, add = T)
rDiff <- rast - test ###difference between map and data
plot(rDiff)
length(rDiff[rDiff == 0])/length(rDiff[!is.na(rDiff)]) ###total prop same

########confusion matrix###############
modVals <- mask(test, rast)

loopMerge <- function(x,y){
  merge(x,y,by = "Model", all.x = T, all.y = T)
}
confMat <- foreach(val = unique(rast[!is.na(rast)]), .combine = loopMerge) %do% {
 tTemp <- as.data.frame(table(modVals[rast == val]))
 colnames(tTemp) <- c("Model",val)
 #tTemp[,2] <- (tTemp[,2]/sum(tTemp[,2]))*100
 tTemp
}

temp <- rast
temp2 <- mask(modVals, temp)

temp <- as.character(modKey$category)
names(temp) <- modKey$ID
confMat$Model <- plyr::revalue(as.character(confMat$Model), temp) ##convert to SS names - ignore warning
rownames(confMat) <- confMat$Model
confMat <- confMat[,-1]
colnames(confMat) <- plyr::revalue(colnames(confMat), temp)
confMat <- confMat[,order(colnames(confMat))]
confMat <- confMat[,-1]
confMat <- confMat[order(rownames(confMat)),]
confMat


##############################################################################
###Extract training points
#########################################################################
###2018 Data 
setwd("C:/Users/Kiri Daust/Desktop/cLHSData")
pnts <- st_read(dsn = "Field Points", layer = "DeceptionAllFieldPts2018")
pnts <- pnts[,c("Site.Serie","MAP_LABEL")]
pnts <- pnts[grep("SBSmc2",pnts$Site.Serie),]
pnts$Site.Serie <- gsub("/","_",pnts$Site.Serie)
pnts <- pnts[pnts$Site.Serie %in% modKey$category,]
pnts <- merge(pnts, modKey, by.x = "Site.Serie", by.y = "category")
pnts <- st_transform(pnts, 3005)
pntsBuff <- st_buffer(pnts, dist = 2.5)
library(fasterize)
library(randomForest)
rPnts <- fasterize(pntsBuff, rTemp, field = "ID")
rAll <- raster::merge(rPnts, rast) ###merge 2018 points and 2019 transects
###or###3
rAll <- rast

setwd("C:/Users/Kiri Daust/Desktop/Deception_SamplingPlan/Stage1Analysis/2.5m")
rastList <- list.files(".")
numAnc <- length(rastList) ## read rasters of ancillary data
ancDat <- stack()
for(i in 1:numAnc){
  temp <- raster(rastList[i])
  ancDat <- stack(ancDat,temp)
}
ancDat <- stack(ancDat, rAll)
ancDat <- crop(ancDat, rAll)
train1 <- as.data.frame(ancDat) ###this can take a minute
train1 <- train1[!is.na(train1$layer),]
train1 <- train1[complete.cases(train1),]
train1$layer <- as.factor(train1$layer)
mod1 <- randomForest(layer ~ ., data=train1, nodesize = 5, do.trace = 10, ###random forest model with all layers
                         ntree=101, na.action=na.fail, importance=TRUE, proximity=FALSE)
imp <- as.data.frame(mod1[["importance"]])
imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = T),] ##extract importance

numLay = 6 ##number of layer to use in cLHS
lay <- rownames(imp)[1:numLay]
numTP <- 200 ###number of training points in each unit

library(clhs)

###loop through each unit and select cLHS points using layers based on RF importance
train2 <- train1[,c(lay,"layer")]
tp <- foreach(ss = unique(train2$layer), .combine = rbind) %do% {
  tSub <- train2[train2$layer == ss,]
  if(nrow(tSub) <= numTP){
    tSub
  }else{
    t1 <- clhs(tSub[,-length(tSub)], size = numTP, iter = 3000, simple = F, progress = T) ### create first slice
    tSub[t1$index_samples,]
  }
}

tp ###training points






####OLD CODE##################################################################
##############################################################################
t2 <- st_read(dsn = "TransectData/SBSmc2_1.1_1_NA", layer = "SBSmc2new (1)")
t2 <- t2["name"]
t2 <- st_transform(t2, 3005)
t2 <- st_zm(t2)

plot(t1[t1$name %in% c("POC","TP1","TP2","POT"),"name"], col = "black", add = T)
plot(t2, add = T)

s1 <- st_snap_points(t1[t1$name %in% c("POC","TP1"),"name"],t2)

s1 <- st_cast(s1, "POINT")
buf <- st_buffer(s1,0.5)
parts = st_collection_extract(lwgeom::st_split(t2$geometry, buf),"LINESTRING")
parts_all = st_as_sf(
  data.frame(
    id = 1:length(parts),
    geometry = parts
  )
)

st_snap_points = function(x, y, max_dist = 1000) {
  
  if (inherits(x, "sf")) n = nrow(x)
  if (inherits(x, "sfc")) n = length(x)
  
  out = do.call(c,
                lapply(seq(n), function(i) {
                  nrst = st_nearest_points(st_geometry(x)[i], y)
                  nrst_len = st_length(nrst)
                  nrst_mn = which.min(nrst_len)
                  if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                  return(st_cast(nrst[nrst_mn], "POINT")[2])
                })
  )
  return(out)
}
