######################
## Random Forest script that includes:
## Extraction of covariates to points
## Prediction interval creation
## Cross Validation
## Most steps parallelized
######################
## Predictions of air dry water content (w15l2)


# Workspace setup
# Install packages if not already installed
required.packages <- c("raster", "sp", "rgdal", "randomForest", "snow", "snowfall", "quantregForest", "doSNOW","dplyr","itertools","foreach")# might need snowfall
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
#memory.limit(500000) # windows only
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)
#options(scipen = 999) #turns off scientific notation
#par(mar=c(0.3,0.3,0.3,0.3)) changes plot margins

######## Load shapefile ##############
#setwd("C:/Models_active_work/UpCo/ECmodel_wLIMS")## FOlder with points
#shp.pts <-readOGR(".", "ec_12pre_ncss_LIMS_UPCO")
#point.proj <- projection(shp.pts)
## If not prj file and you know proj, can specify by name
#shp.proj <- CRS("+proj=longlat +datum=WGS84")

######## Get points for extraction if in table form ###########
setwd("/home/tnaum/data/BLMsoils/pH_h2o_3D_pedonCV")
pts <- read.delim("NCSS17_pH_carbonate_gyp_ttab.txt") # If in delimited file other than csv
### Weed out points with imprecise coordinates ###
pts$latnchar = nchar(abs(pts$latitude_decimal_degrees))
pts$longnchar = nchar(abs(pts$longitude_decimal_degrees))
ptsc = subset(pts, pts$latnchar > 5 & pts$longnchar > 6)
### Turn into spatial file
shp.pts <- ptsc
coordinates(shp.pts) <- ~ longitude_decimal_degrees + latitude_decimal_degrees
temp.proj <- CRS("+proj=longlat +datum=WGS84") ## specify projection
projection(shp.pts) <- temp.proj


######## Load map clip boundary (if needed) ###########
setwd("/home/tnaum/Dropbox/USGS/BLM_projects/Utah_BLM_Salinity/Huc6_boundary")
polybound <- readOGR(".", "CO_River_watershed_Meade_alb")
polybound <- spTransform(polybound, temp.proj)
## Now clip points and check with visualization
shp.pts = shp.pts[polybound,]#clip by outer extent of all polybound features
plot(polybound)
plot(shp.pts, add=TRUE)
shp.pts$depth = (shp.pts$hzn_bot+shp.pts$hzn_top)/2

######### Grid Prep #################
## Make list of grids
setwd("/home/tnaum/data/UCRB_Covariates")
cov.grids <- list.files(pattern=".tif$")
## If points need to be matched up to grids ###
projgrid = raster("ELEVm.tif")
## Or Make a stack of grids to extract all at once (for smaller datasets)
#cov.stack <- stack()
cov.proj <- projection(projgrid)
shp.pts <- spTransform(shp.pts, CRS(cov.proj)) # project to match rasters

## Plot to ensure alignment bw points and rasters
plot(projgrid)
plot(shp.pts, add=TRUE)

## Parallelized extract: (larger datasets)
cpus = 30
sfInit(parallel=TRUE, cpus=cpus)
sfExport("shp.pts", "cov.grids")
sfLibrary(raster)
sfLibrary(rgdal)
ov.lst <- sfLapply(cov.grids, function(i){try( raster::extract(raster(i), shp.pts) )}) 
snowfall::sfStop()
detach(package:snowfall, unload=TRUE)
ov.lst <- as.data.frame(ov.lst)
names(ov.lst) = tools::file_path_sans_ext(basename(cov.grids))
ov.lst$DID <- seq.int(nrow(ov.lst))
shp.pts$DID = seq.int(nrow(shp.pts))
pts.ext = merge(as.data.frame(shp.pts),ov.lst, by="DID")

## Save points
setwd("/home/tnaum/data/BLMsoils/pH_h2o_3D_pedonCV")
write.table(pts.ext, "cop_ncss17pH_h20_covarsc.txt", sep = "\t", row.names = FALSE)
## Pre-prepared points with variables extracted.
pts.ext <- read.delim("cop_ncss17pH_h20_covarsc.txt", stringsAsFactors = F)## File available on repository
pts.ext$depth <- (pts.ext$hzn_bot+pts.ext$hzn_top)/2


## Prep for Random Forest
ptspred.list <- gsub(".tif","", cov.grids)# Take .tif off of the grid list to just get the variable names
prop <- "ph_h2o" ## Dependent variable
depth <- "depth" # to use in 3D soil property mapping
ptspred.list <- c(ptspred.list,"prop",depth) #Add dependent variable
pts.ext$prop <- pts.ext$ph_h2o ## UPDATE EVERY TIME!!!
pts.extc <- pts.ext[c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use 
pts.extc$depth <- as.numeric(as.character(pts.extc$depth))
pts.extc <- na.omit(pts.extc)# Remove any record with NA's (in any column - be careful)
xtrain <- as.matrix(pts.extc[c("depth",gsub(".tif","", cov.grids))])
ytrain <- c(as.matrix(pts.extc$prop))
## If transformation is needed
# logytrain <- log(ytrain)
# sqrtytrain <- sqrt(ytrain)

############### Build quantile Random Forest
Qsoiclass <- quantregForest(x=xtrain, y=ytrain, importance=TRUE, ntree=100, keep.forest=TRUE)
#soiclass = randomForest(ec_12pre ~ ., data = ptsc, importance=TRUE, proximity=FALSE, ntree=100, keep.forest=TRUE)
soiclass <- Qsoiclass
class(soiclass) <- "randomForest"
soiclass## Get oob error
varImpPlot(soiclass)
saveRDS(Qsoiclass,"Qsoiclass_RFmodel_pH_h2o.rds")

##### Raster Prediction Step 
depths <- c(0) # vector of depths to be predicted
for (d in depths){
  #depth <- calc(projgrid, fun=function(x)(ifelse(x>-999,d)), progress="text") #raster to set prediction depth
  beginCluster(30,type='SOCK')
  depth.fn <- function(x) {
    ind <- ifelse(x>-999,d)
    return(ind) 
  }
  depth <- clusterR(projgrid, calc, args=list(fun=depth.fn),progress='text',export='d')
  gc() # Flush out RAM
  #depth <- projgrid
  #depth <- setValues(depth, d) #Second argument is the depth value
  names(depth)<-"depth"
  # Reference covar rasters to use in prediction
  setwd("/home/tnaum/data/UCRB_Covariates")
  rasters <- stack(cov.grids, depth)
  ## Parallelized predict
  rasterOptions(maxmemory = 3e+09, chunksize = 1e+08)
  predl <- clusterR(rasters, predict, args=list(model=Qsoiclass,what=c(0.025)),progress="text")
  predh <- clusterR(rasters, predict, args=list(model=Qsoiclass,what=c(0.975)),progress="text")
  pred <- clusterR(rasters, predict, args=list(model=soiclass),progress="text")
  ## Untransformted PI width prep steps
  s <- stack(predh,predl)
  PIwidth.fn <- function(a,b) {
    ind <- a-b
    return(ind)
  }
  PIwidth <- clusterR(s, overlay, args=list(fun=PIwidth.fn),progress = "text")
  # Determine 95% interquantile range of original training data for horizons that include the depth being predicted
  pts.ext.pred.depth <- subset(pts.ext, as.numeric(pts.ext$hzn_top) <= d & as.numeric(pts.ext$hzn_bot) > d)
  varrange <- as.numeric(quantile(pts.ext.pred.depth$prop, probs=c(0.975), na.rm=T)-quantile(pts.ext.pred.depth$prop, probs=c(0.025),na.rm=T)) ## TRANSFORM IF NEEDED!
  PIrelwidth.fn <- function(a,b) {
    ind <- (a-b)/varrange
    return(ind)
  }
  PIrelwidth <- clusterR(s, overlay, args=list(fun=PIrelwidth.fn),progress = "text",export='varrange')
  bt.fn <- function(x) {
    ind <- (x)^2 #If a backtransform is needed 10^(x) or exp(x) or ^2
    return(ind) 
  }
  #predh_bt = clusterR(predh, calc, args=list(fun=bt.fn),progress='text')
  #predl_bt = clusterR(predl, calc, args=list(fun=bt.fn),progress='text')
  #pred_bt = clusterR(pred, calc, args=list(fun=bt.fn),progress='text')
  #s_bt = stack(predh_bt,predl_bt)
  PIwidth_bt.fn <- function(a,b) {
    ind <- a-b
    return(ind)
  }
  #PIwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIwidth_bt.fn),progress = "text")
  ## If transformed, use the following code for PI width prep steps
  PIrelwidth_bt.fn <- function(a,b) {
    ind <- (a-b)/varrange
    return(ind)
  }
  #PIrelwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIrelwidth_bt.fn),progress = "text", export='varrange')
  endCluster()
  ## Write new geotiff files
  setwd("/home/tnaum/data/BLMsoils/pH_h2o_3D_pedonCV/Outputs")
  ## Untranformed code block
  writeRaster(pred, overwrite=TRUE,filename=paste(prop,d,"cm_QRF.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  writeRaster(predl, overwrite=TRUE,filename=paste(prop,d,"cm_QRF_95PI_l.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  writeRaster(predh, overwrite=TRUE,filename=paste(prop,d,"cm_QRF_95PI_h.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  writeRaster(PIrelwidth, overwrite=TRUE,filename=paste(prop,d,"cm_QRF_95PI_relwidth.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  writeRaster(PIwidth, overwrite=TRUE,filename=paste(prop,d,"cm_QRF_95PI_width.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  ## Transformed code block
  #writeRaster(pred_bt, overwrite=TRUE,filename=paste(prop,d,"cm_QRF_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  #writeRaster(predl_bt, overwrite=TRUE,filename=paste(prop,d,"cm_QRF_95PI_l_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  #writeRaster(predh_bt, overwrite=TRUE,filename=paste(prop,d,"cm_QRF_95PI_h_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  #writeRaster(PIwidth_bt, overwrite=TRUE,filename=paste(prop,d,"cm_QRF_95PI_width_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  #writeRaster(PIrelwidth_bt, overwrite=TRUE,filename=paste(prop,d,"cm_QRF_95PI_relwidth_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  print(paste("Done with ", d))
}


################### Pedon-based Cross validation ################################
pedons <- unique(pts.ext$pedon_key)
pedondf <- as.data.frame(pedons)
nfolds <- 10
pedondf$folds <- sample.int(nfolds,size =length(pedons),replace=T)
ptspred.listcv <- c(ptspred.list,"pedon_key")
pts.extcv <- pts.ext[c(ptspred.listcv)]
pts.extcv$depth <- as.numeric(as.character(pts.extcv$depth))
pts.extcv <- na.omit(pts.extcv)# Remove any record with NA's (in any column - be careful)
pts.extcv <- cbind(pedondf[match(pts.extcv$pedon_key, pedondf$pedons),], pts.extcv) #attaches fold sample to horizon records
pts.extcv$prop_t <- pts.extcv$prop #transformed: if needed
formulaStringCVp <- as.formula(paste('prop_t ~','depth','+', paste(gsub(".tif","", cov.grids), collapse="+")))
pts.extcv$pcvpred <- "NA"
#for (g in seq(nfolds)){
CV_factorRF <- function(g,pts.extcv, formulaStringCVp){
  traindf <- subset(pts.extcv, pts.extcv$folds != g)
  testdf <- subset(pts.extcv, pts.extcv$folds == g)
  rf.pcv <- randomForest(formulaStringCVp, data=traindf, importance=FALSE, proximity=FALSE, ntree=100, keep.forest=TRUE)
  testdf$pcvpred <- as.character(predict(rf.pcv, newdata=testdf))
  return(testdf)
  #pts.extcv$pcvpred <- ifelse(pts.extcv$folds == g, predict(rf.pcv, newdata=pts.extcv),pts.extcv$pcvpred)
  #print(g)
}
snowfall::sfInit(parallel=TRUE, cpus=nfolds)
snowfall::sfExport("pts.extcv","formulaStringCVp","CV_factorRF")
snowfall::sfLibrary(randomForest)
pts.extpcv <- snowfall::sfLapply(1:nfolds, function(g){CV_factorRF(g, pts.extcv=pts.extcv,formulaStringCVp=formulaStringCVp)})
snowfall::sfStop()
pts.extpcv <- plyr::rbind.fill(pts.extpcv)
pts.extpcv$pcvpred = as.numeric(pts.extpcv$pcvpred)
cvp.RMSE = sqrt(mean((pts.extpcv$prop_t - pts.extpcv$pcvpred)^2, na.rm=TRUE))
cvp.Rsquared = 1-var(pts.extpcv$prop_t - pts.extpcv$pcvpred, na.rm=TRUE)/var(pts.extpcv$prop_t, na.rm=TRUE)
# Back transformed: create pcvpred_bt even if not tranformed for cv.depth function
pts.extpcv$pcvpred_bt <- pts.extpcv$pcvpred
cvp.RMSE_bt = sqrt(mean((pts.extpcv$prop - pts.extpcv$pcvpred_bt)^2, na.rm=TRUE))
cvp.Rsquared_bt = 1-var(pts.extpcv$prop - pts.extpcv$pcvpred_bt, na.rm=TRUE)/var(pts.extpcv$prop, na.rm=TRUE)
## Now separate out a CV for one depth to graph
setwd("C:/models/2Dv3D/pH_h2o_3D_pedonCV")
pts.ext.0cm <- subset(pts.ext, as.numeric(pts.ext$hzn_top) <= 0 & as.numeric(pts.ext$hzn_bot) > 0)
pts.ext.0cm$depth <- (pts.ext.0cm$hzn_bot + pts.ext.0cm$hzn_top) / 2
pts.ext.0cm$cvid <- paste(pts.ext.0cm$pedon_key,pts.ext.0cm$depth,pts.ext.0cm$prop, sep="") # specify the predicted variable...
pts.extpcv.cv <- pts.extpcv
pts.extpcv.cv$cvid <- paste(pts.extpcv.cv$pedons,pts.extpcv.cv$depth,pts.extpcv.cv$prop, sep="")
pts.ext.0cm.cv <- subset(pts.ext.0cm, select=c("cvid"))
pts.ext.0cm.cv <- merge(pts.ext.0cm,pts.extpcv.cv, by="cvid")
saveRDS(pts.ext.0cm.cv, "UCRB_ph_h2o_PCV_pts_0cm.rds")
rm(pts.ext.0cm,pts.extpcv.cv,pts.ext.0cm.cv)


## Parse Cross Validation by depth
# need to re-match with original horizons
depths <- c(0,15,30,60,100) # Depths to be evaluated for CIs
cv.depth <- function(d, pts.ext, pts.extpcv){ ## fucntion definition: d is variable iterator name for depths for use in lapply
#for(d in depths){
  pts.ext.depth <- subset(pts.ext, as.numeric(pts.ext$hzn_top) <= d & as.numeric(pts.ext$hzn_bot) > d)
  cvdepth <- d
  pts.ext.depth$depth <- (pts.ext.depth$hzn_bot + pts.ext.depth$hzn_top) / 2
  pts.ext.depth$cvid <- paste(pts.ext.depth$pedon_key,pts.ext.depth$depth,pts.ext.depth$prop, sep="") # specify the predicted variable...
  pts.extpcv.cv <- pts.extpcv
  pts.extpcv.cv$cvid <- paste(pts.extpcv.cv$pedons,pts.extpcv.cv$depth,pts.extpcv.cv$prop, sep="")
  pts.ext.depth.cv <- subset(pts.ext.depth, select=c("cvid"))
  pts.ext.depth.cv <- merge(pts.ext.depth.cv,pts.extpcv.cv, by="cvid")
  RMSE <- sqrt(mean((pts.ext.depth.cv$prop - pts.ext.depth.cv$pcvpred_bt)^2, na.rm=TRUE)) # must specify predicted variable
  Rsq <- 1-var(pts.ext.depth.cv$prop - pts.ext.depth.cv$pcvpred_bt, na.rm=TRUE)/var(pts.ext.depth.cv$prop, na.rm=TRUE)# must specify predicted variable
  n <- length(pts.ext.depth.cv$cvid)
  cvdf <- data.frame(cvdepth,RMSE,Rsq,n)
  names(cvdf) <- c("cvdepth","RMSE","Rsq","n")
  return(cvdf)
  gc()
}
cv.list <- lapply(depths, function(d){try(cv.depth(d,pts.ext,pts.extpcv))})
cv.df <- cv.list[[1]] 
cv.df <- cv.df[FALSE,]
for(i in seq(1:length(cv.list))){
  newrow <- cv.list[[i]]
  if(class(newrow)=="data.frame"){
    cv.df <- rbind(cv.df, newrow)
  }
  print(paste("Done with ", i, sep=""))
}
# Save table to folder
setwd("/home/tnaum/data/BLMsoils/pH_h2o_3D_pedonCV")
write.table(cv.df, "UCRB_ph_h2o_CV_by_depth.txt", sep = "\t", row.names = FALSE)


##### Relative PI Interval statistics for different depths
pred.pts.ext <- read.delim("/home/tnaum/data/BLMsoils/UCRB_Summary_pts/UCRB_ncss17_sample_subset_covarsc.txt") ## prediction summary locations
pred.pts.extdf <- as.data.frame(pred.pts.ext)
depths <- c(0,15,30,60,100) # Depths to be evaluated for PIs
relPI.depth <- function(d, pts.ext, pred.pts.extdf,Qsoiclass){
  num_splits <- 30 # number of cpus to use
  cl <- makeCluster(num_splits)
  registerDoSNOW(cl)
  pred.pts.extdf$depth<-d
  DF_pred_l<-
    foreach(d=isplitRows(pred.pts.extdf, chunks=num_splits),.combine=c, .packages=c("quantregForest", "randomForest")) %dopar% {
      predict(Qsoiclass,newdata=d, what=c(0.025))
    }
  DF_pred_h<-
    foreach(d=isplitRows(pred.pts.extdf, chunks=num_splits),.combine=c, .packages=c("quantregForest", "randomForest")) %dopar% {
      predict(Qsoiclass,newdata=d, what=c(0.975))
    }
  stopCluster(cl)
  ## Now compute stats for PI intervals
  pts.ext.depth <- subset(pts.ext, as.numeric(pts.ext$hzn_top) <= d & as.numeric(pts.ext$hzn_bot) > d)
  ## Untransform
  DF_pred_l_bt <- DF_pred_l # Prepare for all runs, just don't use back transform fn if it wasn't used 
  DF_pred_h_bt <- DF_pred_h # Prepare for all runs, just don't use back transform fn if it wasn't used 
  df.PIwidth=DF_pred_h_bt-DF_pred_l_bt
  df.varrange=as.numeric(quantile(pts.ext.depth$prop, probs=c(0.975),na.rm=T)-quantile(pts.ext.depth$prop, probs=c(0.025),na.rm=T))##Change soil factor name
  df.PIrelwidth=(DF_pred_h_bt-DF_pred_l_bt)/df.varrange
  ## Stats (untransformed)
  PIdepth <- d
  mean <- mean(df.PIwidth)
  median <- median(df.PIwidth)
  q25 <- unname(quantile(df.PIwidth, probs=c(0.25),na.rm=T))
  q75 <- unname(quantile(df.PIwidth, probs=c(0.75),na.rm=T))
  sd <-sd(df.PIwidth)
  max <- max(df.PIwidth)
  min <- min(df.PIwidth)
  relmean <- mean(df.PIrelwidth)
  relmedian <- median(df.PIrelwidth)
  relq25 <- unname(quantile(df.PIrelwidth, probs=c(0.25),na.rm=T))
  relq75 <- unname(quantile(df.PIrelwidth, probs=c(0.75),na.rm=T))
  relsd <-sd(df.PIrelwidth)
  relmax <- max(df.PIrelwidth)
  relmin <- min(df.PIrelwidth)
  PIdf <- data.frame(PIdepth, mean, median, q25,q75,sd,max,min,relmean, relmedian, relq25,relq75,relsd,relmax,relmin)
  names(PIdf) <- c("PIdepth","mean","median", "q25","q75","sd","max","min","relmean","relmedian", "relq25","relq75","relsd","relmax","relmin")
  return(PIdf)
}
relPI.list <- lapply(depths, function(d){try(relPI.depth(d,pts.ext, pred.pts.extdf,Qsoiclass))})
relPI.df <- relPI.list[[1]] 
relPI.df <- relPI.df[FALSE,]
for(i in seq(1:length(relPI.list))){
  newrow <- relPI.list[[i]]
  if(class(newrow)=="data.frame"){
    relPI.df <- rbind(relPI.df, newrow)
  }
  print(paste("Done with ", i, sep=""))
}
# Save table to folder
setwd("/home/tnaum/data/BLMsoils/pH_h2o_3D_pedonCV")
write.table(relPI.df, "UCRB_ph_h2o_relPI_by_depth.txt", sep = "\t", row.names = FALSE)


