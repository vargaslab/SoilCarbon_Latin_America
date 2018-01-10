library(raster)
library(reshape)
library(caret)
library(Metrics)
library(rgdal)
library(car)
library(automap)
library(raster)

#### Extract Data from WOSIS ##################################
### adapted from Eloi Ribero
setwd("WoSIS_2016_July/")
attributes = read.table("wosis_201607_attributes.txt", sep="\t",quote = "", header=TRUE)
profiles = read.table("wosis_201607_profiles.txt", sep="\t",quote = "", header=TRUE)
layers = read.table("wosis_201607_layers.txt", sep="\t", quote="", header=TRUE)

dim(attributes)
dim(profiles)
dim(layers)

colnames(attributes)
colnames(profiles)
colnames(layers)

# # display data
# attributes[1:22, 1:5]
# profiles[1:10, 1:6]
# layers[1:10, 1:6]

# merge profiles with layers
mer <- merge(x = profiles, y = layers, by = "profile_id", all = TRUE)

sel <- mer[mer$country_name %in% c("Argentina", "Belize", "Bolivia (Plurinational State of)",
                                   "Brazil", "Chile", "Costa Rica", "Colombia", "Cuba",
                                   "Dominican Republic", "Ecuador", "Guatemala", "Honduras",
                                   "Jamaica", "Mexico", "Nicaragua", "Panama", "Peru", "Suriname",
                                   "Uruguay", "Venezuela (Bolivarian Republic of)")
           
           & !is.na(mer[,"orgc_value_avg"]),]

sel <- sel[order(sel$profile_id, sel$top, sel$bottom),]
dim(sel)

dat=data.frame(id=as.character(sel$profile_id),
               country = sel$country_name,
               X=sel$longitude,
               Y=sel$latitude,
               top=sel$top,
               bottom=sel$bottom,
               SOC=sel$orgc_value_avg,
               BLD=sel$bdfi_value_avg,
               CRF=sel$cfvo_value_avg)

dat$country <- as.factor(as.character(dat$country))

#### Fill NAs, run mpslines ####################################################


load("DSM_supportfunctions.RData")

## Fixing missing CRF information, assuming no CRF
dat$CRF[is.na(dat$CRF)] <- 0

# BLD estimation from SOC Cookbook
# Available methods:
# Saini_1996, Drew_1973, Jeffrey_1979, Grigal_1989, Adams_1973, 
# Honeyset_Ratkowsky_1989
dat$BLD[is.na(dat$BLD)] <- estimateBD(dat$SOC[is.na(dat$BLD)], 
                                      method="Drew_1973")

summary(dat)
dat <- dat[complete.cases(dat),]

library(aqp)
depths(dat) <- id ~ top + bottom
site(dat) <- ~ X + Y + country
coordinates(dat) <- ~ X + Y

library(GSIF)

## convert 1 horizon profiles into  horizons

## Estimate 0-30 standard horizon usin mass preserving splines
try(SOC <- mpspline(dat, 'SOC', d = t(c(0,30))))
try(BLD <- mpspline(dat, 'BLD', d = t(c(0,30))))
try(CRFVOL <- mpspline(dat, 'CRF', d = t(c(0,30))))

## Prepare final data frame
dat <- data.frame(id = dat@site$id,
                  country = dat@site$country,
                  Y = dat@sp@coords[,2],
                  X = dat@sp@coords[,1],
                  SOC = SOC$var.std[,1],
                  BLD = BLD$var.std[,1],
                  CRFVOL = CRFVOL$var.std[,1])

dat <- dat[complete.cases(dat),]

# Estimate Organic Carbon Stock
# SOC must be in g/kg
# BLD in kg/m3
# CRF in percentage
# BLD error reported by Cuba: 20%
OCSKGM <- OCSKGM(ORCDRC = dat$SOC, BLD = dat$BLD*1000, CRFVOL = dat$CRFVOL, 
                 HSIZE = 30)

dat$OCSKGM <- OCSKGM
dat$meaERROR <- attr(OCSKGM,"measurementError")
dat <- dat[dat$OCSKGM>0,]



table(dat$country)

#write.csv(dat, "LAC_OCSKGM_WOSIS.csv")


res <- data.frame(country=character(), 
		n=numeric(),
		corSVM=numeric(),
		corRF=numeric(),
		corPL=numeric(), 
		corKK=numeric(),
		corRK=numeric(),
		errSVM=numeric(), 
		errRF=numeric(), 
		errPL=numeric(),	
    		errKK=numeric(), 
    		errRK=numeric(),
		idx=character(),
		n1=numeric())	
		res$country <- as.character(res$country)
		res$idx <- as.character(res$idx)

bestCor <- data.frame(country = character(), predictor = character(), 
			correlation = numeric())
bestCor$country <- as.character(bestCor$country)
bestCor$predictor <- as.character(bestCor$predictor)

performance <- data.frame(obs=numeric(), pred=numeric(), method=character(), country=character())
performance$method <- as.character(performance$method)
performance$country <- as.character(performance$country)


 #dat <- read.csv('LAC_OCSKGM_WOSIS.csv')

 targets <- dat[,5:8] 

 lev <- levels(dat$country)
 lev <- lev[-9]
	
countryID <- c('ARG', 'BLZ', 'BOL', 'BRA', 'CHL',
'COL', 'CRI', 'CUB', 'DOM', 'ECU', 'GTM', 'HND', 
'JAM', 'MEX' , 'NIC', 'PAN', 'PER', 'SUR', 'URY', 'VEN')

countryID <- countryID[-9]

 d <- data.frame( y=targets[,3], #######HERE
			country=dat$country, 
				dat[,3:4] )
				

for (i in 1:length(lev)){
#i=2
		train <- d[d$country==lev[i],]

		#train$y[train$y==0]<-NA
		
		n <- dim(train)[1]

		country <- lev[i]

		print(country)
		
		res[i, 1] <- country
		
		res[i, 2] <- n
		
		
grd.lst <- list.files(path = "~/Downloads/regmatspt1/worldGrids5km/", pattern="tif$")
COV1 <- stack(paste0("~/Downloads/regmatspt1/worldGrids5km/", grd.lst[grep(paste0(countryID[i], '_worldgridsCOVS_5km.tif'), grd.lst)]))
lev1 <- readRDS("~/Downloads/regmatspt1/worldGrids5km/worldgridsCOVS_names.rds")
names(COV1) <- lev1

grd.lst <- list.files(path = "~/Downloads/regmatspt1/topo5km/", pattern="tif$")
COV2 <- stack(paste0("~/Downloads/regmatspt1/topo5km/", grd.lst[grep(paste0(countryID[i], 'topo5km.tif'), grd.lst)]))
lev2 <- readRDS("~/Downloads/regmatspt1/topo5km/namesTOPO.rds")
names(COV2) <- lev2

COV <- stack(COV1, COV2)
rm(COV1);rm(COV2)

t <- train
coordinates(t) <- ~X+Y
train <- cbind (train[1], extract(COV, t), X=train$X, Y=train$Y)

		train$ln2dms3a  <- NULL	
		train$lnmdms3a  <- NULL
		
		train$DEM <- NULL
		train$lmtgsh3a <- NULL
		train$lmbgsh3a <- NULL
		train$glwwwf3a <- NULL
		train$wmkmod3a <- NULL
		train$cntgad3a <- NULL
		train$smkisr3a <- NULL
	
		cat1 <- grep('igb', names(train))[1:6]
		cat2 <- grep('esa', names(train))[23]		
		cat <- c(cat1, cat2)

		t <- na.omit(train[-cat])
		#if ( n>=10){
		
		COR <- cor(as.matrix(t[,1]), as.matrix(t[-c(1, 119, 120)]))

		x <- subset(melt(COR), value != 1 | value != NA)
		x <- x[with(x, order(-abs(x$value))),]
		
		names(x)[1] <- 'country'
		names(x)[2] <- 'predictor'
		names(x)[3] <- 'correlation'
		x$country <- country

		bestCor <- rbind (bestCor, x[1:10,])			
	
		idx <- as.character(x$predictor[1:10])
		
		t <- train[c('y')]

		train_sp <- cbind(t, train[c(idx, 'X', 'Y')])
		train_sp <- na.omit(train_sp)
		coordinates(train_sp)=~X+Y

		train <- cbind(t, train[idx])

		n1 <- dim(na.omit(train))[1]
		train <- na.omit(train)
		res[i, 14] <- n1
		 
(fm = as.formula(paste("y ~", paste0(names(train[-1]), collapse = "+"))))	

ctrl <- trainControl(method = "cv", savePred=T)

	  mod1 <- train(fm, data=train, method = "svmLinear", trControl = ctrl)
	  resid1=mod1$pred
	  mod2 <- train(fm, data=train, method = "rf", trControl = ctrl)
	  resid2=mod2$pred
	  mod3 <- train(fm, data=train, method = "pls", trControl = ctrl)
	  resid3=mod3$pred
	  mod4 <- train(fm, data=train, method = "kknn", trControl = ctrl)
	  resid4=mod4$pred

 modelo.MLR <- lm(y ~ ., data = train_sp@data) 
 modelo.MLR.step <- step(modelo.MLR, direction="both")
 
 proj4string(train_sp) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
 train_sp <- spTransform(train_sp, CRS("+init=epsg:3857"))
 train_sp  <- train_sp [which(!duplicated(train_sp@coords)), ]
   

   n <- dim(train_sp)[1]

 err1 <- rmse(resid1$obs, resid1$pred)
	  err2 <- rmse(resid2$obs, resid2$pred)
	  err3 <- rmse(resid3$obs, resid3$pred)
	  err4 <- rmse(resid4$obs, as.numeric(resid4$pred))
	  cor1 <- cor(resid1$obs, resid1$pred)
	  cor2 <- cor(resid2$obs, resid2$pred)
	  cor3 <- cor(resid3$obs, resid3$pred)
	  cor4 <- cor(resid4$obs, as.numeric(resid4$pred))
	  
	  res[i, 3] <- round(cor1, 3)
	  res[i, 4] <- round(cor2, 3)
	  res[i, 5] <- round(cor3, 3)
	  res[i, 6] <- round(cor4, 3)
	  res[i, 8] <- round(err1, 3)
	  res[i, 9] <- round(err2, 3)
	  res[i, 10] <- round(err3, 3)
          res[i, 11] <- round(err4, 3)
	  IDX <- paste0(idx[1],' + ', idx[2],' + ', idx[3],' + ', idx[4], ' + ',idx[5], ' + ',idx[6], ' + ',idx[7], ' + ',idx[8], ' + ',idx[9], ' + ',idx[10])
	  res[i, 13] <- IDX

SVM <- data.frame(obs=resid1$obs,pred=resid1$pred,  method='SVM', country=country)
RF  <- data.frame(obs=resid2$obs,pred=resid2$pred,  method='RF', country=country)
PL  <- data.frame(obs=resid3$obs,pred=resid3$pred,  method='PL', country=country)
KK  <- data.frame(obs=resid4$obs,pred=resid4$pred,  method='KK', country=country)

performance <- rbind(performance, SVM, RF, PL, KK)


beginCluster(3)
pred1 <- clusterR(COV[[idx]], predict, args=list(model=mod1))
pred2 <- clusterR(COV[[idx]], predict, args=list(model=mod2))
pred3 <- clusterR(COV[[idx]], predict, args=list(model=mod3))
pred4 <- clusterR(COV[[idx]], predict, args=list(model=mod4))
s <- stack(pred1, pred2, pred3, pred4)
endCluster()

  
	if ( n>=10){

 OCS.krige.cv <- autoKrige.cv(formula =fm, 
                       input_data = train_sp,  nfold=10)
  	 
	
val <-data.frame(OCS.krige.cv[1][[1]])
	  err5 <- rmse(val$observed, val$var1.pred)
	  cor5 <- cor(val$observed, val$var1.pred)
	  res[i, 7] <- round(cor5, 3)
	  res[i, 12] <- round(err5, 4)
	  
RK <- data.frame(obs=val$observed, pred=val$var1.pred, method='RK', country=country)

performance <- rbind(performance, RK)

COV <- projectRaster(COV[[idx]], crs = CRS(projection(train_sp)), method='ngb')
COV.sp <- as(COV, "SpatialGridDataFrame")

OCS.krige <- autoKrige(formula = as.formula(modelo.MLR.step$call$formula), 
                       input_data = train_sp, 
                       new_data = COV.sp,
                       verbose = TRUE,
                       block = c(1000, 1000))
RKprediction <- expm1(raster(OCS.krige$krige_output[1]))
RKpredsd <- expm1(raster(OCS.krige$krige_output[3]))

s2 <- stack(RKprediction, RKpredsd)

s <- stack (s, projectRaster(s2, s))


}


writeRaster(s, paste0('/home/mario/Downloads/regmatspt1/topo5km/predicted/',country,'_SOC_predictions_10predictors_oct19_crf.tif'), overwrite=TRUE)


}

per <- readRDS('/home/mario/Downloads/data/WoSIS_2016_July/OCS.rds')
performance <-per[[2]]

#performance$BLD <- 'BLD'

performance$model <- performance$method
performance$site <- performance$country
performance$pred <- as.numeric(performance$pred)
performance <- na.omit(performance)
performance$pred[performance$pred<=0] <- 0
performance$obs[performance$obs<=0] <- 0
performance$pred[performance$pred>=max(performance$obs)] <- NA
performance <- na.omit(performance)
#TaylorDiagram(performance, obs = "obs", mod = "pred", 
#                group = "model", normalise=FALSE, type='site', arrow.lwd=5)

TaylorDiagram(performance, obs = "obs", mod = "pred", 
                group = "model", normalise=FALSE,  arrow.lwd=5,  main='CRF (%)')




lev <- levels (performance$country)
par(mfrow=c(5, 4))
for (i in 1 : length (lev)){



n <- paste0(lev[i], "_taylorDiag.pdf")
dev.copy(png,n)

i=4
x <- performance[performance$country==lev[i],]

TaylorDiagram(x, obs = "obs", mod = "pred", 
                group = "model", normalise=TRUE, main=lev[i])


#http://www.openair-project.org/PDF/OpenAir_Manual.pdf
library(openair)
options(digits = 2)
(statSite <- modStats(performance, obs = "obs", mod = "pred", type = "site"))
(statModel <- modStats(performance, obs = "obs", mod = "pred", type = "model"))
l <- list (statSite, statModel)
saveRDS(l, 'StatsByModelandCountryOCS.rds')








