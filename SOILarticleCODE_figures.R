library(raster)
library(mapdata)
library(maps)
library(rgeos)
library(maptools)
library(dichromat)

setwd("~/Downloads/regmatspt1/topo5km/predicted")
e<- readRDS('extent.rds')
mer <- raster('filledSOCmap.tif')

lis <- list.files()[grep('SOC_predictions_10predictors_oct19.tif', list.files())]
r <- raster('spatialReference5x5kmgrid.tif')

e<- readRDS('extent.rds')
lis

li <- list()

for(i in 1:length(lis)){
s <- stack(lis[i])
li[[i]] <- resample(s, r, method='ngb')
print(i)
}

per <- readRDS('/home/mario/Downloads/data/WoSIS_2016_July/OCS.rds')
performance <-per[[3]]

norm <- function (x) (x-min(x))/(max(x)-min(x))
performance <- data.frame(apply(na.omit(performance[,3:12]), 2, norm))
performance$errCorRatioSVM <- performance$errSVM/performance$corSVM
performance$errCorRatioRF <- performance$errRF/performance$corRF
performance$errCorRatioPL<- performance$errPL/performance$corPL
performance$errCorRatioKK<- performance$errKK/performance$corKK
performance$errCorRatioRK <- performance$errRK/performance$corRK
w <- apply(performance[11:15], 2, median)

svmList <- list()
rfList <- list()
plList <- list()
kkList <- list()
rkList <- list()

for(i in 1:length(lis)){

SVM <- stack(lis[i], bands=1)
SVM <- resample(SVM, r)
RF <- stack(lis[i], bands=2)
RF <- resample(RF, r)
PL <- stack(lis[i], bands=3)
PL <- resample(PL, r)
KK <- stack(lis[i], bands=4)
KK <- resample(KK, r)
RK <- stack(lis[i], bands=5)
RK <- resample(RK, r)

svmList[[i]] <- SVM
rfList[[i]] <- RF
plList[[i]] <- PL
kkList[[i]] <- KK
rkList[[i]] <- RK

print(i)

}

l <- svmList
SVMmap <- mosaic(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]], l[[6]], l[[7]], l[[8]], l[[9]], l[[10]], l[[11]],l[[12]], l[[11]], l[[12]], l[[13]], l[[14]], l[[15]], l[[16]], l[[17]], l[[18]], l[[19]], fun=median)
l <- rfList
RFmap <- mosaic(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]], l[[6]], l[[7]], l[[8]], l[[9]], l[[10]], l[[11]],l[[12]], l[[11]], l[[12]], l[[13]], l[[14]], l[[15]], l[[16]], l[[17]], l[[18]], l[[19]], fun=median)
l <- plList
PLmap <- mosaic(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]], l[[6]], l[[7]], l[[8]], l[[9]], l[[10]], l[[11]],l[[12]], l[[11]], l[[12]], l[[13]], l[[14]], l[[15]], l[[16]], l[[17]], l[[18]], l[[19]], fun=median)
l <- kkList
KKmap <- mosaic(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]], l[[6]], l[[7]], l[[8]], l[[9]], l[[10]], l[[11]],l[[12]], l[[11]], l[[12]], l[[13]], l[[14]], l[[15]], l[[16]], l[[17]], l[[18]], l[[19]], fun=median)
l <- rkList
RKmap <- mosaic(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]], l[[6]], l[[7]], l[[8]], l[[9]], l[[10]], l[[11]],l[[12]], l[[11]], l[[12]], l[[13]], l[[14]], l[[15]], l[[16]], l[[17]], l[[18]], l[[19]], fun=median)

 me <- GSIF::merge(RKmap, SVMmap, KKmap, PLmap, RFmap, RMSE.l=w)

me[is.infinite(me)==TRUE]<- NA

x <- data.frame(a =numeric(), b=numeric())
for (i in 1:8){
myRaster <- (me-quantile(me, probs = c(0.025, 0.975), type=7,names = FALSE)[1])/(quantile(me, probs = c(0.025, 0.975), type=7,names = FALSE)[2]-quantile(me, probs = c(0.025, 0.975), type=7,names = FALSE)[1])
#myRaster[myRaster>=10]<-10
myRaster[is.infinite(myRaster)==TRUE]<- NA
myRaster[myRaster>=i]<-i
x[i, 1]<-(cellStats(exp(myRaster), sum)*2.5e+07)*1e-12
x[i, 2]<-i
print(i)
}

par(mar=c(5,5,5,5))
par(mfrow=c(2,1))


names(w) <- c('SVM', 'RF', 'PL', 'KK', 'RK')
 barplot(sort(w), las=2, ylab='rmse/correlation ratio', cex.axis=1.5, cex.lab=1.5, cex.names=1.3, xlab='')
plot(exp(x$b), x$a, ylab='SOC stock (Pg)', xlab='SOC max limit (kg/m/30cm)', cex.axis=1.5, cex.lab=1.5)
lines(x$b, x$a , col = "gray")

 cellStats(exp(mer), sum)
[1] 1246792
 cellStats(exp(mer), sum)*2.5e+07
[1] 3.116979e+13


#myRaster <- mer
ext <- as.vector(extent(myRaster))
boundaries <- map('worldHires', fill=TRUE,
    xlim=ext[1:2], ylim=ext[3:4],
    plot=FALSE)
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                               proj4string=CRS(projection(myRaster)))
	

plot(e, main='', col='gray')
plot(myRaster, add=TRUE, col=dichromat(rev(terrain.colors(15))))

plot(myRaster,zlim=c(1, 3), add=TRUE, legend=FALSE, col=dichromat(rev(terrain.colors(15))))


plot(myRaster,zlim=c(0, 0.5), add=TRUE, legend=FALSE, col=dichromat(rev(terrain.colors(15))))
plot(myRaster,zlim=c(0.5, 0.7), add=TRUE, legend=FALSE, col=dichromat(rev(terrain.colors(15))))
plot(myRaster,zlim=c(0.7, 1), add=TRUE, legend=FALSE, col=dichromat(rev(terrain.colors(15))))
plot(bPols, add=TRUE)

#writeRaster(myRaster, file='SOCpred1.tif', overwrite=TRUE)

newproj <- '+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '

r <- aggregate(r, 5, min)
r <- mask (r, bPols)
ref <- projectRaster(r, crs=newproj)
ref <- as(ref, 'SpatialPixelsDataFrame')
meuse.grid <- ref

samp <- myRaster
#samp[samp>=2] <- NA
x <- sampleRandom(samp, size=15000, na.rm=TRUE, sp=TRUE, asRaster=FALSE)

x <- spTransform(x, CRS=newproj)
meuse <- x
names(meuse)[1] <- 'SOCmosaic'

library(raster)
library(automap)
library('sp')
library('gstat')
library('parallel')
library(automap)

model <- autofitVariogram(SOCmosaic ~ 1, input_data = meuse, verbose=c(T,T), 
GLS.model=TRUE)
m <- model$var_model
plot(model)


no_cores <- detectCores() - 1

cl <- makeCluster(no_cores)

parts <- split(x = 1:length(meuse.grid), f = 1:no_cores)

clusterExport(cl = cl, varlist = c("meuse", "meuse.grid", "m", "parts"), envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))

system.time(parallelX <- parLapply(cl = cl, X = 1:no_cores, fun = function(x) krige(formula = SOCmosaic~1, locations = meuse, newdata = meuse.grid[parts[[x]],], model = m)))

stopCluster(cl)

mergeParallelX <- maptools::spRbind(parallelX[[1]], parallelX[[2]])
mergeParallelX <- maptools::spRbind(mergeParallelX, parallelX[[3]])
mergeParallelX <- maptools::spRbind(mergeParallelX, parallelX[[4]])
mergeParallelX <- maptools::spRbind(mergeParallelX, parallelX[[5]])
mergeParallelX <- maptools::spRbind(mergeParallelX, parallelX[[6]])
mergeParallelX <- maptools::spRbind(mergeParallelX, parallelX[[7]])

# Create SpatialPixelsDataFrame from mergeParallelX
mergeParallelX <- SpatialPixelsDataFrame(points = mergeParallelX, data = mergeParallelX@data)

spplot(raster(mergeParallelX["var1.pred"]), main = "ordinary kriging predictions")

pro <- projectRaster(raster(mergeParallelX["var1.pred"]),myRaster)

mer <- merge(myRaster, pro, RMSE.l=c(1, 2))

mer <- mask (mer, bPols[-55])


library(raster)
library(mapdata)
library(maps)
library(rgeos)
library(maptools)
library(dichromat)

setwd("~/Downloads/regmatspt1/topo5km/predicted")
e<- readRDS('extent.rds')
mer <- raster('filledSOCmap.tif')

lis <- list.files()[grep('SOC_predictions_10predictors_oct19.tif', list.files())]
r <- raster('spatialReference5x5kmgrid.tif')

myRaster <- mer
ext <- as.vector(extent(myRaster))
boundaries <- map('worldHires', fill=TRUE,
    xlim=ext[1:2], ylim=ext[3:4],
    plot=FALSE)
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                               proj4string=CRS(projection(myRaster)))

mer <-mask (mer, bPols[-55])


plot(e, main='', col='gray100')
plot(mer, add=TRUE, legend=FALSE, col=dichromat(rev(terrain.colors(15))))
plot(mer,zlim=c(0, 0.5), add=TRUE, legend=FALSE, col=dichromat(rev(terrain.colors(15))))
plot(mer,zlim=c(0.5, 1), add=TRUE, legend=TRUE, col=dichromat(rev(terrain.colors(15))))
plot(bPols[-55], add=TRUE, cex=2)

l <- list()
for(i in 1:length(lis)){
if (dim(stack(lis[i]))[3]<6){

s <- stack(lis[i])[[1:4]]
u <- calc(s, sd) / calc(s, median)
u <- resample(u, r)
l[[i]] <- u
print(i)

}

else{
s <- stack(lis[i])[[1:5]]
u <- exp(calc(s, sd)) / exp(calc(s, median))
u <- resample(u, r)
l[[i]] <- u
print(i)

}

}

m <- mosaic(l[[1]], l[[2]], l[[3]], l[[4]], l[[5]], l[[6]], l[[7]], l[[8]], l[[9]], l[[10]], l[[11]],l[[12]], l[[11]], l[[12]], l[[13]], l[[14]], l[[15]], l[[16]], l[[17]], l[[18]], l[[19]], fun=median)

m <- m*100
m[is.infinite(m)==TRUE] <- NA
m[m>100] <- 100

par(mar=c(6,6,6,6))
par(mfrow=c(1,2))

plot(e, main='', col='gray100')
plot(mer, add=TRUE, legend=FALSE, col=dichromat(rev(terrain.colors(15))))
plot(mer,zlim=c(0, 0.5), add=TRUE, legend=FALSE, col=dichromat(rev(terrain.colors(15))))
plot(mer,zlim=c(0.5, 1), add=TRUE, legend=TRUE, col=dichromat(rev(terrain.colors(15))))
plot(bPols[-55], add=TRUE, cex=2)


plot(e, main='', col='gray100')
plot(m, add=TRUE,
legend=TRUE, col=rev(heat.colors(15)))
plot(bPols[-55], add=TRUE, cex=2)


plot(bPols[c(25, 29, 40)],add=TRUE, col='black')



c <- weighted.mean(stack(myRaster,  projectRaster(raster(mergeParallelX["var1.pred"]), myRaster)),  w=c(1, 3), na.rm=FALSE)




plot(e, main='', col='gray13')
plot(sc,  add=TRUE,  col='black', legend=FALSE)

for(i in 1:length(l)){

sc <- merge(sc, l[[i]])

}


lim <- quantile(l[[i]], probs = c(0.025, 0.975))
plot(l[[i]],  add=TRUE, legend=FALSE, col=rev(heat.colors(25)), zlim=c(lim[1], lim[2]))
  
  
}
  
library(maps)
map('world', add=TRUE, col='gray', cex=0.5)

> best method?

  country    best-R  best-RMSE
1  Argentina  RK     RK
2  Belize     RF     RK
3  Bolivia    SVM    KK
4  Brazil     RF     RF
5  Chile      PL     PL
6  Colombia   RF     RF
7  Costa Rica SVM    SVM
8  Cuba       PL     PL
9  Ecuador    RK     RK
10 Guatemala  KK     RF
11 Honduras   SVM    KK
12 Jamaica    RF     RF
13 Mexico     RK     RK
14 Nicaragua  RF     RF
15 Panama     PL     KK
16 Peru       KK     KK
17 Suriname   SVM    PL
18 Uruguay    RF     RK
19 Venezuela  RK     RK


res <- data.frame(country=character(), 
                  nPix=numeric(),
                  sumPix=numeric(),
                  sumPixUNC=numeric())
res$country <- as.character(res$country)
for(i in 1:length(l)){
country <- strsplit(lis[i], '_')[[1]][1]
s <- mask(resample(mer, stack(lis[[i]])), stack(lis[i])[[1]])
u <-  mask(resample(m, stack(lis[[i]])), stack(lis[[i]])[[1]])
  res[i, 1] <- country
  res[i, 2] <- dim(na.omit(as.data.frame(s)))[1]
  res[i, 3] <- round(cellStats(s, sum), 2)
  res[i, 4] <- round(cellStats(u, sum), 2)
print(country)  
}

res[3, 1] <- 'Bolivia'
res[19, 1] <- 'Venezuela'

###
###
	par(mar=c(5,5,5,5))
	df <- data.frame(x=log(res$nPix), y=log(res$sumPix))
	mod <- lm(y ~ x, data = df)
	newx <- seq(min(df$x), max(df$x), length.out=100)
	preds <- predict(mod, newdata = data.frame(x=newx), 
		         interval = 'confidence')

#Create a function to generate a continuous color palette


color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=19) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
x <- c((1:19)^2, (19:1)^2)


	plot(y ~ x, data = df, type = 'n',xlab = 'AREA log (n 5x5km pixels per country) ', ylab = 'STOCK log (sum of 5x5km pixels of SOC predicted values)', cex.axis=2, main = '', cex.lab=1.5, col=topo.colors(19))

	polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey95', border = NA)
	points(y~x, pch=20, cex=4, data=df, col=color.gradient(x))

	abline(mod, lty=2)

newdata <- res[order(res$nPix),] 
newdata
x <- c((19:1)^2, (1:19)^2)
legend("bottomright",legend=rev(newdata$country),col =color.gradient(x), pch=20, cex=1.3,pt.cex=1.5, bty='n')

#with(res, text(log(nPix)~log(sumPix), labels = res$country, pos = 3, cex=1.3))


r2 <- summary(mod)[9]
r2 <- as.character(round(unlist(r2), 2))

legend ('bottom', paste0( 'R-squared ', r2), cex=1.3, bty='n')

###
###
res$socArea <- res$sumPix/res$nPix
res$stdUnc <- (res$sumPixUNC/res$sumPix)

par(mar=c(5,5,5,5))
df <- data.frame(x=res$socArea, y=res$stdUnc)
mod <- lm(y ~ x, data = df)
newx <- seq(min(df$x), max(df$x), length.out=100)
preds <- predict(mod, newdata = data.frame(x=newx), 
                 interval = 'confidence')


plot(y ~ x, data = df, type = 'n',xlab='STOCK log (sum of 5x5km pixels of SOC predicted values) / AREA', ylab = 'standardized UNCERTAINTY (%)', cex.axis=2, main = '', cex.lab=1.5)

polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey95', border = NA)
points(y~x, pch=20, cex=4, data=df, col=color.gradient(x))

	abline(mod, lty=2)

#with(res, text(log(uncArea)~log(socArea), labels = res$country, pos = 3, cex=1))
newdata <- res[order(res$socArea),] 
newdata
x <- c((19:1)^2, (1:19)^2)
legend("topright",legend=rev(newdata$country),col =color.gradient(x), pch=20, cex=1.3,pt.cex=1.5, bty='n')


r2 <- summary(mod)[9]
r2 <- as.character(round(unlist(r2), 1))

legend ('bottom', paste0( 'R-squared ', r2),cex=1.3, bty='n')



