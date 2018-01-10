## Establecemos el directorio de trabajo
setwd("~/Proyectos/FAO-RSPCA-Course/")


## Load required packages
library(raster)
library(car)
library(rgdal)
library(gstat)

#### Load covariates #####

#### Covariates preparation ####

COV <- stack("extra/COV1km.tif")
names(COV)
names(COV)[[1]]='temp'
names(COV)[[2]]='soil_sedimentary_deposit_thickness'
names(COV)[[3]]="BDRICM_M_250m_ll"
names(COV)[[4]]="BLDFIE_M_sl1_250m_ll"
names(COV)[[5]]="BLDFIE_M_sl2_250m_ll"
names(COV)[[6]]="BLDFIE_M_sl3_250m_ll"
names(COV)[[7]]="CLYPPT_M_sl1_250m_ll"  
names(COV)[[8]]="CLYPPT_M_sl2_250m_ll"  
names(COV)[[9]]="CLYPPT_M_sl3_250m_ll"  
names(COV)[[10]]="CLYPPT_M_sl4_250m_ll"  
names(COV)[[11]]="CLYPPT_M_sl5_250m_ll"
names(COV)[[12]]="CRFVOL_M_sl1_250m_ll"
names(COV)[[13]]="CRFVOL_M_sl2_250m_ll"
names(COV)[[14]]="CRFVOL_M_sl3_250m_ll"
names(COV)[[15]]="DEMSRE3a"
names(COV)[[16]]='ESA_CCI_land_useBAD'
names(COV)[[17]]='EVMMOD3a'
names(COV)[[18]]='GEAISG3a'
names(COV)[[19]]='HISTPR_250m_ll'
names(COV)[[20]]='OCSTHA_M_30cm_250m_ll'
names(COV)[[21]]='prec'
names(COV)[[22]]='STGHWS1aBAD'
names(COV)[[23]]='tempBAD'
names(COV)[[24]]='upland_valleybottom_lowland_sedimentary_deposit_thickness'
names(COV)[[25]]='ESA_CCI_land_use'
names(COV)[[26]]='STGHWS1a'

names(COV)
COV=COV[[-c(16, 22, 23)]]
names(COV)

grd.lst <- list.files(path = "DEMcov/", pattern="sdat$")
COVdem <- stack(paste0("DEMcov/", grd.lst))
COV <- stack(COV, COVdem)


COV=COV[[-c(29,30, 6,9,11,14)]]

### Preparamos los datos ####


datasets.lst <- list.files("~/Dropbox/SOIL - Central America/Data/", pattern = "csv$")

results <- data.frame(country = character(), top5v = character(), 
                      top3v = character(), r2total = numeric(),
                      r2top5 = numeric(), r2top3 = numeric())




for(i in 1:length(datasets.lst)){
  dat=read.csv(paste0("~/Dropbox/SOIL - Central America/Data/", datasets.lst[i]))
  if(nrow(dat) < 15){next()}
  country_name <- strsplit(datasets.lst[i], "_")[[1]][1]
print(country_name)
  ### Marcamos 
dat$ESA_CCI_land_use <- as.factor(dat$ESA_CCI_land_use)
dat$STGHWS1a <- as.factor(dat$STGHWS1a)

names(dat)
#dat <- dat[,-c(27:28, 34:35)]
dat <- dat[,-c(9,11,12,14,15,16,17,19,34,35,28,27)]
str(dat)

hist(dat$OCSKGM30, breaks=100)

# # Transform data

## Recreamos el objeto con la ubicacion de los puntos
dat_sp <- dat
coordinates(dat_sp) <- ~ longitude + latitude

## Ajustamos un modelo de regresion lineal multiple
## Pruebas de modelos ####

## Ahora, a diferencia del ejercicio 4, el modelo lo ajustamos solo con datos.model
## y no con todos los datos...!
modelo.MLR <- lm(log(OCSKGM30) ~ . -id-X, data = dat_sp@data) 

summary(modelo.MLR)
anova(modelo.MLR)


## Hacemos seleccion de variables por stepwise
modelo.MLR.step <- step(modelo.MLR, direction="both")

summary(modelo.MLR.step)
anova <- anova(modelo.MLR.step)
variables <- names(modelo.MLR.step$coefficients)
variables <- variables[2:length(variables)]
variables <- data.frame(variable=variables, 
                        partialr2=anova$`Sum Sq`[1:length(anova$`Sum Sq`)-1] / sum(anova$`Sum Sq`))

variables <- variables[order(variables$partialr2, decreasing = T),]
top5 <- head(variables, 5)
top3 <- head(variables, 3)

top5v <- paste(as.character(top5$variable), collapse = " + ")
top3v <- paste(as.character(top3$variable), collapse = " + ")

results <- rbind(results, cbind(country_name, top5v, top3v, sum(variables$partialr2), 
                            sum(top5$partialr2), sum(top3$partialr)))

}

boxplot(as.numeric(as.character(results$V4)), as.numeric(as.character(results$V5)), as.numeric(as.character(results$V6)))

write.csv(results, file = "~/Dropbox/SOIL - Central America/Results/topvariables.csv")
