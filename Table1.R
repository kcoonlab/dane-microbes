#Variance partitioning 

#Spatial stuff

library(Imap)

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}
GeoDistanceInMetresMatrix <- function(df.geopoints){  
  GeoDistanceInMetres <- function(g1, g2){
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$latitude, lon.1=g1$longitude, lat.2=g2$latitude, lon.2=g2$longitude, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  n.geopoints <- nrow(df.geopoints)  
  df.geopoints$index <- 1:n.geopoints
  list.geopoints <- by(df.geopoints[,c("index", "latitude", "longitude")], 1:n.geopoints, function(x){return(list(x))})
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  return(mat.distances)
}
coords <- read.csv("coords.csv")
coords$id <- as.character(coords$id)
df.cities <- coords
coords2 <- GeoDistanceInMetresMatrix(df.cities)
space <- pcnm(coords2)
space2 <- space$vectors
write.csv(space2,"pcnm-vectors.csv",row.names = coords[,1])

assem_bac<-qpen(comm = bac, pd = phydist, pd.big.wd = NULL,
                  pd.big.spname = NULL, tree = NULL,
                  bNTI = NULL, RC = NULL, ab.weight = TRUE,
                  meta.ab = NULL, exclude.conspecifics = FALSE,
                  rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
                  nworker = 4, memory.G = 50, project = NA,
                  wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)
bac_assem_ratio<-data.frame(assem_bac[["ratio"]])
bac_assem_result<-data.frame(assem_bac[["result"]])
bac_assem_result['Biotype']= 'Both'

library(gtools)
bac_bnti <- bac_assem_result[,c(1,2,5)]
multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
    do.call(order, c(
        lapply(list(...), function(l){
            if(is.character(l)){
                factor(l, levels=mixedsort(unique(l)))
            } else {
                factor(as.character(l), levels=mixedsort(levels(l)))
            }
        }),
        list(na.last = na.last, decreasing = decreasing)
    ))
}
bac_bnti <- bac_bnti[do.call(multi.mixedorder, c(bac_bnti[1:2], list(decreasing = FALSE))),]
bac_bnti$sample1 <- factor(bac_bnti$sample1, levels = c("S5","S6","S13","S26","S81","S83","S86","S88","S158","S164","S189","S218","S243","S245","S247","S249","S251","S253","S254","S256","S319","S335","S340","S369","S371","S391","S397","S417","S465","S513","S515","S522","S523","S526","S530","S552","S556","S558","S559","S565","S567","S583","S584","S587","S590","S597","S598","S607","S615","S618","S619","S620","S676","S700","S715","S718","S719","S724","S725","S900","S969","S972","S973","S975","S978","S1931","S1932","S2041","S3262","S3333","S3921","S3922","S4022","S4130","S8300","S8301","S8512","S9007","S9009","S9014","S9019","S9098","S9601","S9970","S9971"))

library(tidyr)
bac_bnti  <- spread(bac_bnti , sample1, bNTI)
bac_bnti <- bac_bnti[do.call(multi.mixedorder, c(bac_bnti[1], list(decreasing = FALSE))),]
bac_bnti <- data.frame(bac_bnti, row.names = 1,check.names=FALSE)
bac_bnti <- as.matrix(bac_bnti)
bac_bnti <- cbind(S4 = NA, bac_bnti)
bac_bnti <- rbind(bac_bnti,S9971 = NA )
diag(bac_bnti) <- 0

library(Matrix)
bac_bnti<- as.matrix(forceSymmetric(bac_bnti))
as.dist(bac_bnti)
bac_bnti1 <- decostand(bac_bnti, method="range")

prod<-read.csv("hist-data.csv")
prod_vp<-prod[,c(4:17)]

library(caret)
prod.impute<-preProcess(prod_vp,method="bagImpute")
prod_vp2<-predict(prod.impute,prod_vp)

env<-read.csv("hist-env-data.csv")
env_vp<-env[,c(5:10)]
env_vp$PercSiteDry <- prod$PercSiteDry
env.impute<-preProcess(env_vp,method="bagImpute")
env_vp2<-predict(env.impute,env_vp)

library(GGally)
env_vp3 <- data.frame(scale(env_vp2))
ggpairs(env_vp3)

cap <- capscale(as.dist(bac_bnti1) ~ .,data=env_vp3)
ord <- ordistep(cap) #AvgLSTDay (p = 0.01)
anova(ord) #p = 0.011
anova.cca(ord,by="term") #p = 0.009

prod_vp3 <- data.frame(scale(prod_vp2))
cap <- capscale(as.dist(bac_bnti1) ~ .,data=prod_vp3)
ord <- ordistep(cap) #PercAeFnd (p = 0.050), PercVecFnd (p = 0.035), AvgVecDens (p = 0.035), TotTreat (p = 0.020), AvgAnDens (p = 0.005)
anova(ord) #p = 0.001
anova.cca(ord,by="term") #TotTreat (p = 0.001), PercAeFnd (p = 0.009), PercVecFnd (p = 0.382), AvgAnDens (p = 0.001), and AvgVecDens (p = 0.017) 

cap <- capscale(as.dist(bac_bnti1) ~ .,data=as.data.frame(space2))
ord <- ordistep(cap) #PCNM30 (p = 0.075), PCNM3 (p = 0.075), PCNM33 (p = 0.065), PCNM11 (p = 0.030)
anova(ord) #p = 0.005
anova.cca(ord, by="term") #PCNM3 (p = 0.050), PCNM11 (p = 0.050), PCNM30 (p = 0.085), PCNM33 (p = 0.080) 

env_vp4 <- env_vp3[,1]
prod_vp4 <- prod_vp3[,c(1,4,8,11,14)]
space3 <- space2[,c(3,11,30,33)]
varpart_16s <- varpart(as.dist(bac_bnti1),env_vp4,prod_vp4,space3)
varpart_16s #environmental variables do not contribute individually to variance; discard
varpart_16s <- varpart(as.dist(bac_bnti1),prod_vp4,space3)
varpart_16s
plot(varpart_16s) #prod = 0.23246, space = 0.08019, prod + space = 0.27746, prod | space = 0.19728, space | prod = 0.04501, prod int space = 0.03518, Residuals = 0.72254 

ord1 <- capscale(as.dist(bac_bnti1)~as.matrix(prod_vp4))
anova(ord1) #p = 0.001

ord2 <- capscale(as.dist(bac_bnti1)~as.matrix(prod_vp4)+Condition(as.matrix(space3)))
anova(ord2) #p = 0.001

ord3 <- capscale(as.dist(bac_bnti1)~as.matrix(space3))
anova(ord3) #p = 0.007

ord4 <- capscale(as.dist(bac_bnti1)~as.matrix(space3)+Condition(as.matrix(prod_vp4)))
anova(ord4) #p = 0.013

ord5 <- capscale(as.dist(bac_bnti1)~as.matrix(prod_vp4)+as.matrix(space3))
anova(ord5) #p = 0.001

vars <- cbind(prod_vp4, space3)
ord6 <- capscale(as.dist(bac_bnti1)~vars$TotTreat+vars$PercAeFnd+vars$PercVecFnd+vars$AvgAnDens+vars$AvgVecDens+Condition(vars$PCNM3+vars$PCNM11+vars$PCNM30+vars$PCNM33))
anova(ord6) #p = 0.001
anova.cca(ord6, by="term") #TotTreat (p = 0.005), PercAeFnd (p = 0.024), PercVecFnd (p = 0.304), AvgAnDens (p = 0.001), AvgVecDens (p = 0.019) 

ord7 <- capscale(as.dist(bac_bnti1)~vars$PCNM3+vars$PCNM11+vars$PCNM30+vars$PCNM33+Condition(vars$TotTreat+vars$PercAeFnd+vars$PercVecFnd+vars$AvgAnDens+vars$AvgVecDens))
anova(ord7) #p = 0.013
anova.cca(ord7, by="term") #PCNM3 (p = 0.129), PCNM11 (p = 0.079), PCNM30 (p = 0.171), PCNM33 (p = 0.037)

a <- plot(x=coords$lat, y=coords$lon, data=coords)
ordisurf(coords[,c(2:3)],space2[,33],bubble=4,main='PCNM 33')

coords.new=coords[,c(1,3,2)]
require(sp)
tempsf <- coords[, 2:3]
coordinates(tempsf) <- c("latitude","longitude")
proj4string(tempsf) = "+proj=longlat +ellps=WGS84 +no_defs"
max(spDists(tempsf))
min(spDists(tempsf))
dist <- spDists(as.matrix(coords.new[,2:3]), longlat = TRUE)
write.csv(dist,"pairwise-dist.csv",row.names = coords.new[,1])
write.csv(assem_bac$result,"pairwise-bnti.csv",row.names = FALSE)













#Variance partitioning 

#Spatial stuff

library(Imap)

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}
GeoDistanceInMetresMatrix <- function(df.geopoints){  
  GeoDistanceInMetres <- function(g1, g2){
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$latitude, lon.1=g1$longitude, lat.2=g2$latitude, lon.2=g2$longitude, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  n.geopoints <- nrow(df.geopoints)  
  df.geopoints$index <- 1:n.geopoints
  list.geopoints <- by(df.geopoints[,c("index", "latitude", "longitude")], 1:n.geopoints, function(x){return(list(x))})
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  return(mat.distances)
}
coords <- read.csv("coords.csv")
coords$id <- as.character(coords$id)
df.cities <- coords
coords2 <- GeoDistanceInMetresMatrix(df.cities)
space <- pcnm(coords2)
space2 <- space$vectors
write.csv(space2,"pcnm-vectors.csv",row.names = coords[,1])

assem_bac<-qpen(comm = bac, pd = phydist, pd.big.wd = NULL,
                  pd.big.spname = NULL, tree = NULL,
                  bNTI = NULL, RC = NULL, ab.weight = TRUE,
                  meta.ab = NULL, exclude.conspecifics = FALSE,
                  rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
                  nworker = 4, memory.G = 50, project = NA,
                  wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)
bac_assem_ratio<-data.frame(assem_bac[["ratio"]])
bac_assem_result<-data.frame(assem_bac[["result"]])
bac_assem_result['Biotype']= 'Both'

library(gtools)
bac_bnti <- bac_assem_result[,c(1,2,5)]
multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
    do.call(order, c(
        lapply(list(...), function(l){
            if(is.character(l)){
                factor(l, levels=mixedsort(unique(l)))
            } else {
                factor(as.character(l), levels=mixedsort(levels(l)))
            }
        }),
        list(na.last = na.last, decreasing = decreasing)
    ))
}
bac_bnti <- bac_bnti[do.call(multi.mixedorder, c(bac_bnti[1:2], list(decreasing = FALSE))),]
bac_bnti$sample1 <- factor(bac_bnti$sample1, levels = c("S5","S6","S13","S26","S81","S83","S86","S88","S158","S164","S189","S218","S243","S245","S247","S249","S251","S253","S254","S256","S319","S335","S340","S369","S371","S391","S397","S417","S465","S513","S515","S522","S523","S526","S530","S552","S556","S558","S559","S565","S567","S583","S584","S587","S590","S597","S598","S607","S615","S618","S619","S620","S676","S700","S715","S718","S719","S724","S725","S900","S969","S972","S973","S975","S978","S1931","S1932","S2041","S3262","S3333","S3921","S3922","S4022","S4130","S8300","S8301","S8512","S9007","S9009","S9014","S9019","S9098","S9601","S9970","S9971"))

library(tidyr)
bac_bnti  <- spread(bac_bnti , sample1, bNTI)
bac_bnti <- bac_bnti[do.call(multi.mixedorder, c(bac_bnti[1], list(decreasing = FALSE))),]
bac_bnti <- data.frame(bac_bnti, row.names = 1,check.names=FALSE)
bac_bnti <- as.matrix(bac_bnti)
bac_bnti <- cbind(S4 = NA, bac_bnti)
bac_bnti <- rbind(bac_bnti,S9971 = NA )
diag(bac_bnti) <- 0

library(Matrix)
bac_bnti<- as.matrix(forceSymmetric(bac_bnti))
as.dist(bac_bnti)
bac_bnti1 <- decostand(bac_bnti, method="range")

prod<-read.csv("hist-data.csv")
prod_vp<-prod[,c(4:17)]

library(caret)
prod.impute<-preProcess(prod_vp,method="bagImpute")
prod_vp2<-predict(prod.impute,prod_vp)

env<-read.csv("hist-env-data.csv")
env_vp<-env[,c(5:10)]
env_vp$PercSiteDry <- prod$PercSiteDry
env.impute<-preProcess(env_vp,method="bagImpute")
env_vp2<-predict(env.impute,env_vp)

library(GGally)
env_vp3 <- data.frame(scale(env_vp2))

cap <- capscale(as.dist(bac_bnti1) ~ .,data=env_vp3)
ord <- ordistep(cap) #AvgLSTDay
anova(ord) #p = 0.016

prod_vp3 <- data.frame(scale(prod_vp2))
cap <- capscale(as.dist(bac_bnti1) ~ .,data=prod_vp3)
ord <- ordistep(cap) #AvgCxDens, AvgVecDens, AvgNonVecDens, PercCxFnd, PercVecFnd, TotTreat, PercAeFnd, AvgAnDens
anova(ord) #p = 0.001
vif.cca(ord) #All <10
anova.cca(ord,by="margin") #TotTreat (p = 0.094), PercAeFnd (p = 0.028), PercCxFnd (p = 0.126), PercVecFnd (p = 0.098), AvgAnDens (p = 0.003), AvgCxDens (p = 0.112), AvgNonVecDens (p = 0.111), and AvgVecDens (p = 0.111) 

cap <- capscale(as.dist(bac_bnti1) ~ .,data=as.data.frame(space2))
ord <- ordistep(cap)
anova(ord)
anova.cca(ord, by="margin") #PCNM2 (p = 0.077), PCNM3 (p = 0.041), PCNM4 (p = 0.080), PCNM11 (p = 0.052), PCNM30 (p = 0.064), PCNM33 (p = 0.070) 

env_vp4 <- env_vp3[,1]
prod_vp4 <- prod_vp3[,c(1,4,6,8,11:14)]
space3 <- space2[,c(2:4,11,30,33)]
varpart_16s <- varpart(as.dist(bac_bnti1),env_vp4,prod_vp4,space3)
varpart_16s #environmental variables do not contribute individually to variance; discard
varpart_16s <- varpart(as.dist(bac_bnti1),prod_vp4,space3)
varpart_16s
plot(varpart_16s) #prod = 0.21471, space = 0.11218, prod + space = 0.26930, prod | space = 0.15712, space | prod = 0.05459, prod int space = 0.05758, Residuals = 0.73070 

ord1 <- capscale(as.dist(bac_bnti1)~as.matrix(prod_vp4))
anova(ord1) #p = 0.001

ord2 <- capscale(as.dist(bac_bnti1)~as.matrix(prod_vp4)+Condition(as.matrix(space3)))
anova(ord2) #p = 0.001

ord3 <- capscale(as.dist(bac_bnti1)~as.matrix(space3))
anova(ord3) #p = 0.002

ord4 <- capscale(as.dist(bac_bnti1)~as.matrix(space3)+Condition(as.matrix(prod_vp4)))
anova(ord4) #p = 0.013

ord5 <- capscale(as.dist(bac_bnti1)~as.matrix(prod_vp4)+as.matrix(space3))
anova(ord5) #p = 0.001

vars <- cbind(prod_vp4, space3)
ord6 <- capscale(as.dist(bac_bnti1)~vars$TotTreat+vars$PercAeFnd+vars$PercCxFnd+vars$PercVecFnd+vars$AvgAnDens+vars$AvgCxDens+vars$AvgNonVecDens+vars$AvgVecDens+Condition(vars$PCNM2+vars$PCNM3+vars$PCNM4+vars$PCNM11+vars$PCNM30+vars$PCNM33))
anova(ord6) #p = 0.001
anova.cca(ord6, by="margin") #AvgAnDens (p = 0.001)

ord7 <- capscale(as.dist(bac_bnti1)~vars$PCNM2+vars$PCNM3+vars$PCNM4+vars$PCNM11+vars$PCNM30+vars$PCNM33+Condition(vars$TotTreat+vars$PercAeFnd+vars$PercCxFnd+vars$PercVecFnd+vars$AvgAnDens+vars$AvgCxDens+vars$AvgNonVecDens+vars$AvgVecDens))
anova(ord7) #p = 0.01
anova.cca(ord7, by="margin") #PCNM4 ( p = 0.050), PCNM33 (p = 0.038)

a <- plot(x=coords$lat, y=coords$lon, data=coords)
ordisurf(coords[,c(2:3)],space2[,33],bubble=4,main='PCNM 33')

coords.new=coords[,c(1,3,2)]
require(sp)
tempsf <- coords[, 2:3]
coordinates(tempsf) <- c("latitude","longitude")
proj4string(tempsf) = "+proj=longlat +ellps=WGS84 +no_defs"
max(spDists(tempsf))
min(spDists(tempsf))
dist <- spDists(as.matrix(coords.new[,2:3]), longlat = TRUE)
write.csv(dist,"pairwise-dist.csv",row.names = coords.new[,1])
write.csv(assem_bac$result,"pairwise-bnti.csv",row.names = FALSE)







































prod<-read.csv("hist-data.csv")
prod_vp<-prod[,c(4:17)]

library(caret)
prod.impute<-preProcess(prod_vp,method="bagImpute")
prod_vp2<-predict(prod.impute,prod_vp)

env<-read.csv("hist-env-data.csv")
env_vp<-env[,c(5:10)]
env_vp$PercSiteDry <- prod$PercSiteDry
env.impute<-preProcess(env_vp,method="bagImpute")
env_vp2<-predict(env.impute,env_vp)

library(GGally)
env_vp3 <- data.frame(scale(env_vp2))
#ggpairs(env_vp3)
cor2pcor(cov(env_vp3))

library(corpcor)
cor2pcor(cov(env_vp3)) #AvgEVI and AvgNDVI are colinear; proceed with AvgEVI only

env_vp3 <- data.frame(env_vp2[,c(1:3,5:7)])
cap <- capscale(as.dist(bac_bnti1) ~ .,data=env_vp3)
ord <- ordistep(cap) #AvgLSTDay
anova(ord) #p = 0.01

prod_vp3 <- data.frame(scale(prod_vp2))
#ggpairs(prod_vp3)
cor2pcor(cov(prod_vp3)) #PercCxFnd, PercVecFnd, AvgLarvDens, AvgAeDens, AvgCxDens, AvgNonVecDens, AvgVecDens show evidence of multicollinearity

prod_vp3 <- data.frame(prod_vp3[,c(1:5,8,10,11,14)])
cap <- capscale(as.dist(bac_bnti1) ~ .,data=prod_vp3)
ord <- ordistep(cap) #TotTreat, PercTreat, AvgVecDens, PercAeFnd, AvgAnDens
anova(ord) #p = 0.001
anova.cca(ord,by="term") #TotTreat (p = 0.020), PercTreat (p = 0.018), PercAeFnd (p = 0.006), AvgAnDens (p = 0.001), AvgVecDens (p = 0.007)

cap <- capscale(as.dist(bac_bnti1) ~ .,data=as.data.frame(space2))
ord <- ordistep(cap) #PCNM2, 3, 4, 11, 30
anova(ord) #p = 0.003
anova.cca(ord, by="margin") #PCNM3 (p = 0.031), PCNM11 (p = 0.033) 

env_vp4 <- env_vp3[,1]
prod_vp4 <- prod_vp3[,c(1,2,9,4,5)]
space3 <- space2[,c(3,4,11,30)]
varpart_16s <- varpart(as.dist(bac_bnti1),env_vp4,prod_vp4,space3)
varpart_16s
plot(varpart_16s)

ord1 <- capscale(as.dist(bac_bnti1)~as.matrix(prod_vp4))
anova(ord1) #p = 0.001

ord2 <- capscale(as.dist(bac_bnti1)~as.matrix(prod_vp4)+Condition(as.matrix(space3)))
anova(ord2) #p = 0.001

ord3 <- capscale(as.dist(bac_bnti1)~as.matrix(space3))
anova(ord3) #p = 0.001

ord4 <- capscale(as.dist(bac_bnti1)~as.matrix(space3)+Condition(as.matrix(prod_vp4)))
anova(ord4) #p = 0.001

ord5 <- capscale(as.dist(bac_bnti1)~as.matrix(prod_vp4)+as.matrix(space3))
anova(ord5) #p = 0.001

vars <- cbind(prod_vp4, space3)
ord6 <- capscale(as.dist(bac_bnti1)~vars$TotTreat+vars$PercTreat+vars$PercAeFnd+vars$AvgAnDens+vars$AvgVecDens+Condition(vars$PCNM2+vars$PCNM3+vars$PCNM4+vars$PCNM11+vars$PCNM17+vars$PCNM30+vars$PCNM33))
anova(ord6) #p = 0.001
anova.cca(ord6, by="term") #TotTreat (p = 0.001), PercTreat (p = 0.290), PercAeFnd (p = 0.174), AvgAnDens (p = 0.001), AvgVecDens (p = 0.003) 

ord7 <- capscale(as.dist(bac_bnti1)~vars$PCNM2+vars$PCNM3+vars$PCNM4+vars$PCNM11+vars$PCNM17+vars$PCNM30+vars$PCNM33+Condition(vars$TotTreat+vars$PercTreat+vars$PercAeFnd+vars$AvgAnDens+vars$AvgVecDens))
anova(ord7) #p = 0.003
anova.cca(ord7, by="term") #PCNM2 (p = 0.208), PCNM3 (p = 0.090), PCNM4 (p = 0.022), PCNM11 (p = 0.057), PCNM17 (p = 0.205), PCNM30 (p = 0.212), PCNM33 (p = 0.043)

a <- plot(x=coords$lat, y=coords$lon, data=coords)
ordisurf(coords[,c(2:3)],space2[,4],bubble=4,main='PCNM 4')
ordisurf(coords[,c(2:3)],space2[,33],bubble=4,main='PCNM 33')

coords.new=coords[,c(1,3,2)]
require(sp)
tempsf <- coords[, 2:3]
coordinates(tempsf) <- c("latitude","longitude")
proj4string(tempsf) = "+proj=longlat +ellps=WGS84 +no_defs"
max(spDists(tempsf))
min(spDists(tempsf))
dist <- spDists(as.matrix(coords.new[,2:3]), longlat = TRUE)
write.csv(dist,"pairwise-dist.csv",row.names = coords.new[,1])
write.csv(assem_bac$result,"pairwise-bnti.csv",row.names = FALSE)
