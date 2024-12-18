set.seed(123)

#url <- "https://cran.r-project.org/src/contrib/Archive/Imap/Imap_1.32.tar.gz"
#pkgFile <- "Imap_1.32.tar.gz"
#download.file(url = url, destfile = pkgFile)
# Expand the zip file using whatever system functions are preferred
# Install package
#install.packages(pkgs=pkgFile, type="source", repos=NULL)
# Delete package tarball
#unlink(pkgFile)

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("biomformat")

library(Imap)
library(gtools)
library(tidyr)
library(Matrix)
library(caret)
library(GGally)
library(vegan)
library(iCAMP)
library(biomformat)
library(ape)
library(corpcor)

#Variance partitioning 

#Spatial stuff

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
#write.csv(space2,"pcnm-vectors.csv",row.names = coords[,1])

Bac_otus<-"exported/feature-table.biom"
x1<-read_biom(Bac_otus)
meta<-read.csv("metadata3.csv")
bac<-as(biom_data(x1),"matrix")
bac<-t(bac)
tree<-read.tree("exported-tree/tree.tree")
phydist<-cophenetic(tree)
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

bac_bnti  <- spread(bac_bnti , sample1, bNTI)
bac_bnti <- bac_bnti[do.call(multi.mixedorder, c(bac_bnti[1], list(decreasing = FALSE))),]
bac_bnti <- data.frame(bac_bnti, row.names = 1,check.names=FALSE)
bac_bnti <- as.matrix(bac_bnti)
bac_bnti <- cbind(S4 = NA, bac_bnti)
bac_bnti <- rbind(bac_bnti,S9971 = NA )
diag(bac_bnti) <- 0

bac_bnti<- as.matrix(forceSymmetric(bac_bnti))
as.dist(bac_bnti)
bac_bnti1 <- decostand(bac_bnti, method="range")

prod<-read.csv("GitHub-files/hist-data.csv")
prod_vp<-prod[,c(5,6,9:11,15:17)]

#prod.impute<-preProcess(prod_vp,method="bagImpute")
#prod_vp2<-predict(prod.impute,prod_vp)

env<-read.csv("hist-env-data.csv")
env_vp<-env[,c(5:10)]
env_vp$PercSiteDry <- prod$PercSiteDry
env.impute<-preProcess(env_vp,method="bagImpute")
env_vp2<-predict(env.impute,env_vp)

env_vp3 <- data.frame(scale(env_vp2))
#env_vp3 <- decostand(env_vp2, method = "standardize")
cor2pcor(cov(env_vp3)) #AvgEVI and AvgNDVI are colinear; proceed with AvgEVI only

env_vp4 <- data.frame(env_vp3[,c(1:3,5:7)])
cap <- capscale(as.dist(bac_bnti1) ~ .,data=env_vp4)
ord <- ordistep(cap) #AvgLSTDay
anova(ord) #p = 0.013

prod_vp2 <- prod_vp[,-1]
prod_vp3 <- data.frame(scale(prod_vp2))
#prod_vp3 <- decostand(prod_vp2, method = "standardize")
cor2pcor(cov(prod_vp3)) #PercAnFnd and AvgAnDens are colinear; proceed with PercAnFnd only
prod_vp4 <- prod_vp3[,-3]

cap <- capscale(as.dist(bac_bnti1) ~ .,data=prod_vp4)
ord <- ordistep(cap) #TotTreat, PercAeFnd, PercCxFnd, AvgAnDens, AvgCxDens
anova(ord) #p = 0.001
anova.cca(ord,by="term") #TotTreat (p = 0.002), PercAeFnd (p = 0.007), PercCxFnd (p = 0.508), AvgAnDens (p = 0.001), AvgCxDens (0.008)

cap <- capscale(as.dist(bac_bnti1) ~ .,data=as.data.frame(space2))
ord <- ordistep(cap) #PCNM2, 3, 4, 11, 17, 30, 33
anova(ord) #p = 0.001
anova.cca(ord, by="margin") #PCNM2 (p = 0.070), PCNM3 (p = 0.057), PCNM4 (p = 0.069), PCNM11 (p = 0.034), PCNM17 (p = 0.113), PCNM30 (p = 0.089), PCNM33 (p = 0.061) 

env_vp5 <- env_vp4[,1]
prod_vp5 <- prod_vp4[,c(1:3,5,6)]
space3 <- space2[,c(2:4,11,17,30,33)]
varpart_16s <- varpart(as.dist(bac_bnti1),env_vp5,prod_vp5,space3)
varpart_16s
plot(varpart_16s)

varpart_16s <- varpart(as.dist(bac_bnti1),prod_vp5,space3)
varpart_16s
plot(varpart_16s)

vars <- cbind(prod_vp5, space3)
ord1 <- capscale(as.dist(bac_bnti1)~vars$TotTreat+vars$PercAeFnd+vars$PercCxFnd+vars$AvgAnDens+vars$AvgCxDens)
anova(ord1) #p = 0.001
anova.cca(ord1, by="term") #TotTreat (p = 0.005), PercAeFnd (p = 0.010), PercCxFnd (p = 0.511), AvgAnDens (p = 0.001), AvgCxDens (p = 0.019) 

ord2 <- capscale(as.dist(bac_bnti1)~vars$PCNM2+vars$PCNM3+vars$PCNM4+vars$PCNM11+vars$PCNM17+vars$PCNM30+vars$PCNM33)
anova(ord2) #p = 0.001
anova.cca(ord2, by="term")  #PCNM2 (p = 0.088), PCNM3 (p = 0.036), PCNM4 (p = 0.058), PCNM11 (p = 0.041), PCNM17 (p = 0.103), PCNM30 (p = 0.055), PCNM33 (p = 0.062) 

ord3 <- capscale(as.dist(bac_bnti1)~vars$TotTreat+vars$PercAeFnd+vars$PercCxFnd+vars$AvgAnDens+vars$AvgCxDens+vars$PCNM2+vars$PCNM3+vars$PCNM4+vars$PCNM11+vars$PCNM17+vars$PCNM30+vars$PCNM33)
anova(ord3) #p = 0.001
anova.cca(ord3, by="term") #TotTreat (p = 0.001), PercAeFnd (p = 0.008), PercCxFnd (p = 0.484), AvgAnDens (p = 0.001), AvgCxDens (p = 0.008), #PCNM2 (p = 0.206), #PCNM3 (p = 0.102), #PCNM4 (p = 0.072), #PCNM11 (p = 0.104), #PCNM17 (p = 0.190), #PCNM30 (p = 0.192), #PCNM33 (p = 0.035)       

ord4 <- capscale(as.dist(bac_bnti1)~vars$TotTreat+vars$PercAeFnd+vars$PercCxFnd+vars$AvgAnDens+vars$AvgCxDens+Condition(vars$PCNM2+vars$PCNM3+vars$PCNM4+vars$PCNM11+vars$PCNM17+vars$PCNM30+vars$PCNM33))
anova(ord4) #p = 0.001
anova.cca(ord4, by="term") #TotTreat (p = 0.003), PercAeFnd (p = 0.080), PercCxFnd (p = 0.521), AvgAnDens (p = 0.002), AvgCxDens (p = 0.006) 

ord5 <- capscale(as.dist(bac_bnti1)~vars$PCNM2+vars$PCNM3+vars$PCNM4+vars$PCNM11+vars$PCNM17+vars$PCNM30+vars$PCNM33+Condition(vars$TotTreat+vars$PercAeFnd+vars$PercCxFnd+vars$AvgAnDens+vars$AvgCxDens))
anova(ord5) #p = 0.001
anova.cca(ord5, by="term") #PCNM2 (p = 0.212), PCNM3 (p = 0.067), PCNM4 (p = 0.027), PCNM11 (p = 0.046), PCNM17 (p = 0.186), PCNM30 (p = 0.147), PCNM33 (p = 0.029)
