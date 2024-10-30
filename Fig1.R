setwd("/Users/kerricoon/Library/CloudStorage/Box-Box/Projects/Dane-Redo-Master-Complete/Dane-Redo-Master")

library(ggplot2)
library(Rmisc)
library(cluster)
library(clusterSim)
library(clustertend)
library(factoextra)
#require(devtools)
#install_github("tapj/BiotypeR")
#library(BiotypeR)
#data(Titanium16S)
#Titanium16S.biotypes=biotyper.data.frame(Titanium16S, k=3, manalysis=TRUE)
library(BiotypeR)
library(ade4)

## Fig 1A

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
register_google(key = "AIzaSyBxYNOj8Cv7e25nm0k3MvXy2U3IxjBW18c")

states=map_data("state")
counties=map_data("county")
dane_county=counties[which(counties$region=="wisconsin" & counties$subregion=="dane"),]
us_base=ggplot(data=states,mapping=aes(x=long,y=lat,group=group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color="black",fill="gray")
us_base + theme_nothing()
us_base + theme_nothing() + 
  geom_polygon(data=dane_county,aes(x=long,y=lat,fill="orange"),color="white")

dane.data=read.csv("alpha-div.csv",header=TRUE)
dane.data=subset(dane.data,select=-c(Hclust,faith_pd,observed_features,inv_simpson,cfus_per_ml))
ll_means=sapply(dane.data[2:3],mean)
dane_map=get_map(location=ll_means,maptype="terrain-background",source="google",zoom=10)
ggmap(dane_map) + 
  geom_point(data=dane.data,color="red",size=4) +
  geom_text(data=dane.data,aes(label=paste("  ",as.character(AreaNumb),sep="")),angle=60,hjust=0,color="yellow")

## Fig 1B

data=read.csv("level2-taxa.csv",header=TRUE)
cbPalette=hcl.colors(8,palette="RdYlBu")
data$Taxon=factor(data$Taxon,levels=c("Actinobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Planctomycetes","Proteobacteria","Verrucomicrobia","Other"))
ggplot(data=data,aes(x=as.numeric(as.character(Order)),y=Rel,fill=factor(Taxon)))+geom_bar(stat="identity")+scale_fill_manual(values=cbPalette)

## Fig 1C

data=read.table("rel-level2-analysis-table.txt", header=T, row.names=1, dec=".", sep="\t")
data=data[-1,]
hopkins=hopkins(data,n=nrow(data)-1)
hopkins
JSD=function(x,y)sqrt(0.5*KLD(x,(x+y)/2)+0.5*KLD(y,(x+y)/2))
KLD=function(x,y)sum(x*log(x/y))
dist.JSD=function(inMatrix,pseudocount=0.000001,...){
  KLD=function(x,y)sum(x*log(x/y))
  JSD=function(x,y)sqrt(0.5*KLD(x,(x+y)/2)+0.5*KLD(y,(x+y)/2))
  matrixColSize=length(colnames(inMatrix))
  matrixRowSize=length(rownames(inMatrix))
  colnames=colnames(inMatrix)
  resultsMatrix=matrix(0,matrixColSize,matrixColSize)
  inMatrix=apply(inMatrix,1:2,function(x)ifelse(x==0,pseudocount,x))
    for(i in 1:matrixColSize){
      for(j in 1:matrixColSize){
        resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
		as.vector(inMatrix[,j]))
		}
	  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix,"method")="dist"
  return(resultsMatrix) 
  }
data.dist=dist.JSD(data)
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
                         require(cluster)
                         cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
                         return(cluster)
                        }
nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
    } else {
      data.cluster_temp=pam.clustering(data.dist, k)
      nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
      centrotypes = "medoids")
      }
	}
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
noise.removal <- function(dataframe, percent=1, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
  }
data.denoized=noise.removal(data, percent=1)

Dane16S.biotypes=biotyper.data.frame(data.denoized,k=2,manalysis=TRUE)
s.class(Dane16S.biotypes$PCA$li,fac=as.factor(Dane16S.biotypes$biotypes),grid=F)
Dane16S.biotypes$BET$tab

## Fig 1D

ggplot(data,aes(x=Taxon,y=Rel,fill=as.factor(Hclust)))+geom_boxplot(outlier.shape = NA)
