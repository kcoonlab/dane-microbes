set.seed(123)

library(PMCMRplus)
library(ggplot2)
library(ggpubr)

## Fig S1A

prod <- read.csv("GitHub-files/hist-data.csv")
hist(prod$PercLarvFnd)

## Fig S1B

ggscatter(prod,x="PercLarvFnd",y="AvgLarvDens",add="reg.line",conf.int=FALSE,cor.coef=FALSE,cor.method="spearman")
cor.test(prod$PercLarvFnd, prod$AvgLarvDens, method="spearman")

## Fig S1C
                      
prod2 <- read.csv("Hist-Master-Final.csv",header=TRUE)
prod3 <- prod2[which(prod2$Posttreat==0&prod2$LarvFnd==1),]
prod4 <- prod3[which(prod3$LarvDiv>0),]
hist(prod4$LarvDiv)

## Fig S2

## Fig S3A

## Fig S3B

## Fig S3C

## Fig S3D

## Fig S4A

data=read.csv("alpha-taxa-div2.csv",header=TRUE)
data.all=data[which(data$Taxon=="all"),]
test.all=kruskalTest(Count~as.factor(Hclust),data=data.all) 
test.all #Significant (p = 0.01957)!
data.actino=data[which(data$Taxon=="actino"),]
test.actino=kruskalTest(Count~as.factor(Hclust),data=data.actino)
test.actino #Significant (p = 0.002542)!
data.cyano=data[which(data$Taxon=="cyano"),]
test.cyano=kruskalTest(Count~as.factor(Hclust),data=data.cyano)
test.cyano #Significant (p < 0.0001)!
data.firm=data[which(data$Taxon=="firm"),]
test.firm=kruskalTest(Count~as.factor(Hclust),data=data.firm)
test.firm #Significant (p < 0.0001)!
data.plancto=data[which(data$Taxon=="plancto"),]
test.plancto=kruskalTest(Count~as.factor(Hclust),data=data.plancto)
test.plancto #Significant (p = 0.03892)!
data.proteo=data[which(data$Taxon=="proteo"),]
test.proteo=kruskalTest(Count~as.factor(Hclust),data=data.proteo)
test.proteo #Significant (p = 0.0002541)!
data.verruco=data[which(data$Taxon=="verruco"),]
test.verruco=kruskalTest(Count~as.factor(Hclust),data=data.verruco)
test.verruco #Not significant (p = 0.5094)!
ggplot(data,aes(x=Taxon,y=Count,fill=as.factor(Hclust)))+geom_boxplot(outlier.shape=NA)

## Fig S4B

data=read.csv("alpha-div.csv",header=TRUE)
test.faith=kruskalTest(faith_pd~as.factor(Hclust),data=data) 
test.faith #Not significant (p = 0.1212)!
test.obs=kruskalTest(observed_features~as.factor(Hclust),data=data) 
test.obs #Significant (p = 0.03508)!
test.invsimp=kruskalTest(inv_simpson~as.factor(Hclust),data=data) 
test.invsimp #Not significant (p = 0.08357)!
test.dens=kruskalTest(log10(cfus_per_ml)~as.factor(Hclust),data=data) 
test.dens #Significant (p = 0.008604)!
ggplot(data,aes(x=Hclust,y=log10(cfus_per_ml),fill=as.factor(Hclust)))+geom_boxplot(outlier.shape=NA)

## Fig S5

## TotTreat by biotype

prod <- read.csv("GitHub-files/hist-data.csv")
test <- kruskalTest(TotTreat ~ as.factor(Hclust), data = prod)
test #Significant (p = 0.04583)
ggplot(prod,aes(x=Hclust,y=TotTreat,fill=as.factor(Hclust)))+geom_boxplot(outlier.shape=NA)

## PercAeFnd by biotype

test2 <- kruskalTest(PercAeFnd ~ as.factor(Hclust), data = prod)
test2 #Significant (p = 0.002423)
ggplot(prod,aes(x=Hclust,y=PercAeFnd,fill=as.factor(Hclust)))+geom_boxplot(outlier.shape=NA)

## AvgAnDens by biotype

test3 <- kruskalTest(AvgAnDens ~ as.factor(Hclust), data = prod)
test3 #Not significant (p = 0.06551)
ggplot(prod,aes(x=Hclust,y=AvgAnDens,fill=as.factor(Hclust)))+geom_boxplot(outlier.shape=NA)

## AvgCxDens by biotype

test4 <- kruskalTest(AvgCxDens ~ as.factor(Hclust), data = prod)
test4 #Not significant (p = 0.1417)
ggplot(prod,aes(x=Hclust,y=AvgCxDens,fill=as.factor(Hclust)))+geom_boxplot(outlier.shape=NA)
