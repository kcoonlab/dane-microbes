set.seed(123)

library(PMCMRplus)
library(ggplot2)
library(ggpubr)

## Fig S1A

data=read.csv("Hist-Master-Final.csv",header=TRUE)
data=subset(data,select=-c(DipID,SampDate,Year))
sum.by.site=aggregate(.~AreaNumb,data,function(x) sum(x,na.rm=TRUE),na.action=na.pass)
sum.by.site[is.na(sum.by.site)]=0  
sum.by.site$PercSiteDry=sum.by.site$SiteDry/sum.by.site$SiteVis
sum.by.site=sum.by.site[order(sum.by.site$AreaNumb),]
PercSiteDry=100*(sum.by.site$PercSiteDry)

data=read.csv("Hist-Master-Final.csv",header=TRUE)
data=data[which(data$SiteDry==0),]
data=subset(data,select=-c(DipID,SampDate,Year))
sum.by.site=aggregate(.~AreaNumb,data,function(x) sum(x,na.rm=TRUE),na.action=na.pass)
sum.by.site[is.na(sum.by.site)]=0  
sum.by.site$PercTreat=sum.by.site$Treat/sum.by.site$SiteVis
sum.by.site=sum.by.site[order(sum.by.site$AreaNumb),]
TotTreat=sum.by.site$Treat
PercTreat=100*(sum.by.site$PercTreat)

data=read.csv("Hist-Master-Final.csv",header=TRUE)
data=data[which(data$Posttreat==0&data$SiteDry==0),]
data$AeFnd=ifelse(data$AeCnt>0,1,0)
data$AnFnd=ifelse(data$AnCnt>0,1,0)
data$CxFnd=ifelse(data$CxCnt>0,1,0)
data$NonVecFnd=ifelse((data$CxCnt-data$VecCnt)>0,1,0)
data$VecFnd=ifelse(data$VecCnt>0,1,0)
data$LarvDens=data$TotCnt/data$NumDips
data$AeDens=data$AeCnt/data$NumDips
data$AnDens=data$AnCnt/data$NumDips
data$CxDens=data$CxCnt/data$NumDips
data$NonVecDens=(data$CxCnt-data$VecCnt)/data$NumDips
data$VecDens=data$VecCnt/data$NumDips
data[is.na(data)]=0
data=subset(data,select=-c(DipID,SampDate,Year))
sum.by.site=aggregate(.~AreaNumb,data,function(x) sum(x,na.rm=TRUE),na.action=na.pass)
sum.by.site[is.na(sum.by.site)]=0  
sum.by.site$PercLarvFnd=sum.by.site$LarvFnd/sum.by.site$SiteVis
sum.by.site$PercAeFnd=sum.by.site$AeFnd/sum.by.site$SiteVis
sum.by.site$PercAnFnd=sum.by.site$AnFnd/sum.by.site$SiteVis
sum.by.site$PercCxFnd=sum.by.site$CxFnd/sum.by.site$SiteVis
sum.by.site$PercNonVecFnd=sum.by.site$NonVecFnd/sum.by.site$SiteVis
sum.by.site$PercVecFnd=sum.by.site$VecFnd/sum.by.site$SiteVis
sum.by.site=sum.by.site[order(sum.by.site$AreaNumb),]
avg.by.site=aggregate(.~AreaNumb,data,function(x) mean(x,na.rm=TRUE),na.action=na.pass)
avg.by.site=avg.by.site[order(avg.by.site$AreaNumb),]
PercLarvFnd=100*(sum.by.site$PercLarvFnd)
PercAeFnd=100*(sum.by.site$PercAeFnd)
PercAnFnd=100*(sum.by.site$PercAnFnd)
PercCxFnd=100*(sum.by.site$PercCxFnd)
PercNonVecFnd=100*(sum.by.site$PercNonVecFnd)
PercVecFnd=100*(sum.by.site$PercVecFnd)
AvgLarvDens=avg.by.site$LarvDens
AvgAeDens=avg.by.site$AeDens
AvgAnDens=avg.by.site$AnDens
AvgCxDens=avg.by.site$CxDens
AvgNonVecDens=avg.by.site$NonVecDens
AvgVecDens=avg.by.site$VecDens
AreaNumb=avg.by.site$AreaNumb

df<-data.frame(AreaNumb,PercSiteDry,TotTreat,PercTreat,PercLarvFnd,PercAeFnd,PercAnFnd,PercCxFnd,PercNonVecFnd,PercVecFnd,AvgLarvDens,AvgAeDens,AvgAnDens,AvgCxDens,AvgNonVecDens,AvgVecDens)
clust.data=read.csv("alpha-div.csv",header=TRUE)
data.total=merge(clust.data,df,by="AreaNumb")
data.total=subset(data.total,select=-c(faith_pd,observed_features,inv_simpson,cfus_per_ml))
#write.csv(data.total,"hist-data.csv",row.names = FALSE)

#prod <- read.csv("hist-data.csv")
hist(data.total$PercLarvFnd)

## Fig S1B

ggscatter(data.total,x="PercLarvFnd",y="AvgLarvDens",add="reg.line",conf.int=FALSE,cor.coef=FALSE,cor.method="spearman")
cor.test(data.total$PercLarvFnd, data.total$AvgLarvDens, method="spearman")

## Fig S1C
                      
data=read.csv("Hist-Master-Final.csv",header=TRUE)
data=data[which(data$Posttreat==0&data$LarvFnd==1),]
hist(data$LarvDiv)

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


