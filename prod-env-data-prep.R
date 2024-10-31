set.seed(123)

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
#write.csv(data.total,"input-files/hist-data.csv",row.names = FALSE)
