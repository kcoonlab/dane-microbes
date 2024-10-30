library(picante)
library(iCAMP)
library(biomformat)
library(ape)
library(funrar)
library(caret)
library(phytools)
library(car)
library(ggplot2)

## Fig 2A

#biom convert -i exported/feature-table.biom -o exported/feature-table-json.biom --to-json 

Bac_otus<-"exported/feature-table.biom"
x1<-read_biom(Bac_otus)
meta<-read.csv("metadata3.csv")
bac<-as(biom_data(x1),"matrix")
bac<-t(bac)
tree<-read.tree("exported-tree/tree.tree")
phydist<-cophenetic(tree)

#Testing for phylogenetic signal

prod<-read.csv("hist-data.csv")
prod_vp<-prod[,c(3:17)]
prod.impute<-preProcess(prod_vp,method="bagImpute")
prod_vp2<-predict(prod.impute,prod_vp)
bac2<-make_relative(bac)
prod2<-prod_vp2
bac3<-cbind(bac2,prod2)

PercSiteDry<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$PercSiteDry),check.names=FALSE))
colnames(PercSiteDry)[1]<-"PercSiteDry"
TotTreat<-t(data.frame(lapply(bac3[ ,1:962],weighted.mean,x=bac3$TotTreat),check.names=FALSE))
colnames(TotTreat)[1]<-"TotTreat"
PercTreat<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$PercTreat),check.names=FALSE))
colnames(PercTreat)[1]<-"PercTreat"
PercLarvFnd<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$PercLarvFnd),check.names=FALSE))
colnames(PercLarvFnd)[1]<-"PercLarvFnd"
PercAeFnd<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$PercAeFnd),check.names=FALSE))
colnames(PercAeFnd)[1]<-"PercAeFnd"
PercAnFnd<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$PercAnFnd),check.names=FALSE))
colnames(PercAnFnd)[1]<-"PercAnFnd"
PercCxFnd<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$PercCxFnd),check.names=FALSE))
colnames(PercCxFnd)[1]<-"PercCxFnd"
PercNonVecFnd<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$PercNonVecFnd),check.names=FALSE))
colnames(PercNonVecFnd)[1]<-"PercNonVecFnd"
PercVecFnd<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$PercVecFnd),check.names=FALSE))
colnames(PercVecFnd)[1]<-"PercVecFnd"
AvgLarvDens<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgLarvDens),check.names=FALSE))
colnames(AvgLarvDens)[1]<-"AvgLarvDens"
AvgAeDens<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgAeDens),check.names=FALSE))
colnames(AvgAeDens)[1]<-"AvgAeDens"
AvgAnDens<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgAnDens),check.names=FALSE))
colnames(AvgAnDens)[1]<-"AvgAnDens"
AvgCxDens<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgCxDens),check.names=FALSE))
colnames(AvgCxDens)[1]<-"AvgCxDens"
AvgNonVecDens<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgNonVecDens),check.names=FALSE))
colnames(AvgNonVecDens)[1]<-"AvgNonVecDens"
AvgVecDens<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgVecDens),check.names=FALSE))
colnames(AvgVecDens)[1]<-"AvgVecDens"
prod_otus<-cbind(TotTreat,PercTreat,PercLarvFnd,PercAeFnd,PercAnFnd,PercCxFnd,PercNonVecFnd,PercVecFnd,AvgLarvDens,AvgAeDens,AvgAnDens,AvgCxDens,AvgNonVecDens,AvgVecDens)

trait=prod_otus[,1]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,2]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,3]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,4]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,5]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,6]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,7]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,8]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,9]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,10]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,11]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,12]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,13] 
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=prod_otus[,14]
names(trait)=rownames(prod_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!

env<-read.csv("hist-env-data.csv")
env_vp<-env[,c(5:10)]

env.impute<-preProcess(env_vp,method="bagImpute")
env_vp2<-predict(env.impute,env_vp)
bac2<-make_relative(bac)
env2<-env_vp2
bac3<-cbind(bac2,env2)
AvgLSTDay<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgLSTDay),check.names=FALSE))
colnames(AvgLSTDay)[1]<-"AvgLSTDay"
AvgLSTNight<-t(data.frame(lapply(bac3[ ,1:962],weighted.mean,x=bac3$AvgLSTNight),check.names=FALSE))
colnames(AvgLSTNight)[1]<-"AvgLSTNight"
AvgEVI<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgEVI),check.names=FALSE))
colnames(AvgEVI)[1]<-"AvgEVI"
AvgNDVI<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgNDVI),check.names=FALSE))
colnames(AvgNDVI)[1]<-"AvgNDVI"
AvgGpp<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgGpp),check.names=FALSE))
colnames(AvgGpp)[1]<-"AvgGpp"
AvgNpp<-t(data.frame(lapply(bac3[,1:962],weighted.mean,x=bac3$AvgNpp),check.names=FALSE))
colnames(AvgNpp)[1]<-"AvgNpp"
env_otus<-cbind(PercSiteDry,AvgLSTDay,AvgLSTNight,AvgEVI,AvgNDVI,AvgGpp,AvgNpp)

trait=env_otus[,1]
names(trait)=rownames(env_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=env_otus[,2]
names(trait)=rownames(env_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=env_otus[,3]
names(trait)=rownames(env_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=env_otus[,4]
names(trait)=rownames(env_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=env_otus[,5]
names(trait)=rownames(env_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=env_otus[,6]
names(trait)=rownames(env_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!
trait=env_otus[,7]
names(trait)=rownames(env_otus)
phylosig(tree,trait,method="lambda",test=TRUE,nsim=999) #Significant (corrected p < 0.0001)!

#Community Assembly

#All bacteria

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

#Biotype 1

tree<-read.tree("exported-tree/tree.tree")
meta<-read.csv("metadata3.csv")
bac2<-cbind(bac,meta)
bac_one<-subset(bac2,Hclust==1)
bac_one2<-data.frame(bac_one[,c(1:962)],check.names=FALSE)
bac_one3<-bac_one2[,colSums(bac_one2!=0)>0]
bac_one4<-as.matrix(bac_one3)
onetree<-prune.sample(bac_one4,tree)
distone<-cophenetic(onetree)
assem_one<-qpen(comm = bac_one4, pd = distone, pd.big.wd = NULL,
                  pd.big.spname = NULL, tree = NULL,
                  bNTI = NULL, RC = NULL, ab.weight = TRUE,
                  meta.ab = NULL, exclude.conspecifics = FALSE,
                  rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
                  nworker = 4, memory.G = 50, project = NA,
                  wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)
one_assem_ratio<-data.frame(assem_one[["ratio"]])
one_assem_result<-data.frame(assem_one[["result"]])
one_assem_ratio['Biotype']='One'
one_assem_result['Biotype']='One'

#Biotype 2

bac_two<-subset(bac2,Hclust==2)
bac_two2<-data.frame(bac_two[,c(1:962)],check.names=FALSE)
bac_two3<-bac_two2[,colSums(bac_two2!=0)>0]
bac_two4<-as.matrix(bac_two3)
twotree<-prune.sample(bac_two4, tree)
disttwo<-cophenetic(twotree)
assem_two<-qpen(comm = bac_two4, pd = disttwo, pd.big.wd = NULL,
                  pd.big.spname = NULL, tree = NULL,
                  bNTI = NULL, RC = NULL, ab.weight = TRUE,
                  meta.ab = NULL, exclude.conspecifics = FALSE,
                  rand.time = 1000, sig.bNTI = 1.96, sig.rc = 0.95,
                  nworker = 4, memory.G = 50, project = NA,
                  wd = getwd(), output.detail = FALSE, save.bNTIRC = FALSE)
two_assem_ratio<-data.frame(assem_two[["ratio"]])
two_assem_result<-data.frame(assem_two[["result"]])
two_assem_ratio['Biotype']= 'Two'
two_assem_result['Biotype']= 'Two'

#Violin plot

bnti <- rbind(one_assem_result,two_assem_result)

kruskal.test(bNTI~Biotype,data=bnti) #Significant (p = 0.0002001)!
bnti2<-bnti[!(bnti$bNTI<6),]
bnti3<-bnti[!(bnti$bNTI>6),]

viol <-ggplot(bnti3, aes(x=Biotype, y=bNTI)) + 
  geom_violin(trim=TRUE,fill="gray")+
  labs(x="Taxon", y ="ß-NTI")+
  geom_boxplot(width=0.1,outlier.shape=NA)+
  coord_cartesian(ylim = c(-6,6)) +
  theme_classic()+
  geom_hline(yintercept=2, linetype="dashed", color = "black",size=1.5) +
  geom_hline(yintercept=-2, linetype="dashed", color = "black",size=1.5) +
  theme(axis.title=element_text(size=20)) +
  annotate("text", x=1.85,y=-4.5,label="Biotype: ~italic(P) < 0.001", parse=TRUE, size=5, fontface="bold")+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin = 2, ymax = Inf, alpha = .3, fill = "#0571B0") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin = -Inf, ymax = -2, alpha = .3, fill = "#0571B0") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin = -2, ymax = 2, alpha = .3, fill = "#CA0020") +
  theme(text=element_text(size=20)) 
viol

## Fig 2B

#z tests for assembly processes

#Selection overall

n_s_r=sum((one_assem_result$process == "Homogeneous.Selection")+(one_assem_result$process == "Heterogeneous.Selection"))
t_s_r=nrow(one_assem_result)
n_s_d=sum((two_assem_result$process == "Homogeneous.Selection")+(two_assem_result$process == "Heterogeneous.Selection"))
t_s_d=nrow(two_assem_result)

prop.test(c(n_s_r,n_s_d),c(t_s_r,t_s_d)) #Significant (p < 0.0001; Biotype 1 overall higher than Biotype 2)!

#Hetero selection

n_hs_r=sum(one_assem_result$process == "Heterogeneous.Selection")
t_hs_r=nrow(one_assem_result)
n_hs_d=sum(two_assem_result$process == "Heterogeneous.Selection")
t_hs_d=nrow(two_assem_result)

prop.test(c(n_hs_r,n_hs_d),c(t_hs_r,t_hs_d)) #Significant (p = 0.01302; Biotype 1 higher than Biotype 2)!

#Homo selection

n_hs_r=sum(one_assem_result$process == "Homogeneous.Selection")
t_hs_r=nrow(one_assem_result)
n_hs_d=sum(two_assem_result$process == "Homogeneous.Selection")
t_hs_d=nrow(two_assem_result)

prop.test(c(n_hs_r,n_hs_d),c(t_hs_r,t_hs_d)) #Significant (p = 0.0001444; Biotype 1 overall higher than Biotype 2)!

#Homo dispersal

n_hd_r=sum(one_assem_result$process == "Homogenizing.Dispersal")
t_hd_r=nrow(one_assem_result)
n_hd_d=sum(two_assem_result$process == "Homogenizing.Dispersal")
t_hd_d=nrow(two_assem_result)

library(stats)
A = matrix(c(n_hd_r,n_hs_d,(t_hs_r-n_hs_r), (t_hs_d-n_hs_d)), nrow = 2)
fisher.test(A) #Significant (p = 0.008049)!

#Dispersal Limitation

n_dl_r=sum(one_assem_result$process == "Dispersal.Limitation")
t_dl_r=nrow(one_assem_result)
n_dl_d=sum(two_assem_result$process == "Dispersal.Limitation")
t_dl_d=nrow(two_assem_result)

prop.test(c(n_dl_r,n_dl_d),c(t_dl_r,t_dl_d)) #Not significant (p = 0.4386)

#Drift

n_d_r=sum(one_assem_result$process == "Undominated")
t_d_r=nrow(one_assem_result)
n_d_d=sum(two_assem_result$process == "Undominated")
t_d_d=nrow(two_assem_result)

prop.test(c(n_d_r,n_d_d),c(t_d_r,t_d_d)) #Significant (p = 0.001381; Biotype 2 overall higher than Biotype 1)!

#Bar plots

library(reshape2)
library(ggplot2)

bac_assem<-rbind(one_assem_ratio,two_assem_ratio)
bac_assem2<-melt(bac_assem,id.vars=c("Biotype"))
bac_assem2<-bac_assem2[c(1:10),]
bac_assem2$variable<-as.character(bac_assem2$variable)
bac_assem2[bac_assem2=="Undominated"]<-"Drift"
bac_assem2[bac_assem2=="Dispersal.Limitation"]<-"Dispersal Limitation"
bac_assem2[bac_assem2=="Homogenizing.Dispersal"]<-"Homogenizing Dispersal"
bac_assem2[bac_assem2=="Homogeneous.Selection"]<-"Homogeneous Selection"
bac_assem2[bac_assem2=="Heterogeneous.Selection"]<-"Variable Selection"
bac_assem2$variable<-factor(bac_assem2$variable,levels = c("Homogenizing Dispersal","Dispersal Limitation","Drift","Variable Selection","Homogeneous Selection"))
bac_assem2$variable

bac_bar <-ggplot(bac_assem2,aes(fill=variable,y=value,x=Biotype)) + 
  geom_bar(position="fill",stat="identity") +
  scale_y_continuous("Proportion") +
  scale_fill_manual(values=c("#CA0020","#F4A582","#FDDBC7","#92C5DE","#0571B0")) +
  theme_classic() +
  #theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill="Assembly Process") +
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=20)) +
  annotate("text",x=2,y=.9855,label='*',size=8,fontface="bold")+
  annotate("text",x=2,y=.24,label='*',size=8,fontface="bold")+
  annotate("text",x=1,y=.17,label='**',size=8,fontface="bold")+
  annotate("text",x=1,y=.06,label='***',size=8,fontface="bold")+
  theme(text=element_text(size=20)) 
bac_bar






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















