set.seed(123)

library(picante)
library(iCAMP)
library(biomformat)
library(ape)
library(funrar)
library(caret)
library(phytools)
library(car)
library(ggplot2)
library(reshape2)

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

#Variable selection

n_hs_r=sum(one_assem_result$process == "Heterogeneous.Selection")
t_hs_r=nrow(one_assem_result)
n_hs_d=sum(two_assem_result$process == "Heterogeneous.Selection")
t_hs_d=nrow(two_assem_result)

prop.test(c(n_hs_r,n_hs_d),c(t_hs_r,t_hs_d)) #Significant (p = 0.01302; Biotype 1 higher than Biotype 2)!

#Homogeneous selection

n_hs_r=sum(one_assem_result$process == "Homogeneous.Selection")
t_hs_r=nrow(one_assem_result)
n_hs_d=sum(two_assem_result$process == "Homogeneous.Selection")
t_hs_d=nrow(two_assem_result)

prop.test(c(n_hs_r,n_hs_d),c(t_hs_r,t_hs_d)) #Significant (p = 0.0001444; Biotype 1 overall higher than Biotype 2)!

#Homogenizing dispersal

n_hd_r=sum(one_assem_result$process == "Homogenizing.Dispersal")
t_hd_r=nrow(one_assem_result)
n_hd_d=sum(two_assem_result$process == "Homogenizing.Dispersal")
t_hd_d=nrow(two_assem_result)

library(stats)
A = matrix(c(n_hd_r,n_hs_d,(t_hs_r-n_hs_r), (t_hs_d-n_hs_d)), nrow = 2)
fisher.test(A) #Significant (p = 0.008049)!

#Dispersal limitation

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
