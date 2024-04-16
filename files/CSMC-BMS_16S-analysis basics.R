# install packages for later loading using the 'library()' function, only needs to be done once:
install.packages(c("vegan","phyloseq","ggplot2","dplyr","ape","reshape2")) # for analysis

# analysis:
library(vegan) # for alpha and beta diversity analyses
library(phyloseq) # for rarefying samples to equal numbers of reads
library(lmerTest)
library(lme4)
library(ggplot2)
library(dplyr)
library(ape)
library(reshape2)

outpath<-"~/Downloads/CSMC-BMS_16S_coding_lab/" # output directory for ASV tables and taxonomy files

# read in metadataa:
metadat<-read.csv(paste(outpath,"/metadata.csv",sep=""))
rownames(metadat)<-metadat$Run

# read in ASV table (generated from prior script):

tab<-read.csv(paste(outpath,"/dada2_output/seqtab_md5.csv",sep=""),row.names=1)


# check if all samples in ASV table are represented in row names of metadata:
tab[,!colnames(tab)%in%rownames(metadat)]

# create new metadata and ASV tables that include only samples that are in both files:
tab2<-tab[,colnames(tab)%in%rownames(metadat)]
metadat2<-metadat[rownames(metadat)%in%colnames(tab2),]


# gets rid of taxa that are below 0.001% abundance threshold:
tab2a<-tab2[rowSums(tab2)>sum(tab2)*0.00001,]


# examine read distribution of samples:
colSums(tab2a)[order(colSums(tab2a))]
# merge(metadat2,colSums(tab2a),by=0)

# rarefy reads:
ps <- phyloseq(otu_table(t(tab2a), taxa_are_rows = FALSE))
tab3rar2<-t(rarefy_even_depth(ps,sample.size=min(colSums(tab2a)),replace=F))
colSums(tab3rar2) #check if rarefaction has sub-sampled all to same depth


# calculate alpha diversity:
diversity(tab3rar2, MARGIN = 2, index = "invsimpson")
metadat3<-merge(metadat2,diversity(tab3rar2, MARGIN = 2, index = "invsimpson"),by=0)
colnames(metadat3)[ncol(metadat3)]<-"invsimp"
rownames(metadat3)<-metadat3[,1]

# graph alpha across subject groups:
metadat3$prepostFMT<-NA
metadat3$prepostFMT[metadat3$days_rel_to_FMT<=0]<-"pre"
metadat3$prepostFMT[metadat3$days_rel_to_FMT>0]<-"post"
metadat3$prepostFMT[metadat3$donor_recipient=="D"]<-"donor"
ggplot(metadat3,aes(x=prepostFMT,y=invsimp))+geom_boxplot(outlier.shape=NA,width=0.5)+geom_jitter(size=3,width=0.2, alpha=0.3)+theme_bw()
metadat3$prepostFMT<-factor(metadat3$prepostFMT,levels=unique(metadat3$prepostFMT)[c(3,2,1)])
ggplot(metadat3,aes(x=prepostFMT,y=invsimp))+geom_boxplot(outlier.shape=NA,width=0.5)+geom_jitter(size=3,width=0.2, alpha=0.3)+theme_bw()

# graph alpha diversity over time:
ggplot(metadat3,aes(x=sampling_rel_to_FMT,y=invsimp,color=as.factor(host_subject_id)))+geom_point(size=2)+theme_bw()+geom_line(aes(group=as.factor(host_subject_id)),size=1)

# scico cool palettes: batlow, roma, lipari, tokyo, hawaii
library(scico)
ggplot(metadat3,aes(x=sampling_rel_to_FMT,y=invsimp,color=as.factor(host_subject_id)))+geom_point(size=2)+theme_bw()+geom_line(aes(group=as.factor(host_subject_id)),size=1) + scale_color_manual(values = scico(length(unique(metadat3$host_subject_id)), palette = "roma"))


# remove donors:
ggplot(subset(metadat3,donor_recipient=="R"),aes(x=sampling_rel_to_FMT,y=invsimp))+geom_point(aes(color=as.factor(host_subject_id)),size=4)+theme_bw()+geom_line(size=3,aes(color=as.factor(host_subject_id)))+geom_vline(xintercept = 0,size=2)+geom_boxplot(aes(group=sampling_rel_to_FMT), alpha=0.3,size=0.7,fatten=5)+ guides(color=guide_legend(title="New Title"))+ scale_color_manual(values = scico(length(unique(metadat3$host_subject_id)), palette = "roma"))

# calculate beta diversity matrix and visualize via PCoA plot:
bray<-vegdist(t(tab3rar2),method="bray")
bray2<-as.matrix(bray)
pco1<-pcoa(bray2)
pco2<-pco1$vectors
pco3<-merge(pco2,metadat3,by=0)
ggplot(pco3,aes(x=Axis.1,y=Axis.2,color=prepostFMT))+geom_point(size=6,alpha=0.4)+theme_bw()

ggplot(pco3,aes(x=Axis.1,y=Axis.2,color=prepostFMT))+geom_point(size=6,alpha=0.4)+theme_bw()+stat_ellipse(size=0.75)

#test significance of beta diversity clustering:
bray3<-bray2[rownames(metadat3),rownames(metadat3)]
# line above is essential to ensure that order in the beta diversity matrix is the same as in the metadata file; adonis does not match samples in the beta diversity matrix to the metadat by name, but by order
formu<-formula(paste("bray3 ~ prepostFMT"))
adon<-adonis2(formu, data= metadat3, permutations=10000, by = "margin")
adon
# subset both metadata and beta div matrix to compare only pre vs post samples:
metasub<-subset(metadat3,prepostFMT!="donor")
bray3<-bray2[rownames(metasub),rownames(metasub)]
formu<-formula(paste("bray3 ~ prepostFMT"))
adon<-adonis2(formu, data= metasub, permutations=10000, by = "margin")
adon


# GO TO POWERPOINT:
# identifying taxa with differential abundance between donors and recipients:
# data structure examination for our 16S dataset:
hist(as.numeric(melt(tab3rar2)[,3]),breaks = 40)
# expected data structure for 'parametric' stats (e.g. t-test, Pearson correlation), so-called 'normal distribution':
hist(rnorm(1000))

# instead, can use non-parametric Mann-Whitney U tests (which do not assume normal distribution) to identify taxa in differential abundance between 2 groups:
wilxmat <-matrix(nrow=nrow(tab3rar2),ncol=5)
colnames(wilxmat)<-c("wilx_P","wilx_stat","log2fc","numnonzero","rank")
rownames(wilxmat)<-rownames(tab3rar2)
cases<-subset(metadat3,sampling_rel_to_FMT==(-1))$Run
controls<-subset(metadat3,donor_recipient=="D")$Run
for(j in 1:nrow(tab3rar2))
{
  wilxmat[j,1]<-wilcox.test(x= as.numeric(as.vector(tab3rar2[j,cases])),y= as.numeric(tab3rar2[j,controls]))$p.value
  wilxmat[j,2]<-wilcox.test(x= as.numeric(tab3rar2[j,cases]),y= as.numeric(tab3rar2[j,controls]))$statistic
  wilxmat[j,3]<-log(mean(as.numeric(tab3rar2[j,cases]))+1,base=10)-log(mean(as.numeric(tab3rar2[j,controls]))+1,base=10) #positive means higher abundance in cases
  # wilxmat[j,3]<-tab3rar2[j,cases] %>% as.numeric %>% mean %>% log(base=10) - tab3rar2[j,controls] %>% as.numeric %>% mean %>% log(base=10) #positive means higher abundance in cases
  wilxmat[j,4]<-sum(as.numeric(tab3rar2[j,])>	0)
  wilxmat[j,5]<-j
}
wilxmat2<-cbind(wilxmat,BH_Q=p.adjust(as.numeric(as.vector(wilxmat[,1])),method="BH"))
tax<-read.csv(paste(outpath,"/dada2_output/taxaspec_md5_fasta.csv",sep=""),row.names=1)
wilxmat3<-merge(wilxmat2,tax,by=0)
wilxmat5<-wilxmat3[order(wilxmat3$wilx_P,decreasing=F),]
write.table(wilxmat5,file=paste(outpath,"/mann-whitney_don-v-recipient.csv",sep=""),sep=",",col.names=NA)
wilxmat5[,!colnames(wilxmat5)%in%"fasta"]

ggplot(wilxmat5,aes(x=log2fc,y=-log(wilx_P,base=10)))+geom_point(aes(color=Class),alpha=0.6,size=4)+theme_bw()+geom_hline(yintercept = -log(0.05,base=10),linetype=2)

ggplot(subset(wilxmat5,wilx_P<0.05),aes(x=log2fc,y=-log(wilx_P,base=10)))+geom_point(aes(color=Class,size=-log(wilx_P,base=10)),alpha=0.6)+geom_point(data=subset(wilxmat5,wilx_P>0.05),aes(x=log2fc,y=-log(wilx_P,base=10)),color="gray")+theme_bw()+geom_hline(yintercept = -log(0.05,base=10),linetype=2)+scale_size_continuous(range = c(2, 6))




ttab3<-t(tab3rar2)
merg1<-merge(ttab3,metadat3,by=0)
merg2<-merg1
i<-"acef2cace17f5f8910b2279894d8aa3c"
colnames(merg2)[colnames(merg2)%in%i]<-"taxon"
ggplot(subset(merg2,sampling_rel_to_FMT==(-1)|donor_recipient=="D"),aes(x=donor_recipient,y=taxon))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0,width=0.2)+theme_bw()+ylab(wilxmat5[wilxmat5$Row.names%in%i,"famgenspec"])+scale_y_sqrt()




# sqrt transform with linear mixed effects models to control for confounding variables, to compare pre vs donor:
hist(as.numeric(melt(tab3rar2)[,"value"]))

sqrt3<-sqrt(tab3rar2)
hist(as.numeric(melt(sqrt3)[,"value"]))

hist(log(as.numeric(melt(tab3rar2)[,"value"])+1))

# asinTransform <- function(p) { asin(sqrt(p)) }
# tab3rar2rel<-apply(tab3rar2, 2, function(x) as.numeric(x)/sum(x))
# hist(asinTransform(as.numeric(melt(tab3rar2rel)[,"value"])))
# tab3rar2relasinsqrt<-apply(tab3rar2rel,2,function(x) as.numeric(asinTransform))

sqrt3b<-merge(t(sqrt3),metadat3,by=0)
sqrt3c<-sqrt3b[,2:ncol(sqrt3b)]
sqrt4<-subset(sqrt3c,sampling_rel_to_FMT==(-1)|donor_recipient=="D")

lmresmat<-matrix(ncol=2,nrow=ncol(sqrt4))
rownames(lmresmat)<-colnames(sqrt4)
colnames(lmresmat)<-c("P","t")
for(i in 1:ncol(sqrt4))
{
  if(colnames(sqrt4)[i]!="donor_recipient")
  {
    temp<-cbind(sqrt4[,c("donor_recipient","Host_Age","host_sex")],sqrt4[,i])
    temp<-as.data.frame(temp)
    colnames(temp)<-c("donor_recipient","Host_Age","host_sex","taxon")
    temp$taxon<-as.numeric(((as.numeric(as.vector(temp$taxon)))))
    if(dim(table(temp[,2]))>4)
    {
      
      lmresmat[i,1]<-coef(summary(lmer(formula("taxon~donor_recipient+Host_Age+(1|host_sex)"),data=temp)))[2,5]
      lmresmat[i,2]<-coef(summary(lmer(formula("taxon~donor_recipient+Host_Age+(1|host_sex)"),data=temp)))[2,4]
    }
  }
}
lmresmat2<-lmresmat[order(lmresmat[,1]),]
lmresmat3<-merge(lmresmat2,tax,by=0,all.y=F)
lmresmat4<-cbind(lmresmat3,BHq=p.adjust(lmresmat3$P,method="BH"))
lmresmat5<-lmresmat4[order(lmresmat4$P,decreasing=F),]
write.table(lmresmat5,file=paste(outpath,"/lme_DvR_contforconfound.csv",sep=""),sep=",",col.names=NA)
lmresmat5[,!colnames(lmresmat5)%in%c("fasta","famgenspec")]

ggplot(subset(lmresmat5,P<0.05),aes(x=t,y=-log(P,base=10)))+geom_point(aes(color=Class,size=-log(P,base=10)),alpha=0.6)+geom_point(data=subset(lmresmat5,P>0.05),aes(x=t,y=-log(P,base=10)),color="gray")+theme_bw()+geom_hline(yintercept = -log(0.05,base=10),linetype=2)+scale_size_continuous(range = c(2, 6))

# sqrt4[,c(which(colnames(sqrt4)%in%"9f7307d8375cf75a3ce495abf1e16b05"),which(colnames(sqrt4)%in%"treatment_group"))]

merg1<-merge(ttab3,metadat3,by=0)
merg2<-merg1
i<-"acef2cace17f5f8910b2279894d8aa3c"
colnames(merg2)[colnames(merg2)%in%i]<-"taxon"
ggplot(subset(merg2,sampling_rel_to_FMT==(-1)|donor_recipient=="D"),aes(x=donor_recipient,y=taxon+1))+geom_boxplot(outlier.shape=NA,width=0.64)+geom_jitter(height=0,width=0.2)+theme_bw()+ylab(lmresmat5[lmresmat5$Row.names%in%i,"famgenspec"])+scale_y_log10() #added a 'pseudocount' of 1 to facilitate log transformation


merg2<-sqrt3b
i<-"acef2cace17f5f8910b2279894d8aa3c"
colnames(merg2)[colnames(merg2)%in%i]<-"taxon"
merg2$prepostFMT<-factor(merg2$prepostFMT,levels=unique(merg2$prepostFMT)[c(2,1)])
ggplot(subset(merg2,donor_recipient!="D"),aes(x=prepostFMT,y=taxon))+geom_boxplot(outlier.shape=NA,width=0.64)+geom_jitter(height=0,width=0.2)+theme_bw()+ylab(lmresmat5[lmresmat5$Row.names%in%i,"famgenspec"])#+scale_y_sqrt()






# sqrt with linear mixed effects models to control for intra-individual covariance, identifies taxa that differ in abundance within individuals comparing pre to post-FMT:
sqrt3b<-merge(t(sqrt3),metadat3,by=0)
sqrt4<-subset(sqrt3b,donor_recipient!="D")

lmresmat<-matrix(ncol=2,nrow=ncol(sqrt4))
rownames(lmresmat)<-colnames(sqrt4)
colnames(lmresmat)<-c("P","t")
for(i in 1:ncol(sqrt4))
{
  if(colnames(sqrt4)[i]!="prepostFMT")
  {
    temp<-cbind(sqrt4[,'prepostFMT'],sqrt4[,i],sqrt4[,"host_subject_id"])
    temp<-as.data.frame(temp)
    colnames(temp)<-c("prepostFMT","taxon","host_subject_id")
    temp$taxon<-as.numeric(((as.numeric(as.vector(temp$taxon)))))
    if(dim(table(temp[,2]))>4)
    {
      
      lmresmat[i,1]<-coef(summary(lmer(formula("taxon~prepostFMT+(1|host_subject_id)"),data=temp)))[2,5]
      lmresmat[i,2]<-coef(summary(lmer(formula("taxon~prepostFMT+(1|host_subject_id)"),data=temp)))[2,4]
    }
  }
}
lmresmat2<-lmresmat[order(lmresmat[,1]),]
lmresmat3<-merge(lmresmat2,tax,by=0,all.y=F)
lmresmat4<-cbind(lmresmat3,BHq=p.adjust(lmresmat3$P,method="BH"))
lmresmat5<-lmresmat4[order(lmresmat4$P,decreasing=F),]
write.table(lmresmat5,file=paste(outpath,"/lme_prepostFMT.csv",sep=""),sep=",",col.names=NA)
lmresmat5[,!colnames(lmresmat5)%in%"fasta"]


sqrt4[,c(which(colnames(sqrt4)%in%"9f7307d8375cf75a3ce495abf1e16b05"),which(colnames(sqrt4)%in%"treatment_group"))]

merg1<-merge(ttab3,metadat3,by=0)
merg2<-merg1
i<-"884cb13b27b7dd2a3750d189c988f647"
colnames(merg2)[colnames(merg2)%in%i]<-"taxon"
ggplot(subset(merg2,sampling_rel_to_FMT==(-1)|donor_recipient=="D"),aes(x=donor_recipient,y=taxon))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0,width=0.2)+theme_bw()+ylab(lmresmat5[lmresmat5$Row.names%in%i,"famgenspec"])+scale_y_sqrt()


merg2<-sqrt3b
i<-"dd0234a1d48f74a011f58a58b206e6ad"
colnames(merg2)[colnames(merg2)%in%i]<-"taxon"
ggplot(subset(merg2,donor_recipient!="D"),aes(x=prepostFMT,y=taxon))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0,width=0.2)+theme_bw()+ylab(lmresmat5[lmresmat5$Row.names%in%i,"famgenspec"])+scale_y_sqrt()












