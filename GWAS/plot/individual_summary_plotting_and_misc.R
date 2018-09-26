library(data.table)
library(cowplot)
library(magrittr)
library(cupcake)


SHRINKAGE_METHOD<-'ws_emp'
DEFAULT.SNPSTATS.DIR <- '/home/ob219/rds/hpc-work/as_basis/snpStats/basis_1kg_feb/'

## load data
ref_af_file<-'/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab'
m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_july.tsv'
## dir where all preprocessed gwas files are.
## we expect variants to be reflected in ref_af_file, have there OR aligned and be defined in the manifest file
gwas_data_dir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/all_studies_filtered'
basis.DT<-get_gwas_data(m_file,ref_af_file,gwas_data_dir)
## there appears to be an error with SLE dataset in that 4 SNPs have p==0 which is impossible
## set these to 0.99 for the time being
 basis.DT[p.value==0,p.value:=0.99]
shrink.DT<-compute_shrinkage_metrics(basis.DT)

basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
DATA.DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/individual_proj'

## load in the data

fs <- list.files(path=DATA.DIR,pattern="*.RDS",full.names=TRUE)
all.res <- lapply(fs,function(f){
  cm <- readRDS(f) %>% colMeans %>% t %>% data.table
  cm[, trait:=basename(f) %>% gsub("_projection.RDS","",.)]
  melt(cm,'trait')
}) %>% rbindlist

## plot some of this

plot <- dcast(all.res[variable %in% c('PC1','PC3'),],trait~variable)
library(ggrepel)
ggplot(plot,aes(x=PC1,y=PC3,label=trait)) + geom_point() + geom_text_repel()


plot <- dcast(all.res[variable %in% c('PC1','PC2'),],trait~variable)
library(ggrepel)
ggplot(plot,aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text_repel()

m <- dcast(all.res,trait~variable)
as.matrix(m[,-1],rownames=m$trait) %>% dist %>% hclust %>% plot

## Bartletts test for heterozygisity of variances



all.ind.res <- lapply(fs,function(f){
  cm <- readRDS(f) %>% data.table
  cm[, trait:=basename(f) %>% gsub("_projection.RDS","",.)]
  melt(cm,'trait')
}) %>% rbindlist


all.ind.res[grep("^jia",trait),trait:='jia']
all.ind.res[,trait:=as.factor(trait)]

lapply(paste0('PC',1:11),function(pc){
  DT <- all.ind.res[variable==pc,]
  data.table(pc=pc,p.Bartlett=bartlett.test(value~trait,DT)$p.value,traits=paste(levels(DT$trait),sep=',',collapse=','))
}) %>% rbindlist

names(all.ind.res) <- basename(fs) %>% gsub("_projection.RDS","",.)

jia.idx <- grep("^jia",names(all.ind.res))

jia <- rbindlist(all.ind.res[jia.idx])
jia[,trait:='jia']
single.disease <- rbind(rbindlist(all.ind.res[-jia.idx]),jia)


traits <- c('pid','jdm')
pc <- 'PC3'
idx <- with(single.disease,which(variable==pc & trait %in% traits))
idx <- with(single.disease,which(variable==pc))
tra <- as.data.frame(single.disease[idx,.(trait=as.factor(trait),value)])
plot(value~trait,data=tra)
tra <- split(single.disease[idx,]$value,single.disease[idx,]$trait)


## look at JDM covariates
f <- fs[grep("jdm",fs)]
jdm.pc <- readRDS(f)
jdm.pc <- data.table(sampleid=rownames(jdm.pc),jdm.pc)

ao <- fread('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/deakin/AgeOnset.csv')
fam <- fread('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/jdm.fam')[,lu:=paste('X',V1,sep='')]
comb <- merge(fam,ao,by.x='V2',by.y='ID')
m<-merge(comb,jdm.pc,by.x='lu',by.y='sampleid')
lapply(paste0('PC',1:11),function(pc){
  DT <- m[,.(AgeOnset,pc=get(`pc`))]
  data.table(PC=pc,p.val=summary(lm(DT,formula=AgeOnset~pc))$coefficient['pc',4])
}) %>% rbindlist

p.DT <- data.table(pred.emp)
p.DT[,ind:=rownames(pred.emp)]
setkey(ao,ID)
setkey(p.DT,ind)
# someone has an age of onset less than zero not sure what this means remove it
p.DT<-ao[p.DT][AgeOnset>0]

## fit a simple linear model and return the p.value associated with the coefficent
lapply(paste0('PC',1:11),function(pc){
  DT <- p.DT[,.(AgeOnset,pc=get(`pc`))]
  data.table(PC=pc,p.val=summary(lm(DT,formula=AgeOnset~pc))$coefficient['pc',4])
}) %>% rbindlist




## what about autoantibody data
f <- fs[grep("jdm",fs)]
jdm.pc <- readRDS(f)
jdm.pc <- data.table(sampleid=rownames(jdm.pc),jdm.pc)
ab <- fread('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/jdm_autoantibody_300518.csv')
fam <- fread('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/jdm.fam')[,lu:=paste('X',V1,sep='')]
comb <- merge(fam,ab,by.x='V2',by.y='PatientID')
comb[,ab.positive:=ifelse(Autoantibody!='Negative',1,0)]
m<-merge(comb,jdm.pc,by.x='lu',by.y='sampleid')[!is.na(ab.positive),]
lapply(paste0('PC',1:11),function(pc){
  DT <- m[,.(ab.positive,pc=get(`pc`))]
  data.table(PC=pc,p.val=summary(glm(DT,formula=ab.positive~pc,family="binomial"))$coefficient['pc',4])
}) %>% rbindlist


## try a MANOVA approach

#combine rare categories where count<5 into 'other'

m[,Autoantibody:=make.names(Autoantibody)]
rare.abs<-m[,list(abc=.N),by=Autoantibody][abc<20,]$Autoantibody
m[,ab.category:=ifelse(Autoantibody %in% rare.abs,'Other',Autoantibody) %>% as.factor]

lapply(paste0('PC',1:11),function(pc){
  DT <- m[,.(ab.category,pc=get(`pc`))]
  data.table(PC=pc,p.val=summary(aov(pc~ab.category,data=DT))[[1]][1,5])
}) %>% rbindlist
aov.ex1 = aov(PC1~ab.category,data=m)

pl <- melt(m,id.vars=c('V2','ab.category'),measure.vars=paste0('PC',1:11))
ggplot(pl,aes(x=ab.category,y=value)) + geom_boxplot() + facet_wrap(~variable,scale="free") +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
