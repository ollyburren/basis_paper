## overall file to collate all projection results

## load in UKBB results

BB_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018'
ukbb <- lapply(list.files(path=BB_DIR,pattern="*.RDS",full.names=TRUE),readRDS) %>% rbindlist
ukbb <- melt(ukbb,id.var='trait')
## need to add case and control numbers
get_phenotype_annotation <- function(filter='20002\\_'){
  bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
  pheno <- fread(bb_phenofile)
  setnames(pheno,names(pheno) %>% make.names)

  #med <- pheno[grepl("2000[23]\\_",Phenotype.Code) & Sex=='both_sexes',]
  med <- pheno[grepl(filter,Phenotype.Code) & Sex=='both_sexes',]
  ## load in phenotype file
  P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/phenotypes.both_sexes.tsv')
  P<-P[,.(phenotype,variable_type,non_missing=n_non_missing,cases=n_cases,controls=n_controls)]
  med[,phe:=make.names(Phenotype.Description) %>% gsub("(Non.cancer.illness.code..self.reported..)|(Treatment.medication.code..)|(Cancer.code..self.reported..)","",.)]
  med <- med[,.(Phenotype.Code,phe)]
  med <- merge(P,med,by.x='phenotype',by.y='Phenotype.Code')[,.(code=phenotype,phe,n1=as.numeric(cases),n0=as.numeric(controls))]
  med
}

anno <- get_phenotype_annotation(".*")[!is.na(n1) & !is.na(n0) & phe != 'unclassifiable',]
ukbb <- merge(ukbb,anno[,.(phe,n1,n0)],by.x='trait',by.y='phe')
ukbb <- ukbb[,trait:=paste('bb',trait,sep='_')]


## load in jia sub types

jia <- readRDS("/home/ob219/share/as_basis/GWAS/jia_projections/summary/jia.RDS")
## get sample counts and merge
library(annotSnpStats)
(load('~/share/Data/GWAS/JIA-2017-data/annotsnpstats-22.RData'))
samp.DT <- samples(G) %>% data.table
samp.DT <- samp.DT[,list(n1=.N),by=alt_ilar_code]
samp.DT <- samp.DT[,n0:=5181][!is.na(alt_ilar_code),]
samp.DT[,alt_ilar_code:=paste('jia',alt_ilar_code,sep='_')]
jia <- melt(jia,id.var='trait')
jia <- merge(jia,samp.DT,by.x='trait',by.y='alt_ilar_code')

## load in tian_infectious_disease

tian <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/tian_infectious_disease.RDS")
## get sample counts and merge
tian_samples <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/sample_size.RDS")
tian <- melt(tian,id.var='trait')
tian <- merge(tian,tian_samples,by='trait')

## load in ferreira_asthma

ferriera <- readRDS("/home/ob219/share/as_basis/GWAS/ferreira_projections/ferreira_asthma.RDS")
ferreira_samples <- readRDS("/home/ob219/share/as_basis/GWAS/ferreira_projections/sample_size.RDS")
ferriera <- melt(ferriera,id.var='trait')
ferriera <- merge(ferriera,ferreira_samples,by='trait')



all.proj <- list(ferriera=ferriera,tian=tian,jia=jia,ukbb=ukbb) %>% rbindlist
all.proj[,n:=n1+n0]

## load in basis and variance
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'
var.DT <- readRDS(VARIANCE_FILE)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,value:=(value-control.loading)]
basis.DT[,Z:=sign(value)]
basis.DT[,p.adj:=1]
basis.DT[,control.loading:=NULL]

## compute the variance of a projection
all.DT <- merge(all.proj,var.DT,by.x='variable',by.y='pc')
all.DT <- merge(all.DT,control.DT,by.x='variable',by.y='PC')
all.DT[,variance:=((log(n)-(log(n1) + log(n-n1)))+ log(mfactor)) %>% exp]
## add in control loading
all.DT[,Z:=(value-control.loading)/sqrt(variance)]
all.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
all.DT[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']

## how many traits have a pc score that is less than FDR 5%
keep.traits <- all.DT[p.adj<0.05,]$trait %>% unique
out.DT <- all.DT[trait %in% keep.traits,]
setnames(out.DT,'variable','PC')
out.DT[,delta:=value-control.loading]
## vitamin c is in twice remove !
out.DT <- out.DT[trait!='bb_vitamin.c.product',]
## work out hclust order
deltas <- melt(out.DT,id.vars=c('trait','PC'),measure.vars='delta') %>% dcast(.,trait~PC+variable)
delta.mat <- as.matrix(deltas[,-1])
rownames(delta.mat) <- deltas$trait
delta.hc <- dist(delta.mat) %>% hclust
out.DT[,trait:=factor(trait,levels=delta.hc$labels[delta.hc$order])]
out.DT[,PC:=factor(PC,levels=paste('PC',1:11,sep=""))]
out.DT[,Zplot:=]

library(cowplot)
ggplot(out.DT[p.adj<0.05 & !PC %in% c('PC10','PC11'),],aes(x=trait,y=PC,fill= (pmax(Z,-11) %>% pmin(.,11)),label=signif(delta,digits=1))) + geom_tile(color='black') +
geom_text(size=5,angle=90) + scale_fill_gradientn("Z",colours=c('green','white','red')) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
