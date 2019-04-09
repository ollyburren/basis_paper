## collate PSA basis results - have only done UKBB SRD and JIA

## load in UKBB results

BB_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/psa_ss_shrink_2018'
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
## compose an annotation bar

bb_cod.sr_disease <- fread("/home/ob219/share/as_basis/GWAS/bb_projections/coding6.tsv")
bb_cod.sr_cancer <- fread("/home/ob219/share/as_basis/GWAS/bb_projections/coding3.tsv")
bb_cod.sr_medication <- fread("/home/ob219/share/as_basis/GWAS/bb_projections/coding4.tsv")

ukbb.traits <- ukbb$trait %>% unique %>% gsub("^bb\\_","",.)
anno.bb <- anno[phe %in% ukbb.traits]
anno.bb[,c('clade','coding'):=tstrsplit(code,'_') %>% lapply(.,as.numeric)]
anno.bb[clade=='20002',top.node:='bb_disease']
anno.bb[clade=='20001',top.node:='bb_cancer']
anno.bb[clade=='20003',top.node:='bb_medications']
ukbb <- merge(ukbb,anno.bb[,.(phe,n1,n0,category=top.node)],by.x='trait',by.y='phe')
ukbb <- ukbb[,trait:=paste('bb',trait,sep='_')]

## load in jia sub types

jia <- readRDS("/home/ob219/share/as_basis/GWAS/jia_projections/summary/psa_jia_2019.RDS")
sample.DT <- fread('/home/ob219/share/Data/GWAS/jia-mar-2019/summary-stats-samplecount-mar2019.csv')
setnames(sample.DT,'n','n1')
sub.DT <- data.table(idx=0:9,subtype=c('case','sys','PO','EO','RFneg','RFpos','ERA','PsA','undiff','missing'))
samp.DT <- merge(sub.DT,sample.DT[,.(ilar_pheno,n1)],by.x='idx',by.y='ilar_pheno')
samp.DT[,n0:=9196]
jia <- melt(jia,id.var='trait')
samp.DT[,subtype:=sprintf("jia_%s_19",subtype)]
jia <- merge(jia,samp.DT[,.(subtype,n1,n0)],by.x='trait',by.y='subtype')
jia[,category:='bowes_jia_2019']


all.proj <- list(
  ukbb = ukbb,
  jia = jia
) %>% rbindlist

all.proj[,n:=n1+n0]

## load in basis and variance
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/psa_ss_av_june.RDS'
var.DT <- readRDS(VARIANCE_FILE)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/psa_ss_basis_gwas.RDS'
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
all.DT[,delta:=value-control.loading]

saveRDS(all.DT,'/home/ob219/share/as_basis/GWAS/RESULTS/08_04_19_PSA_summary_results.RDS')

## FOREST PLOT CODE

library(cowplot)

res.DT <- all.DT
all.traits <- traits<-split(res.DT$trait,res.DT$category) %>% lapply(.,unique)
## add in basis traits for comparison read in again just in case changed
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
basis.DT[,category:='basis']

## another idea is to collapse biobank and basis traits where applicable

BB_LU <- list(
  CD = 'crohns.disease',
  CEL = 'malabsorption.coeliac.disease',
  MS = 'multiple.sclerosis',
  RA = 'rheumatoid.arthritis',
  SLE = 'systemic.lupus.erythematosis.sle',
  T1D = 'type.1.diabetes',
  UC = 'ulcerative.colitis',
  PBC = 'PBC',
  PSC = 'PSC',
  asthma = 'asthma',
  PSA = 'PSA'
)

talk.DT <- res.DT
talk.DT<-talk.DT[category %in% talk.DT[p.adj<0.05,]$category,]
at <- talk.DT$category %>% unique
at <- at[!at %in% c('bowes_jia_2019')]
library(RColorBrewer)
cols<-brewer.pal(length(talk.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c('bowes_jia_2019',at,'basis')
cols['bowes_jia_2019'] <- 'deeppink2'
cols['basis'] <- '#7c0799'
cols['bb_disease'] <- '#05af6e'
talk.DT[,trait:=gsub("^bb_","",trait)]
talk.DT[,trait:=gsub("^jia_","",trait)]

## perhaps the other way around

forest_plot_focal_merge <- function(proj.dat,basis.dat=basis.DT,pc,focal,title,cat_levels,fdr_thresh=0.05,theme=NA){
  dat <- proj.dat[variable==pc & (p.adj<fdr_thresh | trait %in% focal | trait %in% unlist(BB_LU)),]
  #dat[trait %in% focal,category:='aa_focal_diseases']
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  ## define a thing called trait label
  dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
  #for(i in seq_along(BB_LU)){
  #  tra <- BB_LU[[i]]
  #  dat[trait==tra,trait:=names(BB_LU)[i]]
  #}
  for(i in seq_along(BB_LU)){
    tra <- names(BB_LU)[i]
    dat[trait==tra,trait:=BB_LU[[i]]]
  }
  dat[,category:=factor(category,levels=names(cat_levels))]
  ## sort focal trait first
  foc.dt <- dat[trait %in% focal,.(trait,delta,category,n1)]
  #foc.dt[,tl:=sprintf("%s(%s)",trait,format(n1,big.mark=",",scientific=FALSE))]
  foc <- foc.dt[order(delta,decreasing=TRUE),.(trait,n1)]
  ## next do the rest
  nf.dt <- dat[!trait %in% focal,.(trait,delta,category,n1)][order(category,decreasing=FALSE),][!duplicated(trait),]
  #nf.dt[,tl:=sprintf("%s(%s)",trait,format(n1,big.mark=",",scientific=FALSE))]
  nf.dt[is.na(n1),tl:=trait]
  nfoc <- nf.dt[order(delta,decreasing=TRUE),.(trait,n1)]
  trait.order <- c(nfoc$trait,foc$trait)
  #dat[!duplicated(trait),][order(trait %in% focal,delta,decreasing=TRUE),]$trait
  #dat[,trait:=factor(trait,levels=dat[order(delta,!trait %in% focal,decreasing=TRUE),]$trait)]
  dat[,trait:=factor(trait,levels=trait.order)]
  #rbind(nfoc,foc)[,tl:=sprintf("%s(%s)",trait,format(n1,big.mark=",",scientific=FALSE))]
  ldat <- rbind(nfoc,foc)[,tl:=sprintf("%s",trait)]
  ldat[is.na(n1),tl:=trait]
  levels(dat$trait) <- ldat$tl
  #ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point(aes(size=log10(n1))) + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  #coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc)
  ## next alter the labels to include sample size
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }
  ggplot(dat,aes(x=trait,y=delta,colour=category,alpha=p.adj<fdr_thresh,lty=p.adj<fdr_thresh)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(title) + theme +
  xlab("Trait") + ylab("Change in basis score from control") + guides(alpha=FALSE,lty=FALSE) +
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.5)) + scale_colour_manual("Category",values=cat_levels,labels=c('JIA','UKBB SR Disease','Basis')) +
  scale_linetype_manual(values=c('TRUE'=1,'FALSE'=2)) +
  theme(axis.text.y=element_text(size=12))

}

pc <- 'PC3'
title <- sprintf("JIA subtypes %s (Variance Explained %0.1f%%) PSA Basis",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
pp1 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits['bowes_jia_2019'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols)
saveRDS(pp1,file="~/tmp/psa_basis_forest_plot_grob.RDS")
