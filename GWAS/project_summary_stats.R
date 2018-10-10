library(devtools)
load_all('~/git/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/jia_projections/summary/jia.RDS'


shrink.DT <- readRDS(SHRINKAGE_FILE)
pc.emp <- readRDS(BASIS_FILE)
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
tmp <- basis.DT[trait=='control',] %>% t
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))
## project jia
proj.traits <- fread(TRAIT_MANIFEST)[grep("^jia",trait),]
proj.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE,proj.traits$trait)
proj.mat.emp<-create_ds_matrix(proj.DT,shrink.DT,SHRINKAGE_METHOD)
pred.emp <- predict(pc.emp,newdata=proj.mat.emp)
pred.DT <- data.table(trait=rownames(pred.emp),pred.emp)
#saveRDS(pred.DT,file=OUT_DIR)
## next compute projection 95% confidence intervals
var.DT <- readRDS(VARIANCE_FILE)
pred.DT <- melt(pred.DT,id.vars=c('trait'))
pred.DT <- merge(pred.DT,proj.traits[,.(trait,n1=cases,n=cases+controls)],by.x='trait',by.y='trait')
pred.DT <- merge(pred.DT,var.DT,by.x='variable',by.y='pc')
pred.DT[,ci.95:=sqrt((n/(n1 * (n-n1))) * mfactor) * 1.96  ]
pred.DT[,c('ci.lo','ci.hi'):=list(value-ci.95,value+ci.95)]
pred.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))
pred.DT <- merge(pred.DT.m,ctrl.DT,by='variable')
pred.DT[,variance:=(n/(n1 * (n-n1))) * mfactor]
pred.DT[,Z:=(value-control.loading)/sqrt(variance)]


pd <- position_dodge(0.3)
pp <- ggplot(pred.DT,aes(x=variable,y=value,group=trait,col=trait)) + geom_point(position=pd) +
geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi), width=.1, position=pd) + geom_line(position=pd)




## do hclust


## previously we had also plotted two bb traits also psoriasis and psoratic arthritis
## which I still need to do

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

p.DT <- get_phenotype_annotation(filter='20002_1477|20002_1453|20002_1464')
bb.files <- list.files(path='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',pattern="*.RDS",full.names=TRUE)
bb.DT<-lapply(bb.files,function(bb){
  readRDS(bb)
}) %>% rbindlist
bb.DT.m <- melt(bb.DT,id.vars='trait')
bb.DT.m <- merge(bb.DT.m,p.DT[,.(trait=phe,n1=n1,n=n1+n0)],by.x='trait',by.y='trait')
bb.DT.m <- merge(bb.DT.m,var.DT,by.x='variable',by.y='pc')
bb.DT.m[,ci.95:=sqrt((n/(n1 * (n-n1))) * mfactor) * 1.96  ]
bb.DT.m [,c('ci.lo','ci.hi'):=list(value-ci.95,value+ci.95)]
bb.DT.m [,variable:=factor(variable,levels=paste0('PC',1:11))]
bb.DT.m <- merge(bb.DT.m,ctrl.DT,by='variable')
bb.DT.m [,variance:=(n/(n1 * (n-n1))) * mfactor]
bb.DT.m [,Z:=(value-control.loading)/sqrt(variance)]



z.DT <- melt(rbind(pred.DT,bb.DT.m),id.vars=c('variable','trait'),measure.vars='Z')
mat.DT <- dcast(z.DT ,trait~variable)
mat <- as.matrix(mat.DT[,paste('PC',1:11,sep=''),with=FALSE])
mat <- rbind(mat,rep(0,ncol(mat)))
rownames(mat) <- c(mat.DT$trait,'control')
dist(mat) %>% hclust %>% plot
