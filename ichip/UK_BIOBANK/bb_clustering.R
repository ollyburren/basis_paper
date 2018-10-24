## look at biobank on ichip

library(devtools)
load_all('~/git/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/ichip/support/shrinkage_ic.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/ichip/support/basis_ic.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/ichip/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/ichip/snp_manifest/ichip_all.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/ichip/trait_manifest/as_manifest_ichip.tsv'
VARIANCE_FILE <- '/home/ob219/share/as_basis/ichip/support/shrink_av_ichip.RDS'
#OUT_DIR <- '/home/ob219/share/as_basis/GWAS/jia_projections/summary/jia.RDS'

shrink.DT <- readRDS(SHRINKAGE_FILE)
pc.emp <- readRDS(BASIS_FILE)
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
tmp <- basis.DT[trait=='control',] %>% t
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))


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

var.DT <- readRDS(VARIANCE_FILE)
#p.DT <- get_phenotype_annotation(filter='20002_1477|20002_1453|20002_1464')
p.DT <- get_phenotype_annotation()
bb.files <- list.files(path='/home/ob219/share/as_basis/ichip/bb_projections/shrink_2018/',pattern="*.RDS",full.names=TRUE)
bb.DT<-lapply(bb.files,function(bb){
  readRDS(bb)
}) %>% rbindlist
bb.DT.m <- melt(bb.DT,id.vars='trait')
bb.DT.m <- merge(bb.DT.m,p.DT[,.(trait=phe,n1=n1,n=n1+n0)],by.x='trait',by.y='trait')
bb.DT.m <- merge(bb.DT.m,var.DT,by.x='variable',by.y='pc')
bb.DT.m[,ci.95:=sqrt((n/(n1 * (n-n1))) * mfactor) * 1.96  ]
bb.DT.m [,c('ci.lo','ci.hi'):=list(value-ci.95,value+ci.95)]
bb.DT.m [,variable:=factor(variable,levels=paste0('PC',1:13))]
bb.DT.m <- merge(bb.DT.m,ctrl.DT,by='variable')
bb.DT.m [,variance:=(n/(n1 * (n-n1))) * mfactor]
bb.DT.m [,Z:=(value-control.loading)/sqrt(variance)]
bb.DT.m[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
#bb.DT.m[,p.adj:=p.adjust(p.value),by='variable']
bb.DT.m[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']
bb.DT.m[,variable:=factor(variable,levels=paste0('PC',1:13))]
bb.DT.m[,short.trait:=substr(trait,1,15),]

traits.of.interest <- bb.DT.m[,list(sum(p.adj<0.01)),by='trait'][V1!=0,]$trait
bb.DT.m <- bb.DT.m[trait %in% traits.of.interest,]

z.DT <- melt(bb.DT.m,id.vars=c('variable','trait'),measure.vars='Z')

mat.DT <- dcast(z.DT ,trait~variable)
mat <- as.matrix(mat.DT[,paste('PC',1:11,sep=''),with=FALSE])
mat <- rbind(mat,rep(0,ncol(mat)))
rownames(mat) <- c(mat.DT$trait,'control')
dist(mat) %>% hclust %>% plot

## what about clustering of disease pairs how does this look ?


BB_LU <- list(
  ATD = 'hyperthyroidism.thyrotoxicosis',
  CEL = 'malabsorption.coeliac.disease',
  CRO = 'crohns.disease',
  MS = 'multiple.sclerosis',
  NAR = 'NAR',
  PBC = 'PBC',
  PSA = 'psoriatic.arthropathy',
  PSO = 'psoriasis',
  RA = 'rheumatoid.arthritis',
  SLE = 'systemic.lupus.erythematosis.sle',
  T1D = 'type.1.diabetes',
  UC = 'ulcerative.colitis',
  PSC = 'PSC',
  asthma = 'asthma'
)

all <- rbind(bb.DT[trait %in% do.call('c',BB_LU),],basis.DT)
all.mat <- all[,-1] %>% as.matrix
rownames(all.mat) <- all$trait
dist(all.mat) %>% hclust %>% plot
