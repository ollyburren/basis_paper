library(ggplot2)
library(cowplot)
library(ggrepel)
library(data.table)
library(magrittr)

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
  asthma = 'asthma'
)

plot_bb_hclust <- function(proj_dir,basis_file,variance_file,bb_lu,ptitle,hc.method="complete"){
  ## create a matrix of actual projections for each pc
  basis <- readRDS(basis_file)
  basis.DT <- data.table(trait=rownames(basis$x),basis$x)
  ## load in biobank projections
  bb.files <- list.files(path=proj_dir,pattern="*.RDS",full.names=TRUE)
  bb.DT<-lapply(bb.files,function(bb){
    readRDS(bb)
  }) %>% rbindlist
  bb.DT.filt <- bb.DT[trait %in% unlist(bb_lu),]
  bb.unl <- unlist(bb_lu)
  midx <- match(bb.DT.filt$trait,bb.unl)
  bb.DT.filt[,trait:=paste('bb',names(bb.unl)[midx],sep='_'),]
  #plot.DT <- rbindlist(list(bb.DT.filt,basis.DT,bb.DT[!trait %in% unlist(bb_lu),]))
  plot.DT <- rbindlist(list(bb.DT.filt,basis.DT))
  #bb.DT.back <- bb.DT[!trait %in% unlist(bb_lu),]

  proj <- plot.DT[trait!='control',paste0('PC',1:11),] %>% as.matrix
  #rownames(proj) <-
  proj<-rbind(proj,rep(0,ncol(proj)))
  rownames(proj) <- c(plot.DT[trait!='control',]$trait,'control')
  dist(proj) %>% hclust(.,method=hc.method) %>% plot(.,main=ptitle)
}


par(mfrow=c(2,2))

plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS.RDS',bb_lu=BB_LU,
ptitle = 'Beta 2018 UKBB Complete',hc.method="complete"
)

plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',bb_lu=BB_LU,
ptitle = 'Shrinkage 2018 UKBB Complete',hc.method="complete"
)

## add in PSA basis projections for comparison

#plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/psa_beta_2018/',
#basis_file='/home/ob219/share/as_basis/GWAS/support/psa_beta_basis_gwas.RDS',
#variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS.RDS',bb_lu=BB_LU,
#ptitle = 'PSA Beta 2018 UKBB Complete',hc.method="complete"
#)

#plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/psa_ss_shrink_2018/',
#basis_file='/home/ob219/share/as_basis/GWAS/support/psa_ss_basis_gwas.RDS',
#variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',bb_lu=BB_LU,
#ptitle = 'PSA Shrinkage 2018 UKBB Complete',hc.method="complete"
#)



plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS.RDS',bb_lu=BB_LU,
ptitle = 'Beta 2018 UKBB Average',hc.method="average"
)

plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',bb_lu=BB_LU,
ptitle = 'Shrinkage 2018 UKBB Average',hc.method="average"
)


plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS.RDS',bb_lu=BB_LU,
ptitle = 'Beta 2018 UKBB Ward',hc.method="ward.D"
)

plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',bb_lu=BB_LU,
ptitle = 'Shrinkage 2018 UKBB Ward',hc.method="ward.D"
)

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


plot_bb_pca_hclust <- function(proj_dir,basis_file,variance_file,ptitle,p.adj.thresh=0.01,hc.method="complete"){
  basis <- readRDS(basis_file)
  basis.DT <- data.table(trait=rownames(basis$x),basis$x)
  adj.vars <- readRDS(variance_file)

  bb.files <- list.files(path=proj_dir,pattern="*.RDS",full.names=TRUE)
  bb.DT<-lapply(bb.files,function(bb){
    readRDS(bb)
  }) %>% rbindlist
  # problematic as two entries in phenotype file with same descriptor
  bb.DT <- bb.DT[trait!='vitamin.c.product',]

  ## need to be able to map trait names onto sample sizes !
  p.DT <- get_phenotype_annotation()

  bb.DT.m <- melt(bb.DT,id.vars='trait')
  bb.DT.m <- merge(bb.DT.m,p.DT,by.x='trait',by.y='phe')
  bb.DT.m <- merge(bb.DT.m,adj.vars,by.x='variable',by.y='pc')
  bb.DT.m[,variance:=((n1+n0)/(n1 * n0)) * mfactor]
  tmp <- basis.DT[trait=='control',] %>% t
  ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))
  bb.DT.m <- merge(bb.DT.m,ctrl.DT,by='variable')

  ## this seems to not work there appears to be significant inflation
  ## this could be due to the fact that we have simulated under the null whereas
  ## all other conditions may not be null and may contain associations
  bb.DT.m[,Z:=(value-control.loading)/sqrt(variance)]
  bb.DT.m[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
  #bb.DT.m[,p.adj:=p.adjust(p.value),by='variable']
  bb.DT.m[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']

  ## we can try and use the distribution to estimate the parameters we need

  #bb.DT.m[,Z:=(value-mean(value))/sd(value),by='variable']
  #bb.DT.m[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
  #bb.DT.m[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']

  traits.of.interest <- bb.DT.m[,list(sum(p.adj<p.adj.thresh)),by='trait'][V1!=0,]$trait
  traits.of.interest <- traits.of.interest[traits.of.interest!='unclassifiable']
  ## this messes up bb.DT.m
  #bb.DT.m[,value:=Z]

  ## not all pc's are equal
  #vexp <- summary(basis)[['importance']][2,]
  #vexp.DT <- data.table(variable=names(vexp),vexp=vexp)
  #bb.DT.m <- merge(bb.DT.m,vexp.DT,by='variable')
  #bb.DT.m[,value:=value * vexp]

  mat.DT <- dcast(bb.DT.m[trait %in% traits.of.interest,],trait~variable)
  mat <- as.matrix(mat.DT[,paste('PC',1:11,sep=''),with=FALSE])
  mat <- rbind(mat,rep(0,ncol(mat)))
  rownames(mat) <- c(mat.DT$trait,'control')
  dist(mat) %>% hclust(.,method=hc.method) %>% plot(.,main=ptitle)
}

#plot_bb_Z_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_noshrink_2018/',
#basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_noshrink_gwas.RDS',
#variance_file='/home/ob219/share/as_basis/GWAS/support/ss_no_shrink_av_june.RDS'
#)
par(mfrow=c(1,3))

plot_bb_pca_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',
ptitle='2018 UKBB PCA Loading complete',
p.adj.thresh=0.01,hc.method="complete"
)

plot_bb_pca_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',
ptitle='2018 UKBB PCA Loading average',
p.adj.thresh=0.01,hc.method="average"
)

plot_bb_pca_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',
ptitle='2018 UKBB PCA Loading Ward',
p.adj.thresh=0.01,hc.method="ward.D"
)
