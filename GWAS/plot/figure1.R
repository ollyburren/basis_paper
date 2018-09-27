library(ggplot2)
library(cowplot)
library(ggrepel)
library(data.table)
library(magrittr)


BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/basis_gwas.RDS'
BASIS_NOSHRINK_FILE <- '/home/ob219/share/as_basis/GWAS/support/basis_noshrink_gwas.RDS'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
VARIANCE_FILE <- '/home/ob219/share/as_basis/ichip/support/analytical_variances_june.RDS'
VARIANCE_NOSHRINK_FILE <- '/home/ob219/share/as_basis/GWAS/support/no_shrink_av_june.RDS'
BB_DIR_PROJ <- '/home/ob219/share/as_basis/GWAS/bb_projections/'
BB_DIR_PROJ_NOSHRINK <- '/home/ob219/share/as_basis/GWAS/bb_projections/noshrink/'



basis <- readRDS(BASIS_FILE)
basis.DT <- data.table(trait=rownames(basis$x),basis$x)
noshrink.basis <- readRDS(BASIS_NOSHRINK_FILE)

## load in biobank projections

bb.files <- list.files(path=BB_DIR_PROJ,pattern="*.RDS",full.names=TRUE)

bb.DT<-lapply(bb.files,function(bb){
  readRDS(bb)
}) %>% rbindlist


## bind basis traits to bb traits where they exist
bb_lu <- list(
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

bb.DT.filt <- bb.DT[trait %in% unlist(bb_lu),]
bb.unl <- unlist(bb_lu)
midx <- match(bb.DT.filt$trait,bb.unl)
bb.DT.filt[,label:=paste('BB',names(bb.unl)[midx],sep=':'),]
bb.DT.filt[,disease:=names(bb.unl)[midx]]
plot.DT <- rbind(bb.DT.filt,basis.DT[,c('label','disease'):=list(trait,trait)])
y.int <- plot.DT[trait=='control',]$PC1
x.int <- plot.DT[trait=='control',]$PC2

p1.var <- signif(summary(basis)[['importance']][2,]['PC1']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC1',.)
p2.var <- signif(summary(basis)[['importance']][2,]['PC2']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC2',.)


 adj.vars <- readRDS(VARIANCE_FILE)

computeEllipse <- function(a,b,x,y,n.points=100){
  rad <- seq(0,2 * pi,length.out=n.points)
  xe <- (a * cos(rad)) + x
  ye <- (b * sin(rad)) + y
  elipse.DT <- data.table(x=xe,y=ye)
}


## test code for figure 1 to draw

get95CI <- function(var.DT,spc,n,n1,n.points=100){
  factor <- n/(n1 * (n-n1))
  mfactor <- var.DT[pc==spc,]$mfactor
  sqrt(mfactor * factor) * 1.96
}

n.total.bb <- 360000
n.cases <- c(300,1000,5000)

ci.DT <- lapply(n.cases,function(nc){
  N <- n.total.bb-nc
   tmp <- computeEllipse(get95CI(adj.vars,'PC1',N,nc),get95CI(adj.vars,'PC2',N,nc),plot.DT[trait=='control',]$PC1,plot.DT[trait=='control',]$PC2,100)
   tmp[,ncases:=nc]
}) %>% rbindlist

## create some labels for these
labs.DT <- lapply(n.cases,function(nc){
  N <- n.total.bb-nc
  rad <- (5*pi)/4
  a <- get95CI(adj.vars,'PC1',N,nc)
  b <- get95CI(adj.vars,'PC2',N,nc)
  x <- plot.DT[trait=='control',]$PC1
  y <- plot.DT[trait=='control',]$PC2
  xe <- (a * cos(rad)) + x
  ye <- (b * sin(rad)) + y
  data.table(x=xe,y=ye,label=format(nc,big.mark=',',trim=TRUE))
}) %>% rbindlist



pp1 <- ggplot(plot.DT,aes(x=PC1,y=PC2,col=disease,label=label)) + geom_point() + geom_text_repel() +
#coord_cartesian(xlim=c(-0.06,0.06),ylim=c(-0.06,0.06),expand=TRUE) +
scale_color_discrete(guide=FALSE) + xlab(p1.var) + ylab(p2.var) +
geom_hline(yintercept = x.int,col='black',size=0.5,lty=2,alpha=0.5) +
geom_vline(xintercept = y.int,col='black',size=0.5,lty=2,alpha=0.5) +
geom_path(data=ci.DT,aes(x=x,y=y,group=ncases),inherit.aes=FALSE,lty=2,alpha=0.5) +
geom_text(data=labs.DT,aes(x=x,y=y,group=NULL,label=label),inherit.aes=FALSE,alpha=0.3,cex=3)


## This stuff to compute p.values based on sample size and
if(FALSE){
  bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
  pheno <- fread(bb_phenofile)
  setnames(pheno,names(pheno) %>% make.names)
  med <- pheno[grepl("20002\\_",Phenotype.Code) & Sex=='both_sexes',]
  ## load in phenotype file
  bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
  pheno <- fread(bb_phenofile)
  setnames(pheno,names(pheno) %>% make.names)
  med <- pheno[grepl("20002\\_",Phenotype.Code) & Sex=='both_sexes',]
  ## load in phenotype file

  P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/phenotypes.both_sexes.tsv')
  P<-P[,.(phenotype,variable_type,non_missing=n_non_missing,cases=n_cases,controls=n_controls)]
  med[,phe:=make.names(Phenotype.Description) %>% gsub("Non.cancer.illness.code..self.reported..","",.)]
  med <- med[,.(Phenotype.Code,phe)]

  med <- merge(P,med,by.x='phenotype',by.y='Phenotype.Code')[,.(phe,cases,controls)]
  fz <- merge(bb.DT.filt,med,by.x='trait',by.y='phe')[,c('cases','controls'):=list(as.numeric(cases),as.numeric(controls))]

  computeP <- function(val,mean,n1,n0,var){
    scale.var <- ((n1+n0)/(n1 * n0)) * var
    pnorm(abs((val-mean)/sqrt(scale.var)),lower.tail=FALSE) * 2
  }
  pc1.control <- plot.DT[trait=='control']$PC1
  pc1.vars <- adj.vars[pc=='PC1',]$mfactor
  pc2.control <- plot.DT[trait=='control']$PC2
  pc2.vars <- adj.vars[pc=='PC2',]$mfactor
  fz[,c('p.PC1','p.PC2'):=list(computeP(PC1,pc1.control,cases,controls,pc1.vars),computeP(PC2,pc2.control,cases,controls,pc2.vars))]
}


## NO SHRINK !!!


basis <- readRDS(BASIS_NOSHRINK_FILE)
basis.DT <- data.table(trait=rownames(basis$x),basis$x)

## load in biobank projections

bb.files <- list.files(path=BB_DIR_PROJ_NOSHRINK,pattern="*.RDS",full.names=TRUE)

bb.DT<-lapply(bb.files,function(bb){
  readRDS(bb)
}) %>% rbindlist


## bind basis traits to bb traits where they exist
bb_lu <- list(
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

bb.DT.filt <- bb.DT[trait %in% unlist(bb_lu),]
bb.unl <- unlist(bb_lu)
midx <- match(bb.DT.filt$trait,bb.unl)
bb.DT.filt[,label:=paste('BB',names(bb.unl)[midx],sep=':'),]
bb.DT.filt[,disease:=names(bb.unl)[midx]]
plot.DT <- rbind(bb.DT.filt,basis.DT[,c('label','disease'):=list(trait,trait)])
y.int <- plot.DT[trait=='control',]$PC1
x.int <- plot.DT[trait=='control',]$PC2

p1.var <- signif(summary(basis)[['importance']][2,]['PC1']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC1',.)
p2.var <- signif(summary(basis)[['importance']][2,]['PC2']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC2',.)


 adj.vars <- readRDS(VARIANCE_NOSHRINK_FILE)

computeEllipse <- function(a,b,x,y,n.points=100){
  rad <- seq(0,2 * pi,length.out=n.points)
  xe <- (a * cos(rad)) + x
  ye <- (b * sin(rad)) + y
  elipse.DT <- data.table(x=xe,y=ye)
}


## test code for figure 1 to draw

get95CI <- function(var.DT,spc,n,n1,n.points=100){
  factor <- n/(n1 * (n-n1))
  mfactor <- var.DT[pc==spc,]$mfactor
  sqrt(mfactor * factor) * 1.96
}

n.total.bb <- 360000
n.cases <- c(300,1000,5000)

ci.DT <- lapply(n.cases,function(nc){
  N <- n.total.bb-nc
   tmp <- computeEllipse(get95CI(adj.vars,'PC1',N,nc),get95CI(adj.vars,'PC2',N,nc),plot.DT[trait=='control',]$PC1,plot.DT[trait=='control',]$PC2,100)
   tmp[,ncases:=nc]
}) %>% rbindlist

## create some labels for these
labs.DT <- lapply(n.cases,function(nc){
  N <- n.total.bb-nc
  rad <- (5*pi)/4
  a <- get95CI(adj.vars,'PC1',N,nc)
  b <- get95CI(adj.vars,'PC2',N,nc)
  x <- plot.DT[trait=='control',]$PC1
  y <- plot.DT[trait=='control',]$PC2
  xe <- (a * cos(rad)) + x
  ye <- (b * sin(rad)) + y
  data.table(x=xe,y=ye,label=format(nc,big.mark=',',trim=TRUE))
}) %>% rbindlist

pp2 <- ggplot(plot.DT,aes(x=PC1,y=PC2,col=disease,label=label)) + geom_point() + geom_text_repel() +
#coord_cartesian(xlim=c(0,0.1),ylim=c(-0.4,-0.1),expand=TRUE) +
scale_color_discrete(guide=FALSE) + xlab(p1.var) + ylab(p2.var) +
geom_hline(yintercept = x.int,col='black',size=0.5,lty=2,alpha=0.5) +
geom_vline(xintercept = y.int,col='black',size=0.5,lty=2,alpha=0.5) +
geom_path(data=ci.DT,aes(x=x,y=y,group=ncases),inherit.aes=FALSE,lty=2,alpha=0.5) #+
#geom_text(data=labs.DT,aes(x=x,y=y,group=NULL,label=label),inherit.aes=FALSE,alpha=0.5,cex=5)


pp2zoom <- ggplot(plot.DT,aes(x=PC1,y=PC2,col=disease,label=label)) + geom_point() + geom_text_repel() +
coord_cartesian(xlim=c(0,0.1),ylim=c(-0.4,-0.1),expand=TRUE) +
scale_color_discrete(guide=FALSE) + xlab(p1.var) + ylab(p2.var) +
geom_hline(yintercept = x.int,col='black',size=0.5,lty=2,alpha=0.5) +
geom_vline(xintercept = y.int,col='black',size=0.5,lty=2,alpha=0.5) +
geom_path(data=ci.DT,aes(x=x,y=y,group=ncases),inherit.aes=FALSE,lty=2,alpha=0.5) +
geom_text(data=labs.DT,aes(x=x,y=y,group=NULL,label=label),inherit.aes=FALSE,alpha=0.5,cex=5)



plot_grid(pp2,pp1)
