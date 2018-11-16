library(ggplot2)
library(cowplot)
library(ggrepel)
library(data.table)
library(magrittr)

plot_bb_pc <- function(pca1='PC1',pca2='PC2',proj_dir,basis_file,variance_file,bb_lu,lims){

  ## currently hard coded size of biobank
  n.total.bb <- 360000
  ## where to draw ellipses of 95% CI for GWAS with that sample size
  n.cases <- c(300,1000,5000)
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
  bb.DT.filt[,label:=paste('bb',names(bb.unl)[midx],sep='_'),]
  bb.DT.filt[,disease:=names(bb.unl)[midx]]
  plot.DT <- rbind(bb.DT.filt,basis.DT[,c('label','disease'):=list(trait,trait)])
  y.int <- plot.DT[trait=='control',][[pca1]]
  x.int <- plot.DT[trait=='control',][[pca2]]
  ## get the other traits so we can plot on top (or behind)
  bb.DT.back <- bb.DT[!trait %in% unlist(bb_lu),]
  p1.var <- signif(summary(basis)[['importance']][2,]['PC1']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC1',.)
  p2.var <- signif(summary(basis)[['importance']][2,]['PC2']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC2',.)
  adj.vars <- readRDS(variance_file)

  computeEllipse <- function(a,b,x,y,n.points=100){
    rad <- seq(0,2 * pi,length.out=n.points)
    xe <- (a * cos(rad)) + x
    ye <- (b * sin(rad)) + y
    elipse.DT <- data.table(x=xe,y=ye)
  }

  get95CI <- function(var.DT,spc,n,n1,n.points=100){
    factor <- n/(n1 * (n-n1))
    mfactor <- var.DT[pc==spc,]$mfactor
    sqrt(mfactor * factor) * 1.96
  }

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
  pp <- ggplot(plot.DT,aes(x=PC1,y=PC2,col=disease,label=label)) + geom_point() + geom_text_repel() +
  scale_color_discrete(guide=FALSE) + xlab(p1.var) + ylab(p2.var) +
  geom_hline(yintercept = x.int,col='black',size=0.5,lty=2,alpha=0.5) +
  geom_vline(xintercept = y.int,col='black',size=0.5,lty=2,alpha=0.5) +
  geom_path(data=ci.DT,aes(x=x,y=y,group=ncases),inherit.aes=FALSE,lty=2,alpha=0.5) +
  geom_text(data=labs.DT,aes(x=x,y=y,group=NULL,label=label),inherit.aes=FALSE,alpha=0.3,cex=3) +
  geom_point(data=bb.DT.back,aes(x=PC1,y=PC2),inherit.aes=FALSE,alpha=0.05)
  if(!missing(lims))
    pp <- pp + coord_cartesian(xlim=lims$x,ylim=lims$y,expand=TRUE)
  pp
}

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


shrink_ss <- plot_bb_pc(
  proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
  basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
  variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',bb_lu=BB_LU)

noshrink_ss <- plot_bb_pc(
    proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_noshrink_2018/',
    basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_noshrink_gwas.RDS',
    variance_file='/home/ob219/share/as_basis/GWAS/support/ss_no_shrink_av_june.RDS',bb_lu=BB_LU)

beta <- plot_bb_pc(
    proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
    basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
    variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS',bb_lu=BB_LU)

LIMS <- list(x=c(0,0.1),y=c(-0.4,-0.1))

noshrink_ss_zoom <- plot_bb_pc(
    proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_noshrink_2018/',
    basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_noshrink_gwas.RDS',
    variance_file='/home/ob219/share/as_basis/GWAS/support/ss_no_shrink_av_june.RDS',bb_lu=BB_LU,lims=LIMS)


save_plot("~/tmp/bb_beta_2018.pdf",beta)
save_plot("~/tmp/bb_shrink_2018.pdf",shrink_ss)

plot_grid(beta + ggtitle("Beta 2018 UKBB"),shrink_ss + ggtitle("Shrinkage 2018 UKBB"))


plot_pc <- function(pca1='PC1',pca2='PC2',proj_dir,basis_file,variance_file,bb_lu,lims){

  ## currently hard coded size of biobank
  n.total.bb <- 360000
  ## where to draw ellipses of 95% CI for GWAS with that sample size
  n.cases <- c(300,1000,5000)
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
  bb.DT.filt[,label:=paste('bb',names(bb.unl)[midx],sep='_'),]
  bb.DT.filt[,disease:=names(bb.unl)[midx]]
  plot.DT <- basis.DT[,c('label','disease'):=list(trait,trait)]
  y.int <- plot.DT[trait=='control',][[pca1]]
  x.int <- plot.DT[trait=='control',][[pca2]]
  ## get the other traits so we can plot on top (or behind)
  bb.DT.back <- bb.DT[!trait %in% unlist(bb_lu),]
  p1.var <- signif(summary(basis)[['importance']][2,]['PC1']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC1',.)
  p2.var <- signif(summary(basis)[['importance']][2,]['PC2']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC2',.)
  adj.vars <- readRDS(variance_file)

  computeEllipse <- function(a,b,x,y,n.points=100){
    rad <- seq(0,2 * pi,length.out=n.points)
    xe <- (a * cos(rad)) + x
    ye <- (b * sin(rad)) + y
    elipse.DT <- data.table(x=xe,y=ye)
  }

  get95CI <- function(var.DT,spc,n,n1,n.points=100){
    factor <- n/(n1 * (n-n1))
    mfactor <- var.DT[pc==spc,]$mfactor
    sqrt(mfactor * factor) * 1.96
  }

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
  pp <- ggplot(plot.DT,aes(x=PC1,y=PC2,col=disease,label=label)) + geom_point() + geom_text_repel() +
  scale_color_discrete(guide=FALSE) + xlab(p1.var) + ylab(p2.var) +
  geom_hline(yintercept = x.int,col='black',size=0.5,lty=2,alpha=0.5) +
  geom_vline(xintercept = y.int,col='black',size=0.5,lty=2,alpha=0.5) #+
  #geom_path(data=ci.DT,aes(x=x,y=y,group=ncases),inherit.aes=FALSE,lty=2,alpha=0.5) +
  #geom_text(data=labs.DT,aes(x=x,y=y,group=NULL,label=label),inherit.aes=FALSE,alpha=0.3,cex=3) #+
  #geom_point(data=bb.DT.back,aes(x=PC1,y=PC2),inherit.aes=FALSE,alpha=0.05)
  if(!missing(lims))
    pp <- pp + coord_cartesian(xlim=lims$x,ylim=lims$y,expand=TRUE)
  pp
}

beta <- plot_pc(
    proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
    basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
    variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS',bb_lu=BB_LU)

save_plot("~/tmp/beta_2018.pdf",beta)
## next the hclust stuff


plot_bb_hclust <- function(proj_dir,basis_file,variance_file,bb_lu,ptitle){
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
  dist(proj) %>% hclust %>% plot(.,main=ptitle)
}


par(mfrow=c(1,2))



#plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_noshrink_2018/',
#basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_noshrink_gwas.RDS',
#variance_file='/home/ob219/share/as_basis/GWAS/support/ss_no_shrink_av_june.RDS',bb_lu=BB_LU
#)

plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS.RDS',bb_lu=BB_LU,
ptitle = 'Beta 2018 UKBB'
)

plot_bb_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',bb_lu=BB_LU,
ptitle = 'Shrinkage 2018 UKBB'
)

dev.print(pdf,"~/tmp/figure1.pdf")

## use the analytical variance estimations to compute Z scores - take traits forward for hclust if they show
## significant loading on at least one principal component after multiple testing


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




plot_bb_Z_hclust <- function(proj_dir,basis_file,variance_file,ptitle,p.adj.thresh=0.05){
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
  ## this messes up bb.DT.m
  bb.DT.m[,value:=Z]

  ## not all pc's are equal
  #vexp <- summary(basis)[['importance']][2,]
  #vexp.DT <- data.table(variable=names(vexp),vexp=vexp)
  #bb.DT.m <- merge(bb.DT.m,vexp.DT,by='variable')
  #bb.DT.m[,value:=value * vexp]

  mat.DT <- dcast(bb.DT.m[trait %in% traits.of.interest,],trait~variable)
  mat <- as.matrix(mat.DT[,paste('PC',1:11,sep=''),with=FALSE])
  mat <- rbind(mat,rep(0,ncol(mat)))
  rownames(mat) <- c(mat.DT$trait,'control')
  dist(mat) %>% hclust %>% plot(.,main=ptitle)
}

#par(mfrow=c(1,1))

#plot_bb_Z_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
#basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
#variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS',
#ptitle='Beta 2018 UKBB'
#)

par(mfrow=c(1,1))

plot_bb_Z_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',
ptitle='Shrinkage 2018 UKBB'
)

plot_bb_Z_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',
ptitle='Shrinkage 2018 UKBB',
p.adj.thresh=0.01
)

dev.print(pdf,"~/tmp/figure1b.pdf")

proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/'
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'
bb_lu=BB_LU

plot_bb_pca_hclust <- function(proj_dir,basis_file,variance_file,ptitle,p.adj.thresh=0.01){
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
  dist(mat) %>% hclust %>% plot(.,main=ptitle)
}

#plot_bb_Z_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_noshrink_2018/',
#basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_noshrink_gwas.RDS',
#variance_file='/home/ob219/share/as_basis/GWAS/support/ss_no_shrink_av_june.RDS'
#)

plot_bb_pca_hclust(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',
ptitle='2018 UKBB PCA Loading',
p.adj.thresh=0.01
)

dev.print(pdf,"~/tmp/figure1bpca.pdf")

# plot for Chris


annotate_loading <- function(proj.DT,basis,var.DT){
  p.DT <- get_phenotype_annotation()
  basis.DT <- data.table(trait=rownames(basis$x),basis$x)
  proj.DT.m <- melt(proj.DT,id.vars='trait')
  proj.DT.m <- merge(proj.DT.m,p.DT,by.x='trait',by.y='phe')
  proj.DT.m <- merge(proj.DT.m,var.DT,by.x='variable',by.y='pc')
  proj.DT.m[,variance:=((n1+n0)/(n1 * n0)) * mfactor]
  ctrl.DT <- basis.DT[trait=='control',][,paste('PC',1:11,sep=""),with=FALSE] %>% t %>% data.table(variable=rownames(.),.)
  setnames(ctrl.DT,'V1','ctrl.loading')
  proj.DT.m <- merge(proj.DT.m,ctrl.DT,by='variable')
  proj.DT.m[,Z:=(value-ctrl.loading)/sqrt(variance)]
  proj.DT.m[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
  vexp.DT <- data.table(summary(basis)[['importance']][2,] %>% names,summary(basis)[['importance']][2,])
  setnames(vexp.DT,c('variable','variance.exp'))
  proj.DT.m <- merge(proj.DT.m,vexp.DT,by='variable')
  proj.DT.m[,.(variable,trait,code,pc.loading=value,n.cases=n1,Z,Z.exp=Z * variance.exp,unit.Z=value/sqrt(mfactor),p.value)]

}

proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/'
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'
BB_LU2=BB_LU[names(BB_LU) != unlist(BB_LU)]


plot_line_proj <- function(basis_file,variance_file,proj_dir,bb_lu){
  basis <- readRDS(basis_file)
  var.DT <- readRDS(variance_file)
  bb.files <- list.files(path=proj_dir,pattern="*.RDS",full.names=TRUE)
  proj.DT<-lapply(bb.files,function(bb){
    readRDS(bb)
  }) %>% rbindlist
  basis.DT <- data.table(trait=rownames(basis$x),basis$x)
  plot.DT <- rbind(proj.DT,basis.DT) %>% melt(.,id.vars='trait')
  #plot.DT.m <- merge(plot.DT,var.DT,by.x='variable',by.y='pc')[,value:=value/sqrt(mfactor)]
  plot.DT.m <- merge(plot.DT,var.DT,by.x='variable',by.y='pc')
  #anno.DT <- annotate_loading(rbind(basis.DT,proji.DT),basis,var.DT)
  plot.DT.m <- plot.DT.m[trait %in% unlist(bb_lu) | trait %in% names(bb_lu) | trait=='control',]
  bb.unl <- unlist(bb_lu)
  midx <- match(plot.DT.m$trait,bb.unl)
  plot.DT.m[,label:=paste('bb',names(bb.unl)[midx],sep='_'),]
  plot.DT.m[,disease:=names(bb.unl)[midx]]
  plot.DT.m[is.na(disease),c('label','disease'):=list(trait,trait)]
  plot.DT.m[,variable:=factor(variable,levels=paste0('PC',1:11))]
  plot.DT.m[,Type:=ifelse(grepl("^bb_",label),'UK Biobank projection','Basis')]
  #return(plot.DT.m)
  c.DT <-plot.DT[trait=='control',]
  ggplot(  plot.DT.m[trait!='control',],aes(x=variable,y=value,group=trait,color=Type)) + geom_point(size=2) +
  geom_line(size=1) +

  xlab("Principal Component") + ylab("Unit PC score") +
  facet_wrap(~disease) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

#plot_line_proj(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
#basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
#variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',
#bb_lu=BB_LU2) %>% saveRDS(.,file='/home/ob219/share/as_basis/GWAS/results/shrink_bb_projection.RDS')

#plot_line_proj(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
#basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
#variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS',
#bb_lu=BB_LU2)  %>% saveRDS(.,file='/home/ob219/share/as_basis/GWAS/results/beta_bb_projection.RDS')


fline_shrink <- plot_line_proj(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',
bb_lu=BB_LU2)

fline_beta <- plot_line_proj(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS',
bb_lu=BB_LU2)


plot_grid(fline_beta + ggtitle("UKBB 2018 Beta"),fline_shrink + ggtitle("UKBB 2018 Shrinkage"))


plot_line_z <- function(basis_file,variance_file,proj_dir,bb_lu){
  basis <- readRDS(basis_file)
  var.DT <- readRDS(variance_file)
  bb.files <- list.files(path=proj_dir,pattern="*.RDS",full.names=TRUE)
  proj.DT<-lapply(bb.files,function(bb){
    readRDS(bb)
  }) %>% rbindlist
  anno.DT <- annotate_loading(proj.DT,basis,var.DT)
  anno.DT.filt <- anno.DT[trait %in% unlist(bb_lu),]
  bb.unl <- unlist(bb_lu)
  midx <- match(anno.DT.filt$trait,bb.unl)
  anno.DT.filt[,label:=paste('bb',names(bb.unl)[midx],sep='_'),]
  anno.DT.filt[,disease:=names(bb.unl)[midx]]
  anno.DT.filt[,variable:=factor(variable,levels=paste0('PC',1:11))]

  ## add on a rectangle to show PC's that are significant
  n.disease <- anno.DT.filt$trait %>% unique %>% length
  ymax <-  qnorm(0.05/(2*n.disease),lower.tail=FALSE)
  ymin <- qnorm(0.05/(2*n.disease),lower.tail=FALSE) * -1
  anno.DT.filt[,Trait:=gsub("\\."," ",trait)]

  ggplot(anno.DT.filt,aes(x=variable,y=Z,group=Trait,color=Trait)) + geom_point(size=2) +
  geom_line(size=1) +
  annotate("rect", ymin=ymin,ymax=ymax,xmin=0,xmax=Inf, alpha=0.5, fill="black") +
  xlab("Principal Component") + ylab("Z score") + geom_hline(yintercept=0,size=1)
}

line_shrink <- plot_line_z(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS',
bb_lu=BB_LU)

line_beta <- plot_line_z(proj_dir='/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/',
basis_file='/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS',
variance_file='/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS',
bb_lu=BB_LU)

plot_grid(line_beta + ggtitle("UKBB 2018 Beta"),line_shrink + ggtitle("UKBB 2018 Shrinkage"),nrow=2)
