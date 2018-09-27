## code to create and store basis and project all summary statistics

library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/basis_gwas.RDS'
BASIS_NOSHRINK_FILE <- '/home/ob219/share/as_basis/GWAS/support/basis_noshrink_gwas.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/av_june.RDS'
VARIANCE_NOSHRINK_FILE <- '/home/ob219/share/as_basis/GWAS/support/no_shrink_av_june.RDS'

## load data
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
## compute various shrinkage methods and store
shrink.DT<-compute_shrinkage_metrics(basis.DT)
saveRDS(shrink.DT,file=SHRINKAGE_FILE)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
saveRDS(pc.emp,file=BASIS_FILE)
## compute the variances
man.DT <- fread(SNP_MANIFEST_FILE)
w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
analytical.vars <- compute_proj_var(man.DT,w.DT,shrink.DT,REF_GT_DIR,SHRINKAGE_METHOD,quiet=FALSE)
adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars)
#adj.vars <- split(adj.vars,names(adj.vars))
setkey(adj.vars,pc)
saveRDS(adj.vars,file=VARIANCE_FILE)

## create a basis just using gamma_hat - need this for figure 1
## for gamma hat log(or) * 1/ss_emp_maf_se so add this to shrinkage
#shrink.DT[,recip.ss_emp_maf_se:=1/ss_emp_maf_se]
basis.mat.noshrink <- create_ds_matrix(basis.DT,shrink.DT,'recip.ss_emp_maf_se')
basis.mat.noshrink <-rbind(basis.mat.noshrink,control=rep(0,ncol(basis.mat.noshrink)))
pc.emp.noshrink <- prcomp(basis.mat.noshrink,center=TRUE,scale=FALSE)
saveRDS(pc.emp.noshrink,file=BASIS_NOSHRINK_FILE)
w.DT.noshrink <- data.table(pid=rownames(pc.emp.noshrink$rotation),pc.emp.noshrink$rotation)
analytical.vars.noshrink <- compute_proj_var(man.DT,w.DT.noshrink,shrink.DT,REF_GT_DIR,'recip.ss_emp_maf_se',quiet=FALSE)
adj.vars.noshrink <- data.table(pc=names(analytical.vars.noshrink),mfactor=analytical.vars.noshrink)
setkey(adj.vars.noshrink,pc)
saveRDS(adj.vars.noshrink,file=VARIANCE_NOSHRINK_FILE)





pc.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
N <- 50000
N1 <- 200
factor <- N/(N1 * (N-N1))

computeEllipse <- function(a,b,x,y,n.points=1000){
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


ctrl.ellipse <- computeEllipse(get95CI(adj.vars,'PC1',50000,200),get95CI(adj.vars,'PC2',50000,200),pc.DT[trait=='control',]$PC1,pc.DT[trait=='control',]$PC2,100)

#ctrl.ellipse <- computeEllipse(sqrt(adj.vars[1,]$mfactor * factor),sqrt(adj.vars[2,]$mfactor * factor),pc.DT[trait=='control',]$PC1,pc.DT[trait=='control',]$PC2,100)

library(ggplot2)
library(cowplot)
library(ggrepel)


ggplot(pc.DT,aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text_repel() +
geom_path(data=ctrl.ellipse,aes(x=x,y=y),inherit.aes=FALSE)
