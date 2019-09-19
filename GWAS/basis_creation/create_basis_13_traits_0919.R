## code to create and store basis and project all summary statistics

library(devtools)
#install_github('ollyburren/cupcake')
load_all('~/git/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)


## when aligning to the basis the counted allele or the allele with which the or is respect to should
## is allele2. This is demonstrated by ptpn22
## A G 1.2 CD
## A G 0.5 T1D
## Here we know that minor allele A is protective from crohn's and risk for T1D therefore
## basis must be with respect to allele2

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas_13_traits_0919.tab'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_gwas_13_traits_0919.RDS'


## sort out lada data sets
if(FALSE){
lada <- readRDS("/home/ob219/share/as_basis/GWAS/for_fdr/cousminer_lada_source.RDS")
man.DT <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab')
M <- merge(lada,man.DT,by='pid')
M<-M[,.(pid,a1=ref_a1,a2=ref_a2,or,p.value)]
write.table(M[!is.na(or),],file=file.path(GWAS_DATA_DIR,'cousminer_lada.tab'),row.names=FALSE,quote=FALSE,sep="\t")
keep.lada <- M[!is.na(or),]$pid

igneph <- readRDS("/home/ob219/share/as_basis/GWAS/for_fdr/IgA_nephropathy_source.RDS")
M <- merge(igneph,man.DT,by='pid')
M<-M[,.(pid,a1=ref_a1,a2=ref_a2,or,p.value)]
write.table(M[!is.na(or),],file=file.path(GWAS_DATA_DIR,'kiryluk_neph.tab'),row.names=FALSE,quote=FALSE,sep="\t")
keep.neph <- M[!is.na(or),]$pid
man.DT <- man.DT[pid %in% intersect(keep.lada,keep.neph),]
write.table(man.DT,file=SNP_MANIFEST_FILE,row.names=FALSE,quote=FALSE,sep="\t")
}

man.DT <- fread(SNP_MANIFEST_FILE)
## load data
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
basis.DT[,p.value:=as.numeric(p.value)]
## compute various shrinkage methods and store ommit T2D from
## creating the weights
shrink.DT<-compute_shrinkage_metrics(basis.DT)
saveRDS(shrink.DT,file=SHRINKAGE_FILE)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)

sbi <- function(pc.emp){
  vexp <- summary(pc.emp)[['importance']][2,]
  PC1.var<-signif(vexp["PC1"]*100,digits=3)
  PC2.var<-signif(vexp["PC2"]*100,digits=3)
  M <- cbind(as.data.table(pc.emp$x),trait=rownames(pc.emp$x))
  scp <- cbind(data.table(vexp),pcs=factor(names(vexp),levels=names(vexp)))
  scp[,cs:=cumsum(vexp)]
  scp$group=1
  text.size <- 20
  ## do a scree plot
  ppl <- ggplot(scp,aes(x=pcs,y=vexp,group=group)) + geom_point() + geom_line() + ylab("Variance Explained") + xlab("Principal Components") +
  theme(axis.text=element_text(size=text.size),axis.title=element_text(size=text.size),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ppr<-ggplot(M,aes(x=PC1,y=PC2,label=trait)) + geom_point(size=3) + geom_text_repel(size=7) + # hjust = 0, nudge_x = 0.005)  +
  scale_color_discrete(guide=FALSE) + scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) +
  xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +  background_grid(major = "xy", minor = "none") +
  theme(axis.text=element_text(size=text.size),axis.title=element_text(size=text.size))
  plot_grid(ppl, ppr, labels = "auto",label_size=text.size)
}

sbi(pc.emp)

library(pheatmap)
pheatmap(pc.emp$x)

saveRDS(pc.emp,file=BASIS_FILE)
## compute the variances

w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
analytical.vars <- compute_proj_var(man.DT,w.DT,shrink.DT,REF_GT_DIR,SHRINKAGE_METHOD,quiet=FALSE)
adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars)
#adj.vars <- split(adj.vars,names(adj.vars))
setkey(adj.vars,pc)
saveRDS(adj.vars,file=VARIANCE_FILE)

## will be required for the paper
if(FALSE){
## create a basis just using gamma_hat - need this for figure 1
## for gamma hat log(or) * 1/ss_emp_maf_se so add this to shrinkage
#shrink.DT[,recip.ss_emp_maf_se:=1/ss_emp_maf_se]
#basis.mat.noshrink <- create_ds_matrix(basis.DT,shrink.DT,'recip.emp_maf_se')
basis.mat.noshrink <- create_ds_matrix(basis.DT,shrink.DT,'recip.ss_emp_maf_se')
basis.mat.noshrink <-rbind(basis.mat.noshrink,control=rep(0,ncol(basis.mat.noshrink)))
pc.emp.noshrink <- prcomp(basis.mat.noshrink,center=TRUE,scale=FALSE)
saveRDS(pc.emp.noshrink,file=BASIS_NOSHRINK_FILE)
w.DT.noshrink <- data.table(pid=rownames(pc.emp.noshrink$rotation),pc.emp.noshrink$rotation)
analytical.vars.noshrink <- compute_proj_var(man.DT,w.DT.noshrink,shrink.DT,REF_GT_DIR,'recip.ss_emp_maf_se',quiet=FALSE)
adj.vars.noshrink <- data.table(pc=names(analytical.vars.noshrink),mfactor=analytical.vars.noshrink)
setkey(adj.vars.noshrink,pc)
saveRDS(adj.vars.noshrink,file=VARIANCE_NOSHRINK_FILE)

## create basis using just beta

basis.mat.beta <- create_ds_matrix(basis.DT,shrink.DT,'none')
basis.mat.beta <-rbind(basis.mat.beta,control=rep(0,ncol(basis.mat.beta)))
pc.emp.beta <- prcomp(basis.mat.beta,center=TRUE,scale=FALSE)
saveRDS(pc.emp.beta,file=BASIS_BETA_FILE)
w.DT.beta <- data.table(pid=rownames(pc.emp.beta$rotation),pc.emp.beta$rotation)
analytical.vars.beta <- compute_proj_var(man.DT,w.DT.beta,shrink.DT,REF_GT_DIR,'none',quiet=FALSE)
adj.vars.beta <- data.table(pc=names(analytical.vars.beta),mfactor=analytical.vars.beta)
setkey(adj.vars.beta,pc)
saveRDS(adj.vars.beta,file=VARIANCE_BETA_FILE)
}
