## code to create and store basis and project all summary statistics

library(devtools)
#install_github('ollyburren/cupcake')
load_all('~/git/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/psa_ss_shrinkage_myositis_gwas.RDS'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/psa_manifest_gwas.tab'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/psa_gwas_myositis.tab'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/psa_myositis_ss_av_june.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/psa_ss_basis_myositis_gwas.RDS'


if(FALSE){
  ## load in myogen data in order to create a basis snp manifest that only has myogen
  ## available SNPs in
  jdm_myo <- readRDS("/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/myogen_myositis/jdm_myositis_source.RDS")
  jdm_myo[is.na(or)]
  ## load in current snp support file and limit by jdm snps
  snps <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_psa.tab')
  snps_dt <- snps[pid %in% jdm_myo[!is.na(or),]$pid,]
  write.table(snps_dt,file='/home/ob219/share/as_basis/GWAS/snp_manifest/psa_gwas_myositis.tab',sep="\t",row.names=FALSE)
}


man.DT <- fread(SNP_MANIFEST_FILE)
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

w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
analytical.vars <- compute_proj_var(man.DT,w.DT,shrink.DT,REF_GT_DIR,SHRINKAGE_METHOD,quiet=FALSE)
adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars)
#adj.vars <- split(adj.vars,names(adj.vars))
setkey(adj.vars,pc)
saveRDS(adj.vars,file=VARIANCE_FILE)

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

psa_myogen <- sbi(pc.emp)
