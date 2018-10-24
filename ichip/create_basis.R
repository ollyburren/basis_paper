## code to create and store basis and project all summary statistics

#library(devtools)
#install_github('ollyburren/cupcake')
load_all('~/git/cupcake')
library(cupcake)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/ichip/support/shrinkage_ic.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/ichip/support/basis_ic.RDS'
ICHIP_DATA_DIR <- '/home/ob219/share/as_basis/ichip/sum_stats'
SNP_MANIFEST_FILE <- '/home/ob219/share/as_basis/ichip/snp_manifest/ichip_all.tab'
MANIFEST <- '/home/ob219/share/as_basis/ichip/trait_manifest/as_manifest_ichip.tsv'
VARIANCE_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/analytical_variances_ichip.RDS'
REF_GT_DIR <- '/home/ob219/share/as_basis/ichip/ctrl_gt/by.chr+psa'
VARIANCE_FILE <- '/home/ob219/share/as_basis/ichip/support/shrink_av_ichip.RDS'

man.DT <- fread(SNP_MANIFEST_FILE)

## load data
basis.DT<-get_gwas_data(MANIFEST,SNP_MANIFEST_FILE,ICHIP_DATA_DIR,filter_snps_by_manifest=TRUE)
## compute various shrinkage methods and store
shrink.DT<-compute_shrinkage_metrics(basis.DT)
saveRDS(shrink.DT,file=SHRINKAGE_FILE)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
saveRDS(pc.emp,file=BASIS_FILE)

w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
analytical.vars <- compute_proj_var(man.DT,w.DT,shrink.DT,REF_GT_DIR,SHRINKAGE_METHOD,quiet=FALSE)
adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars)
#adj.vars <- split(adj.vars,names(adj.vars))
setkey(adj.vars,pc)
saveRDS(adj.vars,file=VARIANCE_FILE)

if(FALSE){
  library(cowplot)
  library(ggrepel)
  pc.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
  ggplot(pc.DT,aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text_repel()
}


## compute the variances
gw.DT <- basis.DT[trait==head(sample(unique(trait)),n=1),]
## set the standard error to the null
gw.DT[,emp_se:=se_null(n,n1,maf)]
## data table of PC snp loadings
gw.DT[,c('chr','position'):=tstrsplit(pid,':')]
w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
