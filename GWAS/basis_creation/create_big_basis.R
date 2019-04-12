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
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/big_ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/big_ss_basis_gwas.RDS'
#BASIS_NOSHRINK_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_noshrink_gwas.RDS'
#BASIS_BETA_FILE <- '/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS'
GWAS_DATA_DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/all_imputed_studies_filtered'
SNP_MANIFEST_FILE <-'/home/ob219/rds/hpc-work/as_basis/support_tab/imputed_as_basis_snp_support_uk10k.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/big_as_manifest_gwas.tab'
#VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'
#VARIANCE_NOSHRINK_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_no_shrink_av_june.RDS'
#VARIANCE_BETA_FILE <- '/home/ob219/share/as_basis/GWAS/support/beta_av_june.RDS'

if(FALSE){
  ## manifest still needs some filtration
  man.DT <- fread(SNP_MANIFEST_FILE)
  ## remove SNPs with MAF <1%
  man.DT[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
  man.DT <- man.DT[maf>0.01,]
  ## palindromic SNPs
  man.DT[,alleles:=paste(ref_a1,ref_a2,sep='/')]
  idx <- which(man.DT$alleles %in% c('A/T','T/A','G/C','C/G'))
  man.DT<-man.DT[-idx,]
  ## remove the MHC
  man.DT[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
  mhc.pid <- man.DT[chr==6 & between(pos,20e6,40e6),]$pid
  man.DT <- man.DT[!pid %in% mhc.pid,.(pid,ref_a1,ref_a2,ref_a1.af,ld.block)]
  write.table(man.DT,file='/home/ob219/rds/hpc-work/as_basis/support_tab/imputed_as_basis_snp_support_uk10k.tab',sep="\t",quote=FALSE,row.names=FALSE)
}
man.DT <- fread(SNP_MANIFEST_FILE)

## load data
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
## remove instances where p==1 or p==0
remove.pid <- basis.DT[(p.value==0 | p.value==1 | or==0),]$pid %>% unique
basis.DT <- basis.DT[!pid %in% remove.pid,]

## compute various shrinkage methods and store
shrink.DT<-compute_shrinkage_metrics(basis.DT)
saveRDS(shrink.DT,file=SHRINKAGE_FILE)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
saveRDS(pc.emp,file=BASIS_FILE)

pc.emp.big <- readRDS(BASIS_FILE)

library("cowplot")
library("ggrepel")

## plot both scree and biplot

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

big <- sbi(pc.emp.big)

## how does this compare with the snps  we have used to construct the main basis

#SNP_MANIFEST_FILE <-'/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_uk10k.tab'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
man.DT <- fread(SNP_MANIFEST_FILE)
## load data
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
## remove instances where p==1 or p==0
remove.pid <- basis.DT[(p.value==0 | p.value==1 | or==0),]$pid %>% unique
basis.DT <- basis.DT[!pid %in% remove.pid,]
## compute various shrinkage methods and store
shrink.DT<-compute_shrinkage_metrics(basis.DT)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp.small <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
small <- sbi(pc.emp.small)

plot_grid(big,small,nrow=2)

#similar although PSC has swapped sides 
## better to look at hclust across all components

big.mat <- pc.emp.big$x
small.mat <- pc.emp.small$x
par(mfrow=c(1,2))
pcs <- paste('PC',1:6,sep="")
dist(big.mat[,pcs]) %>% hclust %>% plot(.,main="BIG")
dist(small.mat[,pcs]) %>% hclust %>% plot(.,main="SMALL")
library(pheatmap)
par(mfrow=c(1,2))
p1 <- dist(big.mat) %>% as.matrix %>% pheatmap
p2 <- dist(small.mat) %>% as.matrix %>% pheatmap

g <- grid.arrange(arrangeGrob(grobs= list(p1[[4]],p2[[4]]),ncol=2))


## plot shrinkages against each other

shrink.small.DT <- shrink.DT
shrink.big.DT <- readRDS(SHRINKAGE_FILE)
M <- merge(shrink.small.DT[,.(pid,small=ws_ppi,small.ld=ld.block)],shrink.big.DT[,.(pid,big=ws_ppi)],by='pid')
big.block<-shrink.big.DT[ws_ppi>0.01,]$ld.block %>% unique
small.block<-shrink.small.DT[ws_ppi>0.01,]$ld.block %>% unique
length(small.block)
length(big.block)
intersect(small.block,big.block) %>% length
setdiff(small.block,big.block) %>% length
setdiff(big.block,small.block) %>% length

big.sum <- shrink.big.DT[,list(count=.N,thresh.count=sum(ws_ppi>0.01)),by=ld.block]
big.sum[,norm.count:=thresh.count/count]
small.sum <- shrink.small.DT[,list(count=.N,thresh.count=sum(ws_ppi>0.01)),by=ld.block]
small.sum[,norm.count:=thresh.count/count]

big.sum[ld.block %in% setdiff(big.block,small.block),]$thresh.count %>% table
small.sum[ld.block %in% setdiff(small.block,big.block),]$thresh.count %>% table
