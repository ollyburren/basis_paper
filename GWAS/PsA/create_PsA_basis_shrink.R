## process PsA data
## N appears to be 12,800
## N1 appears to be 3609
## N0 appears to be 9192

SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
PsA.dir <- '/home/ob219/share/Data/GWAS-summary/psa-2019-unpublished'
psa.files <- list.files(path=PsA.dir,pattern="*.out",full.names=TRUE)


DT <- lapply(psa.files,fread) %>% rbindlist
#DT <- fread(psa.files[1])
## what prop cases and controls
DT[,ctrl_a1_maf:=(controls_AA*2 + controls_AB)/(2*controls_total)]


DT[,list(ncases=(cases_AA + cases_AB + cases_BB)/2,nctrl=(controls_AA + controls_AB + controls_BB)/2)]
## have we flipped things
# from GWAS catalog rs76956521-C allele is the risk allele


DT.f <- DT[,.(rsid,chr=chromosome,pos=position,a1=alleleA,a2=alleleB,a1_maf=ctrl_a1_maf,p.value=frequentist_add_pvalue,beta=frequentist_add_beta_1,se.beta=frequentist_add_se_1)]
DT.f[,pid:=paste(chr,pos,sep=':')]
## filter for basis variants
snp.DT <- fread(SNP_MANIFEST_FILE)
## looks as if they are all here
#snp.DT[,c('chrm','positionm'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
#M <- merge(snp.DT[chrm==1,],DT.f,by='pid',all.x=TRUE)
M <- merge(snp.DT,DT.f,by='pid')
## there are a few SNPs missing so build a new manifest
out <- M[order(chr,pos),.(pid,a1,a2,or=exp(beta),p.value)]
write.table(out,file="/home/ob219/share/as_basis/GWAS/sum_stats/psa_bowes.tab",sep="\t",row.names=FALSE,quote=FALSE)
## next build a new basis with all basis traits and psa
library(cupcake)

## create a new snp.manifest to use

new.snp.DT <- snp.DT[pid %in% out$pid,]
write.table(new.snp.DT,file="/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_psa.tab",col.names=TRUE,row.names=FALSE,quote=FALSE)
## create a new trait.manifest to use

library(cupcake)
SHRINKAGE_METHOD<-'ws_emp_shrinkage'
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_psa_shrinkage_gwas.RDS'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/psa_manifest_gwas.tab'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_psa.tab'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
man.DT <- fread(SNP_MANIFEST_FILE)
trait.DT <- fread(TRAIT_MANIFEST)

basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
saveRDS(shrink.DT,file=SHRINKAGE_FILE)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
saveRDS(pc.emp,file="/home/ob219/share/as_basis/GWAS/support/psa_ss_basis_gwas.RDS")

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

psa <- sbi(pc.emp)

## compute variance with projection onto basis

w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
analytical.vars <- compute_proj_var(man.DT,w.DT,shrink.DT,REF_GT_DIR,SHRINKAGE_METHOD,quiet=FALSE)
adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars)
#adj.vars <- split(adj.vars,names(adj.vars))
setkey(adj.vars,pc)
saveRDS(adj.vars,file='/home/ob219/share/as_basis/GWAS/support/psa_ss_av_june.RDS')



##
