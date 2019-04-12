## code to create and store basis and project all summary statistics

library(devtools)
#install_github('ollyburren/cupcake')
load_all('~/git/cupcake')
library(cupcake)

## attempts to ask the question as to what threshold on the shrinkage
## causes the PCA to break down. Here we start with the full basis
## and then remove SNPs based on shrinkage score and see how this effects
## basis trait positions

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'

man.DT <- fread(SNP_MANIFEST_FILE)

## load data
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
## compute various shrinkage methods and store
shrink.DT<-compute_shrinkage_metrics(basis.DT)

quants <- quantile(shrink.DT$ws_emp_shrinkage,probs=seq(0,1,0.01))
quants <- quants[-1]
#quants <- quants[-length(quants)]
library(parallel)
all.res <- mclapply(seq_along(quants),function(i){
  message(quants[i])
  tshrink <- shrink.DT[ws_emp_shrinkage<=quants[i],]
  tbasis <- basis.DT[pid %in% tshrink$pid,]
  basis.mat.emp <- create_ds_matrix(tbasis,tshrink,SHRINKAGE_METHOD)
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
  rdt <- data.table(trait=rownames(pc.emp$x),centile=i,pc.emp$x)
},mc.cores=8) %>% rbindlist

M <- melt(all.res,id.vars=c('trait','centile'))
## compare to actual basis
act.basis <- readRDS(BASIS_FILE)
act.basis <- data.table(trait=rownames(act.basis$x),act.basis$x)
act.basis <- melt(act.basis,id.vars='trait')
setnames(act.basis,'value','act.basis.value')
Mb <- merge(M,act.basis,by=c('trait','variable'))
Mb[,delta:=abs(act.basis.value-value)]

library(cowplot)

ggplot(Mb[trait!='control' & centile>95],aes(x=centile,y=delta,colour=trait,group=trait)) + geom_point() + geom_line() + facet_wrap(~variable)

## all action is in the last centile so look carefully at these

quants <- quantile(shrink.DT$ws_emp_shrinkage,probs=seq(0.99,1,0.0001))
quants <- quants[-1]

all.res <- mclapply(seq_along(quants),function(i){
  message(quants[i])
  tshrink <- shrink.DT[ws_emp_shrinkage<=quants[i] & ws_emp_shrinkage>6.759033e-04,]
  tbasis <- basis.DT[pid %in% tshrink$pid,]
  basis.mat.emp <- create_ds_matrix(tbasis,tshrink,SHRINKAGE_METHOD)
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
  rdt <- data.table(trait=rownames(pc.emp$x),centile=i,pc.emp$x)
},mc.cores=8) %>% rbindlist

M <- melt(all.res,id.vars=c('trait','centile'))
## compare to actual basis
act.basis <- readRDS(BASIS_FILE)
act.basis <- data.table(trait=rownames(act.basis$x),act.basis$x)
act.basis <- melt(act.basis,id.vars='trait')
setnames(act.basis,'value','act.basis.value')
Mb <- merge(M,act.basis,by=c('trait','variable'))
Mb[,delta:=abs(act.basis.value-value)]

library(cowplot)

ggplot(Mb[trait!='control' & centile>75],aes(x=centile,y=delta,colour=trait,group=trait)) + geom_point() + geom_line() + facet_wrap(~variable)

## a cut off is ws_emp_shrinkage > 0.01 (which similar to what we have used previously)

tshrink <- shrink.DT[ws_emp_shrinkage>=0.01,]
tbasis <- basis.DT[pid %in% tshrink$pid,]
basis.mat.emp <- create_ds_matrix(tbasis,tshrink,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
rdt <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
M <- melt(rdt,id.vars=c('trait'))
act.basis <- readRDS(BASIS_FILE)
act.basis <- data.table(trait=rownames(act.basis$x),act.basis$x)
act.basis <- melt(act.basis,id.vars='trait')
setnames(act.basis,'value','act.basis.value')
Mb <- merge(M,act.basis,by=c('trait','variable'))
Mb[,delta:=abs(act.basis.value-value)]

ggplot(Mb[!variable %in% c('PC8','PC9'),],aes(x=trait,y=variable,fill=delta)) + geom_tile()
ggplot(Mb,aes(x=trait,y=variable,fill=delta)) + geom_tile()

library(cowplot)
library(ggrepel)
ggplot(Mb[trait!='control' & centile>75],aes(x=centile,y=delta,colour=trait,group=trait)) + geom_point() + geom_line() + facet_wrap(~variable)
#tshrink <- shrink.DT[ws_emp_shrinkage>0.0110710793,]
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
