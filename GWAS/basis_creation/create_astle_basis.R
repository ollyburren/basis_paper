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
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/astle_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/astle_basis_gwas.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/astle_manifest_gwas.tab'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/astle_av.RDS'

man.DT <- fread(SNP_MANIFEST_FILE)

## load data
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
basis.DT[,p.value:=as.numeric(p.value)]
basis.DT[p.value==0,p.value:=3e-308]
basis.DT[,or:=as.numeric(or)]

## compute various shrinkage methods and store
shrink.DT<-compute_shrinkage_metrics(basis.DT,method="quant")
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

sbi(pc.emp)

library(pheatmap)
pheatmap(pc.emp$x[,paste('PC',1:14,sep="")],main="Astle")

if(FALSE){
  pc.emp <- readRDS(BASIS_FILE)
  shrink.DT <- readRDS(SHRINKAGE_FILE)
}

## project on IMD

imd.DT <- get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE,trait_list=c('UC','CD','CEL','T1D','MS','PBC','PSC','RA','SLE','asthma'))
#imd.DT[p.value==0,p.value:=3e-308]
#imd.DT[or==1,or:=1.000001]
imd.mat.emp <- create_ds_matrix(imd.DT,shrink.DT,SHRINKAGE_METHOD)
imd.proj <- predict(pc.emp,newdata=imd.mat.emp)

rbind(pc.emp$x,imd.proj) %>% pheatmap

## Chris wants to see the results of using varimax rotation

keep.pc <- paste0('PC',1:14)

X <- pc.emp$rotation[,keep.pc]
vm <- varimax(X,normalize=FALSE) # I understand normalize equalises relative contributions from different components, which may not be what we want
x.pca <- pc.emp$x[,keep.pc]
x.rot <- x.pca %*% vm$rotmat
colnames(x.rot) <- keep.pc

x.proj.rot <- imd.proj[,keep.pc] %*% vm$rotmat

M <- rbind(x.rot,x.proj.rot)
colnames(M) <- paste('PC',1:ncol(M),sep='')
pheatmap(M)

## work out the variance associated with the projection of onto varimax
load <- unclass(vm$loadings)
w.DT <- data.table(pid=rownames(load),load)
analytical.vars <- compute_proj_var(man.DT,w.DT,shrink.DT,REF_GT_DIR,SHRINKAGE_METHOD,quiet=FALSE)
adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars)
#adj.vars <- split(adj.vars,names(adj.vars))
setkey(adj.vars,pc)
saveRDS(adj.vars,file='/home/ob219/share/as_basis/GWAS/support/astle_av_varimax.RDS')

control.DT <- data.table(PC=names( x.rot["control",]),control.loading= x.rot["control",])

## compute confidence intervals for IMDs

colnames(x.proj.rot)<-keep.pc
x.proj.rot.dt <- data.table(trait=rownames(x.proj.rot),x.proj.rot)
x.proj.rot.dt <- melt(x.proj.rot.dt,id.vars='trait')
trait.manifest <- fread(TRAIT_MANIFEST)
x.proj.rot.dt <- merge(x.proj.rot.dt,trait.manifest[,.(trait,n1=cases,n=cases+controls)])
x.proj.rot.dt <- merge(x.proj.rot.dt,adj.vars,by.x='variable',by.y='pc')
x.proj.rot.dt <- merge(x.proj.rot.dt,control.DT,by.x='variable',by.y='PC')
x.proj.rot.dt[,variance:=((log(n)-(log(n1) + log(n-n1)))+ log(mfactor)) %>% exp]
## add in control loading
x.proj.rot.dt[,Z:=(value-control.loading)/sqrt(variance)]
x.proj.rot.dt[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
x.proj.rot.dt[,p.adj:=p.adjust(p.value,method="bonferroni"),by='variable']
x.proj.rot.dt[,delta:=value-control.loading]
x.proj.rot.dt[,category:='IMD']

basis.DT <- data.table(trait=rownames(x.rot),x.rot) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
basis.DT[,Z:=sign(value)]
basis.DT[,p.adj:=1]
basis.DT[,control.loading:=NULL]
basis.DT[,category:='basis']

all.DT <- rbind(x.proj.rot.dt,basis.DT,fill=TRUE)


all.DT[,ci:=1.96 * sqrt(variance)]
all.DT[,c('lower','upper'):=list(delta-ci,delta+ci)]


forp <- function(pdat,pc,p.thresh=0.05,theme=NA){
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }
  dat <- pdat[variable==pc,]
  dat[,trait:=factor(trait,levels=dat[order(category,delta,decreasing=TRUE),]$trait)]
  ggplot(dat,aes(x=trait,y=delta,colour=category,alpha=p.adj<p.thresh)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc) + theme +
  xlab("Trait") + ylab("Change in basis score from control") + guides(alpha=FALSE) + scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.3))
}

pdf(file="~/tmp/astle_varimax_results.pdf",paper="a4r",width=14,height=8)
lapply(keep.pc,function(pc){
  forp(all.DT,pc=pc)
})
dev.off()

saveRDS(all.DT,"/home/ob219/share/as_basis/GWAS/RESULTS/24_04_19_astle_varimax_summary_results.RDS")
