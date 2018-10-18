library(devtools)
load_all('~/git/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'


shrink.DT <- readRDS(SHRINKAGE_FILE)
pc.emp <- readRDS(BASIS_FILE)
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
tmp <- basis.DT[trait=='control',] %>% t
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))
## project jia
proj.traits <- fread(TRAIT_MANIFEST)[grep("^myositis_myogen",trait),]
proj.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE,proj.traits$trait)
proj.mat.emp<-create_ds_matrix(proj.DT,shrink.DT,SHRINKAGE_METHOD)
pred.emp <- predict(pc.emp,newdata=proj.mat.emp)
pred.DT <- data.table(trait=rownames(pred.emp),pred.emp)
var.DT <- readRDS(VARIANCE_FILE)
pred.DT <- melt(pred.DT,id.vars=c('trait'))
pred.DT <- merge(pred.DT,proj.traits[,.(trait,n1=cases,n=cases+controls)],by.x='trait',by.y='trait')
pred.DT <- merge(pred.DT,var.DT,by.x='variable',by.y='pc')
pred.DT[,ci.95:=sqrt((n/(n1 * (n-n1))) * mfactor) * 1.96  ]
pred.DT[,c('ci.lo','ci.hi'):=list(value-ci.95,value+ci.95)]
pred.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))
pred.DT <- merge(pred.DT,ctrl.DT,by='variable')
pred.DT[,variance:=(n/(n1 * (n-n1))) * mfactor]
pred.DT[,Z:=(value-control.loading)/sqrt(variance)]
pred.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
#bb.DT.m[,p.adj:=p.adjust(p.value),by='variable']
pred.DT[,p.adj:=p.adjust(p.value),by='variable']
pred.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
pred.DT[,short.trait:=substr(trait,1,15),]
pd <- position_dodge(0.1)
ggplot(all,aes(x=variable,y=value-control.loading,group=trait,col=trait,pch=p.value<0.05)) + geom_point(position=pd,aes(size=-log10(p.value))) +
geom_line(position=pd)
all<-rbind(melt(basis.DT,id.vars='trait'),pred.DT,fill=TRUE)
ggplot(all[trait %in% c('control','myositis_myogen'),],aes(x=variable,y=value,group=trait,col=trait)) + geom_point(position=pd) +
geom_line(position=pd)
