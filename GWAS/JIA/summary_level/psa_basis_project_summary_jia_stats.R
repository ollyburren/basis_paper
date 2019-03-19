library(devtools)
#load_all('~/git/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_psa_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/psa_ss_basis_gwas.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_psa.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/psa_manifest_gwas.tab'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/psa_ss_av_june.RDS'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/jia_projections/summary/psa_jia.RDS'

shrink.DT <- readRDS(SHRINKAGE_FILE)
pc.emp <- readRDS(BASIS_FILE)
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
tmp <- basis.DT[trait=='control',] %>% t
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))
## project jia
proj.traits <- fread(TRAIT_MANIFEST)[grep("^jia",trait),]
proj.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE,proj.traits$trait)
proj.mat.emp<-create_ds_matrix(proj.DT,shrink.DT,SHRINKAGE_METHOD)
pred.emp <- predict(pc.emp,newdata=proj.mat.emp)
pred.DT <- data.table(trait=rownames(pred.emp),pred.emp)
#saveRDS(pred.DT,file=OUT_DIR)
## next compute projection 95% confidence intervals
var.DT <- readRDS(VARIANCE_FILE)
pred.DT <- melt(pred.DT,id.vars=c('trait'))
pred.DT <- merge(pred.DT,proj.traits[,.(trait,n1=cases,n=cases+controls)],by.x='trait',by.y='trait')
pred.DT <- merge(pred.DT,var.DT,by.x='variable',by.y='pc')
pred.DT[,ci.95:=sqrt((n/(n1 * (n-n1))) * mfactor) * 1.96  ]
pred.DT[,variable:=factor(variable,levels=paste0('PC',1:12))]
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))
pred.DT <- merge(pred.DT,ctrl.DT,by='variable')
pred.DT[,delta:=value-control.loading]
pred.DT[,c('ci.lo','ci.hi'):=list(delta-ci.95,delta+ci.95)]
pred.DT[,variance:=(n/(n1 * (n-n1))) * mfactor]
pred.DT[,Z:=(value-control.loading)/sqrt(variance)]
pred.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
#bb.DT.m[,p.adj:=p.adjust(p.value),by='variable']
pred.DT[,p.adj:=p.adjust(p.value),by='variable']
pred.DT[,variable:=factor(variable,levels=paste0('PC',1:12))]
pred.DT[,short.trait:=substr(trait,1,15),]

saveRDS(pred.DT,file="~/share/as_basis/GWAS/tmp/psa_jia_plot.RDS")
pd <- position_dodge(0.1)
pa <- ggplot(pred.DT,aes(x=variable,y=delta,group=trait,col=trait)) + geom_point(position=pd) +
geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi), width=.1, position=pd) + geom_line(position=pd) +
ylab(expression(Delta*"Control Loading")) + geom_hline(yintercept=0,color="black")
## here we use 0.05 rather than doing BF over all tests for all PC's
pred.DT[,Subtype:=gsub("jia\\_","",trait)]
#pb <- ggplot(pred.DT,aes(x=variable,y=value-control.loading,group=Subtype,col=Subtype,pch=p.adj<0.05)) + geom_point(position=pd,aes(size=-log10(p.value))) +
pb <- ggplot(pred.DT,aes(x=variable,y=value-control.loading,group=Subtype,col=Subtype,pch=p.adj<0.05)) + geom_point(size=2,position=pd) +
geom_line(position=pd) + ylab(expression(Delta*"Control Loading")) + xlab("Principal Component") + geom_hline(yintercept=0,color="black") +
background_grid(major = "xy", minor = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

save_plot(pb,file="~/tmp/line_jia.pdf",base_aspect=1.3)

forhc <- melt(pred.DT[,.(Subtype,pc=variable,Z,value)],id.vars=c('Subtype','pc'),measure.vars=c('Z','value'))

## do hclust on Z first "Z",colours=c('green','white','red')

bZ <- dcast(forhc[variable=='Z'],Subtype~pc)
bZm <- bZ[,-1] %>% as.matrix
rownames(bZm) <- bZ$Subtype
hc <- dist(bZm) %>% hclust
pred.DT[,Subtype:=factor(Subtype,levels=hc$labels[hc$order])]
pred.DT[,plotZ:=0]
pred.DT[p.adj<0.05,plotZ:=Z]
pc <- ggplot(pred.DT,aes(x=variable,y=Subtype,fill=plotZ,label=signif(value,digits=1))) +
geom_tile(color='black') + geom_text(cex=2.1) + scale_fill_gradient2("Z score") + xlab("Principal Component") +
ylab("Trait") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

save_plot(pc,file="~/tmp/heatmap_jia.pdf",base_aspect=1.3)

## plot the cluster

pdf("~/tmp/hclust_jia.pdf")
plot(hc)
dev.off()

plot_grid(pb,pc)


## compare with hinks et al.

hinks.DT <- fread("~/tmp/hinks_hla_correlation.csv")
hinks.m <- as.matrix(hinks.DT[,-1])
rownames(hinks.m) <- hinks.DT$trait
dist(hinks.m) %>% hclust(.,method="ward.D2") %>% plot


#save_plot("~/tmp/jia_2018.pdf",pb)
dev.print(pdf,"~/tmp/jia_2018.pdf")



## do hclust


## previously we had also plotted two bb traits also psoriasis and psoratic arthritis
## which I still need to do

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

#p.DT <- get_phenotype_annotation(filter='20002_1477|20002_1453|20002_1464')
p.DT <- get_phenotype_annotation()
bb.files <- list.files(path='/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018/',pattern="*.RDS",full.names=TRUE)
bb.DT<-lapply(bb.files,function(bb){
  readRDS(bb)
}) %>% rbindlist
bb.DT.m <- melt(bb.DT,id.vars='trait')
bb.DT.m <- merge(bb.DT.m,p.DT[,.(trait=phe,n1=n1,n=n1+n0)],by.x='trait',by.y='trait')
bb.DT.m <- merge(bb.DT.m,var.DT,by.x='variable',by.y='pc')
bb.DT.m[,ci.95:=sqrt((n/(n1 * (n-n1))) * mfactor) * 1.96  ]
bb.DT.m [,c('ci.lo','ci.hi'):=list(value-ci.95,value+ci.95)]
bb.DT.m [,variable:=factor(variable,levels=paste0('PC',1:11))]
bb.DT.m <- merge(bb.DT.m,ctrl.DT,by='variable')
bb.DT.m [,variance:=(n/(n1 * (n-n1))) * mfactor]
bb.DT.m [,Z:=(value-control.loading)/sqrt(variance)]
bb.DT.m[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
#bb.DT.m[,p.adj:=p.adjust(p.value),by='variable']
bb.DT.m[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']
bb.DT.m[,variable:=factor(variable,levels=paste0('PC',1:11))]
bb.DT.m[,short.trait:=substr(trait,1,15),]

## we can try and use the distribution to estimate the parameters we need

#bb.DT.m[,Z:=(value-mean(value))/sd(value),by='variable']
#bb.DT.m[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
#bb.DT.m[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']



traits.of.interest <- bb.DT.m[,list(sum(p.adj<0.01)),by='trait'][V1!=0,]$trait
bb.DT.m <- bb.DT.m[trait %in% traits.of.interest,]
## drop unclassifiable - what is this ?

## join with jia stuff and save as an RDS as useful for pooled variance approach
saveRDS(rbind(bb.DT.m,pred.DT),file="/home/ob219/share/as_basis/GWAS/tmp/jia_bb_summary.RDS")



z.DT <- melt(rbind(pred.DT,bb.DT.m),id.vars=c('variable','trait'),measure.vars='Z')

mat.DT <- dcast(z.DT ,trait~variable)
mat <- as.matrix(mat.DT[,paste('PC',1:11,sep=''),with=FALSE])
mat <- rbind(mat,rep(0,ncol(mat)))
rownames(mat) <- c(mat.DT$trait,'control')
dist(mat) %>% hclust %>% plot


groups <- dist(mat) %>% hclust %>% cutree(h=10)
groups <- split(names(groups),groups)
## plot loadings for group 11 and group 7
all <- rbind(bb.DT.m,pred.DT)

## plot a line plot of cluster to see which PC's are driving the cluster.
g1 <- ggplot(all[trait %in% c(groups[['7']],'rheumatoid.arthritis'),],aes(x=variable,y=value-control.loading,group=trait,col=short.trait,pch=p.value<0.05)) + geom_point(position=pd,aes(size=-log10(p.value))) +
geom_line(position=pd) + guides(size=FALSE)
g2 <- ggplot(all[trait %in% groups[['11']],],aes(x=variable,y=value-control.loading,group=trait,col=short.trait,pch=p.value<0.05)) + geom_point(position=pd,aes(size=-log10(p.value))) +
geom_line(position=pd) + guides(size=FALSE)
plot_grid(g1,g2,nrow=2)

## compare between PC's to see which are significantly different between the two clusters

c1 <- all[trait %in% groups[['7']],]
c1 <- split(c1$value-c1$control.loading,c1$variable)
c2 <- all[trait %in% groups[['11']],]
c2 <- split(c2$value-c2$control.loading,c2$variable)

res <- lapply(names(c2),function(x){
  tres <- t.test(c1[[x]],c2[[x]])
  data.table(p.value=tres$p.value,pc=x)
}) %>% rbindlist
