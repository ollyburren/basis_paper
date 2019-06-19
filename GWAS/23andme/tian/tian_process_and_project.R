## process 23andme data


SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
DATA.DIR <- '/home/ob219/share/Data/GWAS-summary/aav_limy_wong'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- "/home/ob219/share/as_basis/GWAS/tian_projections/tian_infectious_disease_0619.RDS"

## Tian
anno.file <- '/home/ob219/share/Data/GWAS-summary/23andme/annotation/all_snp_info-4.1.txt'

anno.DT <- fread(anno.file)
anno.DT <- anno.DT[,.(assay.name,chr=sub("chr","",scaffold),position,alleles,all.data.id,name=assay.name)]
anno.DT[,pid:=paste(chr,position,sep=':')]

## load in SNP manifest

snp.man <- fread(SNP_MANIFEST)
anno.DT.f <- anno.DT[pid %in% snp.man$pid,]
## there are some duplicates which we need to resolve

anno.DT.f[,c('a1','a2'):=tstrsplit(alleles,'/')]
M <- merge(anno.DT.f,snp.man,by='pid')

## first off attempt to align alleles - at this stage not sure
## which is the effect allele

library(annotSnpStats)
alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=M$alleles)
#alleles <- alleles[!duplicated(pid),]
#alleles <- M[,list(al.x=paste(uk10_A1,uk10_A2,sep='/'),al.y=paste(a1,a2,sep='/')),by='pid']
## to make quick
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$pid
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
alleles[,g.class:=align.class]
idx<-which(alleles$g.class=='impossible')
M <- M[-idx,]
alleles <- alleles[-idx,]
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
tian.anno <- M[,.(all.data.id,pid,g.class,alleles)]


DATA_DIR <- '/home/ob219/share/Data/GWAS-summary/23andme/tian/'
files <- list.files(path=DATA_DIR,pattern="*.gz",full.names=TRUE)

library(parallel)
all.DT <- mclapply(files,function(file){
  message(file)
  DT <- sprintf("zcat %s",file) %>% fread
  ## only take forward those that pass
  #DT <- DT[pass=='Y',]
  DT <- merge(DT,tian.anno,by='all.data.id')
  DT[,or:=exp(effect)]
  DT[g.class=='rev',or:=exp(effect * -1)]
  trait <- basename(file) %>% sub("-4.1.dat.gz","",.)
  topr <- DT[,.(trait=trait,pid,or,n0=im.num.0,n1=im.num.1)]
  topr
},mc.cores=8) %>% rbindlist
setkey(all.DT,pid)
all.DT <- all.DT[!is.na(or),]

## stripped down code from cupcake to build matrix
## want to experiment with dcast fill=0 to help with missing
## variants

sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
tmp<-all.DT[stmp]
tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
B <- dcast(tmp,pid ~ trait,value.var='metric',fill=0)
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)


all.proj.DT <- data.table(trait=rownames(all.proj),all.proj)
saveRDS(all.proj.DT,file="/home/ob219/share/as_basis/GWAS/tian_projections/tian_infectious_disease_0619.RDS")


all.proj.m <- melt(all.proj.DT,id.var='trait')
all.proj.m[,variable:=factor(variable,levels=paste('PC',1:11,sep=''))]

library(ggplot2)
ggplot(all.proj.m,aes(x=variable,y=value,color=trait,group=trait)) +
geom_point() + geom_line()

VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'

var.DT <- readRDS(VARIANCE_FILE)

proj.traits <- all.DT[!is.na(n0) & !is.na(n1),list(controls=max(n0),cases=max(n1)),by='trait']
tian_samples <- all.DT[!is.na(n0) & !is.na(n1),list(n1=max(n1),n0=max(n0)),by='trait']
saveRDS(tian_samples,file="/home/ob219/share/as_basis/GWAS/tian_projections/sample_size.RDS")

pred.DT <- merge(all.proj.m,proj.traits[,.(trait,n1=cases,n=cases+controls)],by.x='trait',by.y='trait')
pred.DT <- merge(pred.DT,var.DT,by.x='variable',by.y='pc')
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
pred.DT <- merge(pred.DT,control.DT,by.x='variable',by.y='PC')
pred.DT[,variance:=((log(n)-(log(n1) + log(n-n1)))+ log(mfactor)) %>% exp]
## add in control loading
pred.DT[,Z:=(value-control.loading)/sqrt(variance)]
pred.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]

basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,value:=(value-control.loading)]
basis.DT[,Z:=sign(value)]
basis.DT[,p.adj:=1]
basis.DT[,control.loading:=NULL]

sig.traits <- pred.DT[p.adj<0.05,]$trait %>% unique
plot.DT <- rbind(pred.DT[trait %in% sig.traits,.(trait,variable,value=value-control.loading,Z,p.adj)],basis.DT)

library(cowplot)
plot.DT[,variable:=factor(variable,levels=paste('PC',1:11,sep=""))]
plot.DT[p.adj>0.05,Z:=sign(value)]
plot.DT[,trait:=factor(trait,levels=plot.DT$trait %>% unique)]

ggplot(plot.DT,aes(x=trait,y=variable,fill=Z,label=signif(value,digits=2))) + geom_tile() +
geom_text() + scale_fill_gradientn("Z",colours=c('green','white','red')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))





## get significance
