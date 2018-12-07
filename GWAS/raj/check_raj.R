## testing to see where I have screwed up

## step one take a region of the genome 50MB either side of PTPN22 and select all the genes that
## have expression data
library(annotSnpStats)
fs <- list.files(path="/home/ob219/share/as_basis/GWAS/individual_data/filtered_gt/raj/cd4",pattern="^chr[0-9]+.RDS",full.names=TRUE)
all.gt<-lapply(fs,readRDS)

## bind manually
## super slow !!!
# snpBind <- function(x,y){
#   snps <- rbind(snps(x),snps(y))
#   gt <- cbind(x@.Data,y@.Data) %>% sm
#   samples <- samples(x)
#   new("aSnpMatrix",
#      .Data=gt,
#      snps=snps,
#      samples=samples,
#      alleles=c("allele.1","allele.2"),
#      phenotype="affected")
# }
# merged.gt <- do.call('snpBind',all.gt)


## much faster ?? Why ??
gt <- lapply(all.gt,function(x) x@.Data) %>% do.call('cbind',.) %>% sm
snps <- lapply(all.gt,snps) %>% do.call('rbind',.)
samples <- samples(all.gt[[1]])

merged.gt <- new("aSnpMatrix",
   .Data=gt,
   snps=snps,
   samples=samples,
   alleles=c("allele.1","allele.2"),
   phenotype="affected")

## fit the linear model for these using aligned data so that we have some indicative summary statistics where
## we know that there is signal.

## i created a list of genes in a 50MB region around PTPN22 that we will use as a trial to see if summary and
## individual level projections are the same

#genes <- readRDS("~/tmp/ptpn22_gene_neighbourhood.RDS")
genes <- scan("~/tmp/IL6R.txt","character")
(load("/home/ob219/share/Projects/twas/raj-cd4-expression.RData"))
texpr <- t(expr)
translate <- readRDS('/home/ob219/share/Projects/twas/model_output/trans_raj-cd4.rds')
idxy<-match(rownames(texpr),names(translate))
rownames(texpr) <- translate[idxy]
(load("/home/ob219/share/Projects/twas/raj-cd4-probes.RData"))
probes <- data.table(probes)
probes2 <- probes[,.(ID,ensg=as.character(mrna_assignment) %>% gsub(".*(ENSG[0-9]+).*","\\1",.))]
probes2 <- probes2[ensg %in% genes,]

## so we need to do 737 linear regressions - actually after filtering 707

texpr.f <- texpr[,colnames(texpr) %in% probes2$ID]
#check that PTPN22 is present otherwise have an issue
#(colnames(texpr) == probes2[ensg=='ENSG00000134242',]$ID) %>% sum

if(FALSE){
## looking at CD4 PTPN22 does not seem to be a strong eqtl according to the paper - try instead
## the list of genes that they find on chromosome one - this allows me to double check my regressions
## with their list also as a secondary check
(load("/home/ob219/share/Projects/twas/raj-cd4-expression.RData"))
texpr <- t(expr)
translate <- readRDS('/home/ob219/share/Projects/twas/model_output/trans_raj-cd4.rds')
idxy<-match(rownames(texpr),names(translate))
rownames(texpr) <- translate[idxy]

## this is from their supp table
pap <- fread("~/tmp/tableS4_eu_cd4T_cis_fdr05.tsv")
texpr.f <- texpr[,colnames(texpr) %in% (pap[GENE_CHR==1,]$PROBESET_ID %>% unique)]

}


## how long does one gene take to compute over all the variants ?


samp <- samples(merged.gt)
idxz <- match(rownames(samp),rownames(texpr.f))
colnames(texpr.f) <- colnames(texpr.f) %>% paste("P",.,sep='_')
samp <- cbind(samp,texpr.f[idxz,])
prob <- 'P_7922309'


all.probes <- colnames(texpr.f)
library(parallel)
all.lm <- mclapply(all.probes,function(prob){
  message(prob)
  res <- snp.rhs.estimates(sprintf("%s~sex",prob) %>% formula,family="gaussian",data=samp,snp.data=sm(merged.gt))
  all.beta <- sapply(res,'[[','beta')
  all.vbeta <- sapply(res,'[[','Var.beta')
  res.DT <- data.table(beta=all.beta,vbeta=all.vbeta,probe=prob,variant=names(all.beta),Z=all.beta/sqrt(all.vbeta))
  res.DT[,variant:=gsub("([^\\.]+)\\..*","\\1",variant)]
  res.DT[,p:=pnorm(abs(Z),lower.tail=FALSE) * 2]
  res.DT
},mc.cores=8)

## next we have to project onto the basis for each probe and thus get the loadings

library(cupcake)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
CAL_FILE <- '/home/ob219/rds/hpc-work/as_basis/support//por_2500_2.0_0.01.RDS'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'

pc.emp <- readRDS(BASIS_FILE)
shrink.DT <- readRDS(SHRINKAGE_FILE)

## which genotypes are missing

missing <- shrink.DT[!pid %in% snps(merged.gt)$snp.name,]
missing.DT <- data.table(pid=missing$pid,trait='DUMMY',or=1)

## want data structure with pid,trait(probe),value,or
sample.proj.DT <- lapply(all.lm,function(x){
  probe <- x$probe[1]
  missing.DT[,trait:=probe]
  DT <- rbind(missing.DT,x[,.(pid=snps(merged.gt)$snp.name,trait=probe,or=exp(beta))])
}) %>% rbindlist
setkey(sample.proj.DT,pid)
mat.emp <- create_ds_matrix(sample.proj.DT,shrink.DT,SHRINKAGE_METHOD)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)

all.proj.DT <- data.table(probe=rownames(all.proj),all.proj)
all.proj.m <- melt(all.proj.DT,id.var='probe')

## compute the mean projection across a component

meany <- all.proj.m[,list(mval=mean(value)),by=variable]
pc.emp.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
ctrl.DT <- pc.emp.DT[trait=='control',]

## next compute empirical Z scores

pc.emp.DT[,Z:=(value-mean(value))/sd(value),by=variable]





data.path <- '/home/ob219/share/as_basis/GWAS/raj/cd4'
saveRDS(all.proj,file.path(data.path,'summ_chr1_proj_IL6.RDS'))

## actually cannot compare directly we load in the regression
## coefficients from regression of gene expression onto individual
## pojections

all.lms <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd4/regression_models.RDS")
names(all.lms)<-paste('PC',1:11,sep="")
for.reg <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd4/cd4.RDS")
pc.t <- lapply(seq_along(names(all.lms)),function(i){
  pc<-names(all.lms)[i]
  all.p <- sapply(all.lms[[i]],function(x){
    summary(x)$coefficient["value",3]
  })
  DT <- data.table(PC=pc,probe=names(for.reg)[grep('^P\\_',names(for.reg))],t.stat=all.p)
})

pc.t <- rbindlist(pc.t)
setnames(all.proj.m,'variable','PC')
M <- merge(pc.t,all.proj.m,by=c('probe','PC'))

M[,PC:=factor(PC,levels=paste('PC',1:11,sep=''))]

ggplot(M,aes(x=value,y=t.stat)) + geom_point() + facet_wrap(~PC) + xlab("summary stats proj") + ylab("individual stats t stat")



## load in individual data
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'
dat <- readRDS("/home/ob219/share/as_basis/GWAS/individual_data/individual_proj/raj_cd4_projection.RDS")
t.DT <- data.table(ind=rownames(dat),dat)
res <- melt(t.DT,id.vars='ind')
pc.emp <- readRDS(BASIS_FILE)
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
tmp <- basis.DT[trait=='control',] %>% t
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))
res <- merge(res,ctrl.DT,by='variable')
res[,delta:=value-control.loading]
res[,variable:=factor(variable,levels=paste('PC',1:11,sep=''))]

## we can then project these onto the basis and compare these projections with those that we get from using individual genotype
## data

library(annotSnpStats)

(load("/home/ob219/share/Projects/twas/raj-cd4-expression.RData"))
pap <- fread("~/tmp/tableS4_eu_cd4T_cis_fdr05.tsv")

library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
snps<-SNPlocs.Hsapiens.dbSNP144.GRCh37
lu<-snpsById(snps,pap$SNP,ifnotfound='drop')
lu<-data.table(as.data.frame(lu))

lu[,pid:=sub("ch","",seqnames) %>% paste(.,pos,sep=':')]
lu <- lu[,.(RefSNP_id,pid)]
M <- merge(pap,lu,by.x='SNP',by.y='RefSNP_id')

## load in genotypes that we have realigned

fs <- list.files(path="/home/ob219/share/as_basis/GWAS/individual_data/filtered_gt/raj/cd4",pattern="^chr[0-9]+.RDS",full.names=TRUE)
all.gt<-lapply(fs,readRDS)

filt <- lapply(all.gt,function(as){
  sdt <- snps(as) %>% data.table
  idx <- which(sdt$snp.name %in% M$pid)
  as[,idx]
}) %>% do.call('rbind2',.)


filt[[1]] %>% snps


texpr <- t(expr)
translate <- readRDS('/home/ob219/share/Projects/twas/model_output/trans_raj-cd4.rds')
idxy<-match(rownames(texpr),names(translate))
rownames(texpr) <- translate[idxy]



samp <- samples(test)
idxz <- match(rownames(samp),rownames(texpr))
colnames(texpr) <- colnames(texpr) %>% paste("P",.,sep='_')
samp <- cbind(samp,texpr[idxz,])

## first take a look at known probe hits on the chromosome

pro <- M[SNP %in% snps(test)$ID,]$PROBESET_ID %>% unique %>% paste("P",.,sep='_')

## for each probeset only test relevant SNPs

bp <- M[SNP %in% snps(test)$ID,.(SNP,PROBESET_ID)]
bp <- split(bp$SNP,paste('P',bp$PROBESET_ID,sep='_'))
all.res <- lapply(seq_along(bp),function(i){
  prob <- names(bp)[i]
  snps <- bp[[i]]
  snpid <- snps(test)$ID %in% snps
  res <- snp.rhs.estimates(sprintf("%s~sex",prob) %>% formula,family="gaussian",data=samp,snp.data=sm(test),sets=snpid)
  foo<-lapply(res,function(x){
    tdt <- data.table(x %>% t)
    tdt[,snp:=x$beta %>% names]
  }) %>% rbindlist
  foo[,Z:=as.numeric(beta)/sqrt(as.numeric(Var.beta))]
}) %>% rbindlist

M[,prob:=paste('P',PROBESET_ID,sep='_')]
all.res[,Y.var:=unlist(Y.var)]

fm <- merge(M,all.res,by.x=c('SNP','prob'),by.y=c('snp','Y.var'))

fdt <- readRDS("~/tmp/raj_chr1_flip.RDS")
fm<-merge(fm,fdt,by.x='SNP',by.y='snp')
setnames(fm,'P-VALUE','p.value')
library(cowplot)
ggplot(fm,aes(x=qnorm(p.value/2,lower.tail=FALSE) * sign(RHO),y=Z,color=sclass)) + geom_point()


## alternative is to load

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'

pc.emp <- readRDS(BASIS_FILE)

## pc3 max rotation

pc3 <- pc.emp$rotation[,"PC3"]

## what could cause this ? Perhaps the effect sizes are so small and local such that
## even with weighting they disappear into the background ?

## there is no signal in bulk cd4 - you need to stimulate them or something ?
