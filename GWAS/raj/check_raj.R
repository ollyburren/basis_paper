## testing to see where I have screwed up

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
}) %>% do.call('cbind',.)


filt[[1]] %>% snps

snp.rhs.estimates

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
