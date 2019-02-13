## look at correlation

BB_LU <- list(
  CD = 'crohns.disease',
  CEL = 'malabsorption.coeliac.disease',
  MS = 'multiple.sclerosis',
  RA = 'rheumatoid.arthritis',
  SLE = 'systemic.lupus.erythematosis.sle',
  T1D = 'type.1.diabetes',
  UC = 'ulcerative.colitis',
  asthma = 'asthma'
)

# trial with T1D

library(cupcake)

GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
shrink.DT <- readRDS(SHRINKAGE_FILE)

## get quantiles

quants <- quantile(shrink.DT$ws_emp_shrinkage,prob=seq(0,1,length.out=11))

## decile bins
bins<-split(shrink.DT$pid,cut(shrink.DT$ws_emp_shrinkage,quants))



## load in bb data

traits <- names(BB_LU)

res.DT <- lapply(traits,function(t){
  dat <- get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE,c(t,sprintf("bb_%s",t)))
  dat[trait==t,trait:='t1']
  dat[trait==sprintf("bb_%s",t),trait:='t2']
  dat[,Z:=qnorm(p.value/2,lower.tail=FALSE) * sign(log(or))]
  M <- melt(dat,id.vars=c('pid','trait'),measure.vars='Z') %>% dcast(.,pid~trait+variable)
  lapply(seq_along(bins),function(i){
    bin <- names(bins)[i]
    data.table(trait=t,bin=bin,cor=with(M[pid %in% bins[[i]],],cor(t1_Z,t2_Z)))
  }) %>% rbindlist
}) %>% rbindlist

res.DT[,bin:=factor(bin,levels=res.DT$bin %>% unique)]

ggplot(res.DT,aes(x=bin,y=cor,color=trait,group=trait)) + geom_line()



## whatcauses the uptick on the left is it MAF ?

library(cowplot)
ggplot(res.DT,aes(x=bin,y=cor)) + geom_point()

maf.strat <- lapply(seq_along(bins),function(i){
  DT <- basis.DT[pid %in% bins[[i]],.(pid,maf)] %>% unique
  DT[,bin:=names(bins)[i]]
  DT
}) %>% rbindlist

maf.strat[,bin:=factor(bin,levels=maf.strat$bin %>% unique)]

ggplot(maf.strat,aes(x=bin,y=maf)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


## perhaps LD region size ?

ld.dat <- basis.DT[,.(pid,ld.block)] %>% unique
ld.dat[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
size.dat <- ld.dat[,list(size=max(pos)-min(pos)),by=ld.block]

ld.strat <- lapply(seq_along(bins),function(i){
  DT <- basis.DT[pid %in% bins[[i]],.(pid,ld.block)] %>% unique
  DT[,bin:=names(bins)[i]]
  DT
}) %>% rbindlist

ld.strat <- merge(ld.strat,size.dat,by='ld.block')

ld.strat[,bin:=factor(bin,levels=maf.strat$bin %>% unique)]
ggplot(ld.strat,aes(x=bin,y=size)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
