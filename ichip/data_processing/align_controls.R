library(annotSnpStats)
library(rtracklayer)

## processing IIM data so we can predict posterior odds ratios
SNP_MANIFEST_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/snp_manifest/ichip_september.tab'
OUT_DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/individual_gt/ctrl'
Y <- load('/home/ob219/rds/rds-cew54-wallace-share/Projects/coloccc/allele.flip.aligned9.build37.RData') %>% get
Y<-Y[samples(Y)$affected==1,]
snps.DT <- snps(Y) %>% data.table
snps.DT[,pid:=rownames(snps(X))]
snps.DF <- as.data.frame(snps.DT)
rownames(snps.DF) <- snps.DF$pid
snps(Y) <- snps.DF
snps.DT <- snps(Y) %>% data.table
man.DT <- fread(SNP_MANIFEST_FILE)
keep.idx <- which(snps.DT$pid %in% man.DT$pid)
Y <- Y[,keep.idx]
snps.DT <- snps(Y) %>% data.table

## split by chromosome
by.chr <- split(1:nrow(snps.DT),paste0('chr',snps.DT$chromosome))

for(chr in names(by.chr)){
  message(chr)
  idx <- by.chr[[chr]]
  X <- Y[,idx]
  fname <- file.path(OUT_DIR,sprintf("%s.RDS",chr))
  saveRDS(X,file=fname)
}

# library(xlsx)
# jia.samples <- read.xlsx('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/IChip/JIA/JIA_Immunochip_phaseIItotaldataset_passedQC_UKONLY.xlsx',1) %>% data.table
# jia.samples[,NA.:=NULL]
# saveRDS(jia.samples,file="/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/ind_sample_info/jia_samples_info.RDS")
