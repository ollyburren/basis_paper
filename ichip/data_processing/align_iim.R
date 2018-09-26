## processing IIM data so we can predict posterior odds ratios
SNP_MANIFEST_FILE <- '/home/ob219/share/as_basis/ichip/snp_manifest/ichip_summary_stats.tab'
OUT_DIR <- '/home/ob219/share/as_basis/ichip/individual_gt/iim'

library(annotSnpStats)
X <- annot.read.plink('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/IChip/iim_subgroups/iim_subgroups_total')
Y <- X
snps.DT <- snps(Y) %>% data.table
snps.DT[,pid:=paste(chromosome,position,sep=':')]

man.DT <- fread(SNP_MANIFEST_FILE)

keep.idx <- which(snps.DT$pid %in% man.DT$pid)
Y <- Y[,keep.idx]
snps.DT <- snps(Y) %>% data.table
snps.DT[,pid:=paste(chromosome,position,sep=':')]

dup.idx <- which(duplicated(snps.DT$pid))
Y <- Y[,-dup.idx]
snps.DT <- snps.DT[-dup.idx,]

snps.DF <- as.data.frame(snps.DT)
rownames(snps.DF) <- snps.DF$pid
snps(Y) <- snps.DF
colnames(Y) <- snps.DF$pid
rownames(Y) <- samples(Y)$member
rownames(samples(Y)) <- samples(Y)$member

## load in control data for alignment

(load('/home/ob219/share/as_basis/ichip/ctrl_gt/aligned13_filt_b37.RData'))
X<-X[,colnames(X) %in% snps.DF$pid]
X.bar <- X[,match(colnames(X),colnames(Y))]
## getting memory error try downsample of X.bar
X.bar.sample <- X.bar[sample.int(nrow(samples(X.bar)),2000),]
Z<-align.alleles(Y,X.bar.sample)
snps.DT <- snps(Z) %>% data.table
snps.DT[,pid:=paste(chromosome,position,sep=':')]

rownames(Z) <- samples(Z)$member
rownames(samples(Z)) <- samples(Z)$member

## split by chromosome
by.chr <- split(1:nrow(snps.DT),paste0('chr',snps.DT$chromosome))

for(chr in names(by.chr)){
  message(chr)
  idx <- by.chr[[chr]]
  Y <- Z[,idx]
  fname <- file.path(OUT_DIR,sprintf("%s.RDS",chr))
  saveRDS(Y,file=fname)
}

## process clinical subgroups

library(xlsx)
iim.samples <- read.xlsx('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/IChip/iim_subgroups/iim_clincial_subgroups.xls',1) %>% data.table
saveRDS(iim.samples,file="/home/ob219/share/as_basis/ichip/ind_sample_info/iim_samples_info.RDS")
