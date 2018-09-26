## processing IIM data so we can predict posterior odds ratios
SNP_MANIFEST_FILE <- '/home/ob219/share/as_basis/ichip/snp_manifest/ichip_summary_stats.tab'
OUT_DIR <- '/home/ob219/share/as_basis/ichip/individual_gt/jia'

library(annotSnpStats)
library(rtracklayer)
X <- annot.read.plink('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/IChip/JIA/plink_dataset_files/Immunochip_JIA_PhaseIInew_QCgp_All_SNPQC_UKonly')
Y <- X
Y<-Y[samples(Y)$affected==2,]
snps.DT <- snps(Y) %>% data.table
snps.DT[,pid:=rownames(snps(X))]
snps.gr <- with(snps.DT,GRanges(seqnames=Rle(paste0('chr',chromosome)),ranges=IRanges(start=position,width=1L),pid=pid))
## liftover - 3 lines ;)
chain <- import.chain("/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain")
snps.gr.37 <- liftOver(snps.gr,chain) %>% unlist
pos.37 <- data.table(pid=snps.gr.37$pid,position.37=start(snps.gr.37))
snps.DT[,'uid':=1:.N]
snps.DT <- merge(snps.DT,pos.37,by.x='pid',by.y='pid',all.x=TRUE)[order(uid),]
snps.DT[,pid:=paste(chromosome,position.37,sep=":")]
snps.DF.37 <- snps.DT[,.(chromosome,snp.name,cM,position=position.37,allele.1,allele.2,pid=paste(chromosome,position.37,sep=':'))] %>% as.data.frame
miss.pos.idx <- which(is.na(snps.DT$position.37))
Y <- Y[,-miss.pos.idx]
snps.DF.37 <- snps.DT[-miss.pos.idx,.(chromosome,snp.name,cM,position=position.37,allele.1,allele.2,pid=paste(chromosome,position.37,sep=':'))] %>% as.data.frame
rownames(snps.DF.37) <- snps.DF.37$pid
snps(Y) <- snps.DF.37
colnames(Y) <- snps.DF.37$pid
snps.DT <- snps(Y) %>% data.table

man.DT <- fread(SNP_MANIFEST_FILE)

keep.idx <- which(snps.DT$pid %in% man.DT$pid)
Y <- Y[,keep.idx]
snps.DT <- snps(Y) %>% data.table
snps.DT[,pid:=paste(chromosome,position,sep=':')]

## load in control data for alignment

(load('/home/ob219/share/as_basis/ichip/ctrl_gt/aligned13_filt_b37.RData'))
X<-X[,colnames(X) %in% snps.DT$pid]
X.bar <- X[,match(colnames(X),colnames(Y))]
## getting memory error try downsample of X.bar
X.bar.sample <- X.bar[sample.int(nrow(samples(X.bar)),2000),]
Z<-align.alleles(Y,X.bar.sample)
snps.DT <- snps(Z) %>% data.table
snps.DT[,pid:=paste(chromosome,position,sep=':')]

## split by chromosome
by.chr <- split(1:nrow(snps.DT),paste0('chr',snps.DT$chromosome))

for(chr in names(by.chr)){
  message(chr)
  idx <- by.chr[[chr]]
  Y <- Z[,idx]
  fname <- file.path(OUT_DIR,sprintf("%s.RDS",chr))
  saveRDS(Y,file=fname)
}

library(xlsx)
jia.samples <- read.xlsx('/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/IChip/JIA/JIA_Immunochip_phaseIItotaldataset_passedQC_UKONLY.xlsx',1) %>% data.table
jia.samples[,NA.:=NULL]
saveRDS(jia.samples,file="/home/ob219/share/as_basis/ichip/ind_sample_info/jia_samples_info.RDS")
