library(annotSnpStats)
library(GenomicRanges)
library(rtracklayer)
#se,or etc.

LD_FILE <- "/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37_ordered.bed"
#SNP_MANIFEST_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/snp_manifest/ichip_september.tab'
SNP_MANIFEST_FILE <- '/home/ob219/share/as_basis/ichip/snp_manifest/ichip_september.tab'
ICHIP_DIR <- '/home/ob219/share/as_basis/ichip/sum_stats'
GT_DIR <- '/home/ob219/share/as_basis/ichip/ctrl_gt'

## aligned odds ratios
ss <- fread('~/share/Projects/coloccc/aligned13.csv')
## build an overall manifest
traits <- names(ss)[grep("beta",names(ss))] %>% gsub("beta\\.","",.)
  ## genotype data
(load('~/share/Projects/coloccc/aligned13.RData'))
snps.DT <- snps(X) %>% data.table
snps.DT[,pid:=rownames(snps(X))]
## Chris says to use PID as the uid
keep <- which(snps.DT$pid %in% ss$pid)
Y <- X[,keep]
snps.DT <- snps.DT[keep,]
## convert to build 37
new.pid <- merge(snps.DT[,.(pid,order=1:.N)],ss[,.(pid,pid37)],by='pid')[order(order),]$pid37
colnames(Y) <- new.pid
rownames(snps(Y)) <- new.pid
snps.DT[,c('pid36','pid'):=list(pid,new.pid)]
snps.DT[,c('chromosome','position'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
## create manifest by removing MHC and adding LD blocks and allele freq
snps.gr <- with(snps.DT,GRanges(seqnames=Rle(chromosome),ranges=IRanges(start=position,width=1L),pid=pid))
LD_FILE <- "/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37_ordered.bed"
ld.gr<-import.bed(LD_FILE)
ol <- findOverlaps(snps.gr,ld.gr) %>% as.matrix
snps.DT[ol[,1],region:=ol[,2]]
snps.DT[,af.wrt.a2:=col.summary(Y)$RAF]
mhc.gr <- GRanges(seqnames=Rle(6),ranges=IRanges(start=25e6,end=40e6))
ol <- findOverlaps(snps.gr,mhc.gr) %>% as.matrix
man.DT <- snps.DT[,.(pid,ref_a1=allele.1,ref_a2=allele.2,ref_a1.af=1-af.wrt.a2,ld.block=region)]
man.DT <- man.DT[-ol[,1],]
man.DT <- man.DT[!is.na(ld.block),]
keep <- which(rownames(snps(Y)) %in% man.DT$pid)
Y <- Y[,keep]
## keep a copy moved to build 37
X<-Y
save(X,file="/home/ob219/share/as_basis/ichip/ctrl_gt/aligned13_filt_b37.RData")
## check things are OK
#identical(snps(Y) %>% rownames,man.DT$pid)
#[1] TRUE
## filter aligned odds ratio so match the manifest file
## note that further filtering is required such that all traits have data for all snps
ss <- ss[pid37 %in% man.DT$pid,]
## next create summary stat source files and save
for(trait in traits){
  message(sprintf("Processing %s",trait))
  cols <- c('pid37','RAF','mm',paste("beta",trait,sep='.'),paste("se",trait,sep='.'))
  tmp <- ss[,cols,with=FALSE]
  setnames(tmp,c('pid','RAF','alleles','lor','se_lor'))
  tmp <- tmp[!(is.na(lor)|is.na(se_lor)|is.na(RAF)),][,c('a1','a2'):=tstrsplit(alleles,'/')]
  ## add allele 1 and 2
  tmp[,p.value:=pnorm(abs(lor/se_lor),lower.tail=FALSE) * 2 %>% signif(.,digits=3)]
  tmp <- tmp[!is.na(p.value),]
  fname=file.path(ICHIP_DIR,sprintf("%s.tab",trait))
  ## to combat very low pvalues for now just add in a very small number
  sprintf("Warning %d variant(s) have pvalue of 0",sum(tmp$p.value==0)) %>% message
  tmp[p.value==0,p.value:=1e-300]
  write.table(tmp[,.(pid,a1,a2,or=exp(lor) %>% signif(.,digits=3),p.value)],file=fname,quote=FALSE,sep="\t",row.names=FALSE)
  summary(-log10(tmp$p.value)) %>% print
}
# save manifest file
write.table(man.DT,file=SNP_MANIFEST_FILE,row.names=FALSE,quote=FALSE)

## finally dump out genotype information for variance estimation purposes

man.DT[,c('chr','bp'):=tstrsplit(pid,':')]
by.chr <- split(1:nrow(man.DT),man.DT$chr)
Y.sm <- sm(Y)
lapply(seq_along(by.chr),function(i){
  message(sprintf("Processing chr%s",names(by.chr)[i]))
  chr.keep <- by.chr[[i]]
  sm <- sm(Y.sm[,chr.keep])
  fname <- file.path(GT_DIR,sprintf("%s.RDS",names(by.chr)[i]))
  saveRDS(sm,file=fname)
})
