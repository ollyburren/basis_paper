library(data.table)
library(magrittr)
library(rtracklayer)

## using a study manifest file we create an overall snp manifest that encompasses the SNPs that
## present in all studies.
## Individual BCF's are created by create_mock_BCF(_individual).R these are then merged with bcftools
## to create a merged bcfile. This is subsequently manipulated to add in missing genotypes and carry
## allele counting. Commands for all these steps can be found in create_mock_BCF etc.

MANIFEST_FILE <- '/home/ob219/git/as_basis/manifest/as_manifest_july.tsv'
BCF_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/mockVCF/latest.bcf'
SNP_MANIFEST_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab'
LD_FILE <- "/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37_ordered.bed"
UK10K_FILE <- '/home/ob219/rds/hpc-work/DATA/UK10K/UK10K_0.005_MAF.RDS'

addLDBlock<-function(DT,ld.gr){
  dt.gr<-with(DT,GRanges(seqnames=Rle(CHR),ranges=IRanges(start=BP,width=1L),idx=1:nrow(DT)))
  ol<-as.matrix(findOverlaps(dt.gr,ld.gr))
  DT[ol[,1],ld:=ol[,2]]
}


## tested with bcftools 1.6
m.DT <- fread(MANIFEST_FILE)[include=='Y',]
traits <- paste(m.DT$disease,sep=',',collapse=',')
tot_traits <- nrow(m.DT)*2
cmd <- sprintf("bcftools view -c %d -s %s -Ov -G %s ",tot_traits,traits,BCF_FILE)
message(cmd)
DT <- fread(cmd)


## load in UK10K reference genotype summaries
uk10 <- readRDS(UK10K_FILE)
setnames(uk10,'CHROM','CHR')
uk10[CHR=='X',CHR:='23']
uk10[,CHR:=as.numeric(CHR)]
uk10 <- uk10[order(CHR,POS),]
uk10m <- uk10[,.(CHR,BP=POS,uk10_A1=REF,uk10_A2=ALT,uk10_A2_AF=AF)]
uk10m[,pid:=paste(CHR,BP,sep=':')]
man.DT <- uk10m[pid %in% DT$ID,]
## add ld information
ld.gr<-import.bed(LD_FILE)
addLDBlock(man.DT,ld.gr)
man.DT <- man.DT[,.(pid,ref_a1=uk10_A1,ref_a2=uk10_A2,ref_a1.af=1-uk10_A2_AF,ld.block=ld)]
write.table(man.DT,file=SNP_MANIFEST_FILE,row.names=FALSE,quote=FALSE)

## now run the filter software using SNP_MANIFEST_FILE !
