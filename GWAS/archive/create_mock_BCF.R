library(data.table)
library(parallel)
## build dossier of all SNPs across all traits

mock_header <- "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##contig=<ID=10,length=2147483647>
##contig=<ID=11,length=2147483647>
##contig=<ID=12,length=2147483647>
##contig=<ID=13,length=2147483647>
##contig=<ID=14,length=2147483647>
##contig=<ID=15,length=2147483647>
##contig=<ID=16,length=2147483647>
##contig=<ID=17,length=2147483647>
##contig=<ID=18,length=2147483647>
##contig=<ID=19,length=2147483647>
##contig=<ID=1,length=2147483647>
##contig=<ID=20,length=2147483647>
##contig=<ID=21,length=2147483647>
##contig=<ID=22,length=2147483647>
##contig=<ID=2,length=2147483647>
##contig=<ID=3,length=2147483647>
##contig=<ID=4,length=2147483647>
##contig=<ID=5,length=2147483647>
##contig=<ID=6,length=2147483647>
##contig=<ID=7,length=2147483647>
##contig=<ID=8,length=2147483647>
##contig=<ID=9,length=2147483647>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s"


in.dir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/'
out.dir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/mockVCF'
files <- list.files(path=in.dir,pattern='*.tab',full.names=TRUE)
BCFTOOLS <- '/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools'

library(optparse)

TEST <- TRUE

option_list = list(
  make_option(c("-f", "--file"), type="character",default='',
              help="ld block to process", metavar="character")
)

if(TEST){
  args <- list(file='/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/sun_KIR2DS2.10428.1.3.tab')
}else{
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}

DT <- fread(args$file)[class!='impossible',]
trait <- gsub(".tab","",basename(args$file))
DT[,c('CHROM','POS'):=tstrsplit(pid,':',fixed=TRUE)]
#DT.f <- DT[,.(CHROM,POS,ID=pid,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO=sprintf("CLASS=%s",class),FORMAT='GT',study='1|1')]
DT.f <- DT[,.(CHROM,POS,ID=pid,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO='.',FORMAT='GT',study='1|1')]
fname <- file.path(out.dir,sprintf("%s.vcf",trait))

write(sprintf(mock_header,trait),file=fname)
write.table(DT.f,file=fname,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
bcfname <- gsub("vcf$","bcf",fname)
cmd <- sprintf("%s convert %s -Ob -o %s",BCFTOOLS,fname,bcfname)
system(cmd)
unlink(fname)
cmd <- sprintf("%s index %s",BCFTOOLS,bcfname)
system(cmd)

if(FALSE){
  library(magrittr)
  Rscript <- '/home/ob219/git/as_basis/R/Individual_projection/create_mock_BCF.R'
  files <- list.files(path=in.dir,pattern="*.tab$",full.names=TRUE)
  all.cmds <- lapply(files,function(f){
    cmd <- sprintf("Rscript %s -f %s",Rscript,f)
  }) %>% do.call('c',.)
  write(all.cmds,file="~/tmp/qsub/create_mock_bcf.txt")
  ## merge all into a giant bcf file !
  files <- list.files(path=out.dir,pattern="*.bcf$",full.names=TRUE)
  write(files,file='~/tmp/as_basis_merged.txt')
  system('bcftools merge -l ~/tmp/as_basis_merged.txt -Ob -o /home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/mockVCF/all_studies.bcf -m snps')
  ## this add tags
  # export BCFTOOLS_PLUGINS=/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/bcftools-1.6-6p3lqyearbsrree33gn7blgmuxjiybc5/libexec/bcftools/
  #bcftools plugin missing2ref  -Ob -o all_studies_missing.bcf all_studies.bcf
  #bcftools plugin fill-AN-AC -Ob -o all_studies_tagged.bcf all_studies.bcf
  ## to get a list of SNPs that are covered by all studies we can do

}

## how to use

if(FALSE){
  ## Want selection of SNPs that occur in all studies
  library(data.table)
  library(rtracklayer)
  #BCF_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/mockVCF/all_studies_tagged.bcf'
  BCF_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/mockVCF/all_studies.bcf'
  TOTAL_STUDY_COUNT <- 238
  cmd <- sprintf("bcftools view -G -c %d  -Ov %s",TOTAL_STUDY_COUNT,BCF_FILE)
  message(cmd)
  DT <- fread(cmd)
  LD_FILE <- "/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37_ordered.bed"

  ## load in UK10K reference genotype summaries
  uk10 <- readRDS("/home/ob219/rds/hpc-work/DATA/UK10K/UK10K_0.005_MAF.RDS")
  setnames(uk10,'CHROM','CHR')
  uk10[CHR=='X',CHR:='23']
  uk10[,CHR:=as.numeric(CHR)]
  uk10 <- uk10[order(CHR,POS),]
  uk10m <- uk10[,.(CHR,BP=POS,uk10_A1=REF,uk10_A2=ALT,uk10_A2_AF=AF)]
  uk10m[,pid:=paste(CHR,BP,sep=':')]
  man.DT <- uk10m[pid %in% DT$ID,]
  ## add ld information
  ld.gr<-import.bed(LD_FILE)
  addLDBlock<-function(DT,ld.gr){
    dt.gr<-with(DT,GRanges(seqnames=Rle(CHR),ranges=IRanges(start=BP,width=1L),idx=1:nrow(DT)))
    ol<-as.matrix(findOverlaps(dt.gr,ld.gr))
    DT[ol[,1],ld:=ol[,2]]
  }
  addLDBlock(man.DT,ld.gr)
  man.DT <- man.DT[,.(pid,ref_a1=uk10_A1,ref_a2=uk10_A2,ref_a1.af=1-uk10_A2_AF,ld.block=ld)]
  write.table(man.DT,file='/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_uk10k.tab',row.names=FALSE,quote=FALSE)
}
