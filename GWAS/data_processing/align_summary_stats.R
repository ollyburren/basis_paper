library(annotSnpStats)
library(data.table)
library(magrittr)
library(rtracklayer)
library(optparse)

TEST <- FALSE

option_list = list(
  make_option(c("-f", "--file"), type="character",default='',
              help="File to process", metavar="character"),
  make_option(c("-o", "--outdir"), type="character",default='',
              help="Output directory", metavar="character")
)

if(TEST){
  args <- list(outdir="~/tmp/",
  file = "/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new/ra_okada.tab")
}else{
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}
print(args)

LD_FILE <- "/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37_ordered.bed"

## load in UK10K reference genotype summaries
uk10 <- readRDS("/home/ob219/rds/hpc-work/DATA/UK10K/UK10K_0.005_MAF.RDS")
setnames(uk10,'CHROM','CHR')
uk10[CHR=='X',CHR:='23']
uk10[,CHR:=as.numeric(CHR)]
uk10 <- uk10[order(CHR,POS),]
uk10m <- uk10[,.(CHR,BP=POS,uk10_A1=REF,uk10_A2=ALT,uk10_A2_AF=AF)]
uk10m[,pid:=paste(CHR,BP,sep=':')]



checkallele<-function(f){
  #DT <- fread(f)[1:1e5,]
  DT <- fread(f)
  M <- merge(uk10m[,.(CHR,BP,uk10_A1,uk10_A2,uk10_A2_AF,pid)],DT,by.x='pid',by.y='pid')
  ## remove duplicates as not sure what to do with these
  dups <- M[duplicated(pid),]$pid
  M <- M[!pid %in% dups,]
  alleles <- data.table(al.x = paste(M$uk10_A1,M$uk10_A2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
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
  M[,g.class:=align.class]
  M[g.class=='comp',c('a1','a2'):=list(uk10_A1,uk10_A2)]
  #M[,flipped:=FALSE]
  sw <- align.class %in% c("rev","revcomp","ambig")
  M[sw,c('a1','a2','or'):=list(uk10_A1,uk10_A2,1/or)]
  ## fix comp
  #M[align.class=='comp',c('a1','a2','g.class'):=list(uk10_A1,uk10_A2,'match')]
  M[align.class=='impossible',c('or','p.value'):=list(NA,NA)]
  #ld.gr<-import.bed(LD_FILE)
  #addLDBlock<-function(DT,ld.gr){
  # dt.gr<-with(DT,GRanges(seqnames=Rle(CHR),ranges=IRanges(start=BP,width=1L),idx=1:nrow(DT)))
  # ol<-as.matrix(findOverlaps(dt.gr,ld.gr))
  # DT[ol[,1],ld:=ol[,2]]
  #}
  #addLDBlock(M,ld.gr)
  M
}

ret<-checkallele(args$file)
out <- ret[order(CHR,BP),.(pid,a1,a2,or,p.value,class=g.class,maf_a2=signif(uk10_A2_AF,digits=3),ld)]

write.table(out,file=file.path(args$outdir,basename(args$file)),row.names=FALSE,quote=FALSE)


if(FALSE){
  data.dir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new'
  files <- list.files(path=data.dir,pattern="*.tab",full.names=TRUE)
  script <- '/home/ob219/git/as_basis/R/Individual_projection/align_summary_stats.R'
  odir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k'
  all.cmds<-sapply(files,function(f){
    cmd <- sprintf("Rscript %s -f %s -o %s",script,f,odir)
  })
  write(all.cmds,file="~/tmp/qsub/align_sum.txt")
}
