library(data.table)
library(magrittr)
library(annotSnpStats)
library(cowplot)
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


## PROCESS INDIVIDUAL DATA FIRST
## NOTE different formats

ind.out.dir <- '/home/ob219/share/as_basis/ichip/BCF/ind_stats'
#MANIFEST_FILE <- '/home/ob219/rds/hpc-work/as_basis/support_tab/ind_proj_manifest.tab'
BCFTOOLS <- 'bcftools'

test.DT <- fread('/home/ob219/share/as_basis/ichip/snp_manifest/ichip_summary_stats.tab')[,.(pid,a1=ref_a1,a2=ref_a2)]

## assumes a dir of annotSnpStats objects
createBCFInd <- function(ind_dir,out_dir,trait,keep.ambig=FALSE){
  ind.files <- list.files(path=ind_dir,pattern="*.RDS",full.names=TRUE)
  DT <- lapply(ind.files,function(f){
    X <- readRDS(f)
    snps(X) %>% data.table
  }) %>% rbindlist
  DT.ind <- DT[,.(pid,ind_a1=allele.1,ind_a2=allele.2)]
  tmp.DT <- merge(DT.ind,test.DT,by.x='pid',by.y='pid')
  y.alleles <- tmp.DT[,list(al=paste(ind_a1,ind_a2,sep="/"))]$al
  x.alleles <- tmp.DT[,list(al=paste(a1,a2,sep="/"))]$al
  names(y.alleles)<-names(x.alleles)
  message("These are the allele codes as currently defined, before any switching:")
  print(tt <- as.matrix(table(x.alleles, y.alleles)))
  tmp.DT[,g.class:=g.class(x.alleles,y.alleles)]
  #trait <- 'ind:ind'
  if(keep.ambig){
    DT <- tmp.DT[!g.class!='impossible',]
  }else{
    DT <- tmp.DT[!g.class %in% c('impossible','ambig'),]
  }

  DT[,c('chromosome','position'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
  DT <- DT[chromosome!=23,]
  DT[,trait:=trait]
  #DT.f <- DT[,.(CHROM,POS,ID=pid,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO=sprintf("CLASS=%s",class),FORMAT='GT',study='1|1')]
  DT.f <- DT[,.(CHROM=chromosome,POS=position,ID=pid,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO='.',FORMAT='GT',study='1|1')]
  DT.f <- DT.f[order(CHROM,POS),]
  fname <- file.path(out_dir,sprintf("%s.vcf",trait))
  write(sprintf(mock_header,trait),file=fname)
  write.table(DT.f,file=fname,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  bcfname <- gsub("vcf$","bcf",fname)
  cmd <- sprintf("%s convert %s -Ob -o %s",BCFTOOLS,fname,bcfname)
  system(cmd)
  unlink(fname)
  cmd <- sprintf("%s index %s",BCFTOOLS,bcfname)
  system(cmd)
}

createBCFInd('/home/ob219/share/as_basis/ichip/individual_gt/iim',ind.out.dir,'ind:imm')
createBCFInd('/home/ob219/share/as_basis/ichip/individual_gt/jia',ind.out.dir,'ind:jia')
# gtfiles <- list.files(path='/home/ob219/share/as_basis/ichip/individual_gt/gtex/archive',pattern='*.RData',full.names=TRUE)
# for(f in gtfiles){
#   load(f)
#   fname <- file.path('/home/ob219/share/as_basis/ichip/individual_gt/gtex',gsub("RData","RDS",basename(f)))
#   saveRDS(Y,file=fname)
#   message(sprintf("Saved %s",fname))
# }
createBCFInd('/home/ob219/share/as_basis/ichip/individual_gt/gtex',ind.out.dir,'ind:gtex')


## summary statistics

createBCFSum <- function(sum_file,out_dir,trait,keep.ambig=FALSE){
  DT <- fread(sum_file)
  DT.sum <- DT[,.(pid,sum_a1=a1,sum_a2=a2)]
  tmp.DT <- merge(DT.sum,test.DT,by.x='pid',by.y='pid')
  y.alleles <- tmp.DT[,list(al=paste(sum_a1,sum_a2,sep="/"))]$al
  x.alleles <- tmp.DT[,list(al=paste(a1,a2,sep="/"))]$al
  names(y.alleles)<-names(x.alleles)
  message("These are the allele codes as currently defined, before any switching:")
  print(tt <- as.matrix(table(x.alleles, y.alleles)))
  tmp.DT[,g.class:=g.class(x.alleles,y.alleles)]
  #trait <- 'ind:ind'
  if(keep.ambig){
    DT <- tmp.DT[!g.class!='impossible',]
  }else{
    DT <- tmp.DT[!g.class %in% c('impossible','ambig'),]
  }
  DT[,c('chromosome','position'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
  DT <- DT[chromosome!=23,]
  DT[,trait:=trait]
  #DT.f <- DT[,.(CHROM,POS,ID=pid,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO=sprintf("CLASS=%s",class),FORMAT='GT',study='1|1')]
  DT.f <- DT[,.(CHROM=chromosome,POS=position,ID=pid,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO='.',FORMAT='GT',study='1|1')]
  DT.f <- DT.f[order(CHROM,POS),]
  fname <- file.path(out_dir,sprintf("%s.vcf",trait))
  write(sprintf(mock_header,trait),file=fname)
  write.table(DT.f,file=fname,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  bcfname <- gsub("vcf$","bcf",fname)
  cmd <- sprintf("%s convert %s -Ob -o %s",BCFTOOLS,fname,bcfname)
  system(cmd)
  unlink(fname)
  cmd <- sprintf("%s index %s",BCFTOOLS,bcfname)
  system(cmd)
}


in.dir <- '/home/ob219/share/as_basis/ichip/sum_stats/'
sum.files <- list.files(path=in.dir,pattern='*.tab',full.names=TRUE)
sum.out.dir <- '/home/ob219/share/as_basis/ichip/BCF/sum_stats'

for(f in sum.files){
  trait <- basename(f) %>% gsub(".tab","",.)
  message(sprintf("Processing %s",trait))
  createBCFSum(f,sum.out.dir,trait)
}




if(FALSE){

  ## merge all into a giant bcf file !
  files.summary <- list.files(path=sum.out.dir,pattern="*.bcf$",full.names=TRUE)
  files.ind <- list.files(path=ind.out.dir,pattern="*.bcf$",full.names=TRUE)
  files <- c(files.summary,files.ind)
  write(files,file='~/tmp/as_basis_merged.txt')
  ## hand edit to remove generated files !!!
  system('bcftools merge -l ~/tmp/as_basis_merged.txt -Ob -o /home/ob219/share/as_basis/ichip/BCF/all.bcf -m snps')
  system('export BCFTOOLS_PLUGINS=/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/bcftools-1.6-6p3lqyearbsrree33gn7blgmuxjiybc5/libexec/bcftools/')
  system('bcftools plugin missing2ref  -Ob -o /home/ob219/share/as_basis/ichip/BCF/all_missing.bcf /home/ob219/share/as_basis/ichip/BCF/all.bcf')
  system('bcftools plugin fill-AN-AC -Ob -o   /home/ob219/share/as_basis/ichip/BCF/all_tagged.bcf  /home/ob219/share/as_basis/ichip/BCF/all_missing.bcf')
  ## this add tags
  # export BCFTOOLS_PLUGINS=/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/bcftools-1.6-6p3lqyearbsrree33gn7blgmuxjiybc5/libexec/bcftools/
  #bcftools plugin missing2ref  -Ob -o all_missing.bcf all.bcf
  #bcftools plugin fill-AN-AC -Ob -o   all_tagged.bcf  all_missing.bcf
  ## to get a list of SNPs that are covered by all studies we can do

}

## how to use

if(FALSE){
  library(data.table)
  library(rtracklayer)
  ## load in manifest
  #m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
  #m.DT <- fread(m_file)
  ## get basis traits
  #basis.traits <- m.DT[basis_trait==1,]
  #tfile<-tempfile()
  #traits <- paste(c(basis.traits$disease,'ind:jia'),sep=',',collapse=',')
  #write(traits,file=tfile)
  tot_traits <- length(files) * 2
  BCF_FILE <- '/home/ob219/share/as_basis/ichip/BCF/all_tagged.bcf'
  #cmd <- sprintf("bcftools view -c %d -s %s -Ov -G %s ",tot_traits,traits,BCF_FILE)
  cmd <- sprintf("bcftools view -c %d -Ov -G %s ",tot_traits,BCF_FILE)

  message(cmd)
  DT <- fread(cmd)
  # build manifest - this file comes from 'process_basis_traits.R in this dir'
  # for the time stick to using old LD block designations (i.e. not the ones for ichip)
  #snps <- fread('/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/snp_manifest/ichip_september.tab')
  man.DT <- fread('/home/ob219/share/as_basis/ichip/snp_manifest/ichip_summary_stats.tab')
  man.DT <- man.DT[pid %in% DT$ID,]
  write.table(man.DT,file='/home/ob219/share/as_basis/ichip/snp_manifest/ichip_all.tab',row.names=FALSE,quote=FALSE)
}
