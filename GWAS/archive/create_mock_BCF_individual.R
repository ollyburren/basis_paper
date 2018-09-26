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


out.dir <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/mockVCF'
MANIFEST_FILE <- '/home/ob219/rds/hpc-work/as_basis/support_tab/ind_proj_manifest.tab'
BCFTOOLS <- 'bcftools'

man.DT <- fread(MANIFEST_FILE)

##############
## check PBC #
##############
pbc.dir <- man.DT[trait=='PBC',]$munged
ilmn.dat <- readRDS(file.path(pbc.dir,'ilmn.summary.RDS'))
ilmn.dat <- ilmn.dat$snps
## are there duplicates by position ?
ilmn.dat[,pid:=paste(chromosome,position,sep=':')]
alleles <- data.table(al.x = paste(ilmn.dat$uk10_A1,ilmn.dat$uk10_A2,sep='/'),al.y=paste(ilmn.dat$allele.1,ilmn.dat$allele.2,sep='/'))
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$pid
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
ilmn.dat[,g.class:=align.class]
ggplot(ilmn.dat,aes(x=uk10_A2_AF,y=RAF,col=g.class)) + geom_point() + ggtitle("PBC_ilmn")

gencall.dat <- readRDS(file.path(pbc.dir,'gc.summary.RDS'))
gencall.dat <- gencall.dat$snps
## are there duplicates by position ?
gencall.dat[,pid:=paste(chromosome,position,sep=':')]
alleles <- data.table(al.x = paste(gencall.dat$uk10_A1,gencall.dat$uk10_A2,sep='/'),al.y=paste(gencall.dat$allele.1,gencall.dat$allele.2,sep='/'))
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$pid
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
gencall.dat[,g.class:=align.class]
ggplot(gencall.dat,aes(x=uk10_A2_AF,y=RAF,col=g.class)) + geom_point() + ggtitle("PBC_gencall")

## merge and plot against each other

gc <- gencall.dat[,.(gc.pid=pid,gc.RAF=RAF,gc.class=g.class)]
setkey(gc,gc.pid)
ilmn <- ilmn.dat[,.(ilmn.pid=pid,ilmn.RAF=RAF,ilmn.class=g.class)]
setkey(ilmn,ilmn.pid)
M <- gc[ilmn]
ggplot(M,aes(x=ilmn.RAF,y=gc.RAF,col=gc.class)) + geom_point() + ggtitle("ilmn vs gc")

## USING GENCALL dataset
trait <- 'ind:pbc'
DT <- gencall.dat[!g.class %in% c('impossible','ambig'),]
DT <- DT[chromosome!=23,]
DT[,trait:=trait]
#DT.f <- DT[,.(CHROM,POS,ID=pid,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO=sprintf("CLASS=%s",class),FORMAT='GT',study='1|1')]
DT.f <- DT[,.(CHROM=chromosome,POS=position,ID=pid,REF=uk10_A1,ALT=uk10_A2,QUAL='.',FILTER='.',INFO='.',FORMAT='GT',study='1|1')]
fname <- file.path(out.dir,sprintf("%s.vcf",trait))

write(sprintf(mock_header,trait),file=fname)
write.table(DT.f,file=fname,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
bcfname <- gsub("vcf$","bcf",fname)
cmd <- sprintf("%s convert %s -Ob -o %s",BCFTOOLS,fname,bcfname)
system(cmd)
unlink(fname)
cmd <- sprintf("%s index %s",BCFTOOLS,bcfname)
system(cmd)






## for GTEX create a unified summary of SNPs across all chr
gtex <- list.files(path=man.DT[trait=='GTeX',]$munged,pattern="^summ*",full.names=TRUE)
gtex.dat <- lapply(gtex,function(f){
  obj <- readRDS(f)
  obj$snps
}) %>% rbindlist

gtex.dat[,pid:=paste(chromosome,position,sep=':')]
alleles <- data.table(al.x = paste(gtex.dat$uk10_A1,gtex.dat$uk10_A2,sep='/'),al.y=paste(gtex.dat$allele.1,gtex.dat$allele.2,sep='/'))
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$pid
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
gtex.dat[,g.class:=align.class]
trait <- 'ind:gtex'
DT <- gtex.dat[!g.class %in% c('impossible','ambig'),]
DT <- DT[chromosome!=23,]
DT[,trait:=trait]

ggplot(DT,aes(x=uk10_A2_AF,y=RAF,col=g.class)) + geom_point() + ggtitle("GTEX")

#DT.f <- DT[,.(CHROM,POS,ID=pid,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO=sprintf("CLASS=%s",class),FORMAT='GT',study='1|1')]
DT.f <- DT[,.(CHROM=chromosome,POS=position,ID=pid,REF=uk10_A1,ALT=uk10_A2,QUAL='.',FILTER='.',INFO='.',FORMAT='GT',study='1|1')]
fname <- file.path(out.dir,sprintf("%s.vcf",trait))

write(sprintf(mock_header,trait),file=fname)
write.table(DT.f,file=fname,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
bcfname <- gsub("vcf$","bcf",fname)
cmd <- sprintf("%s convert %s -Ob -o %s",BCFTOOLS,fname,bcfname)
system(cmd)
unlink(fname)
cmd <- sprintf("%s index %s",BCFTOOLS,bcfname)
system(cmd)


## pid

pid <- list.files(path=man.DT[trait=='PID',]$munged,pattern=".RDS$",full.names=TRUE)
pid.dat <- lapply(pid,function(f){
  obj <- readRDS(f)
  obj$snps
}) %>% rbindlist
pid.dat[,pid:=paste(chromosome,position,sep=':')]
alleles <- data.table(al.x = paste(pid.dat$uk10_A1,pid.dat$uk10_A2,sep='/'),al.y=paste(pid.dat$allele.1,pid.dat$allele.2,sep='/'))
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$pid
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
pid.dat[,g.class:=align.class]
trait <- 'ind:pid'
DT <- pid.dat[!g.class %in% c('impossible','ambig'),]
DT <- DT[chromosome!=23,]
DT[,trait:=trait]

ggplot(DT,aes(x=uk10_A2_AF,y=RAF,col=g.class)) + geom_point() + ggtitle("pid")

#DT.f <- DT[,.(CHROM,POS,ID=pid,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO=sprintf("CLASS=%s",class),FORMAT='GT',study='1|1')]
DT.f <- DT[,.(CHROM=chromosome,POS=position,ID=pid,REF=uk10_A1,ALT=uk10_A2,QUAL='.',FILTER='.',INFO='.',FORMAT='GT',study='1|1')]
fname <- file.path(out.dir,sprintf("%s.vcf",trait))

write(sprintf(mock_header,trait),file=fname)
write.table(DT.f,file=fname,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
bcfname <- gsub("vcf$","bcf",fname)
cmd <- sprintf("%s convert %s -Ob -o %s",BCFTOOLS,fname,bcfname)
system(cmd)
unlink(fname)
cmd <- sprintf("%s index %s",BCFTOOLS,bcfname)
system(cmd)

## JDM
jdm <- list.files(path=man.DT[trait=='JDM',]$munged,pattern=".RDS$",full.names=TRUE)
jdm.dat <- lapply(jdm,function(f){
  obj <- readRDS(f)
  obj$snps
}) %>% rbindlist
jdm.dat[,jdm:=paste(CHR,position,sep=':')]
alleles <- data.table(al.x = paste(jdm.dat$uk10_A1,jdm.dat$uk10_A2,sep='/'),al.y=paste(jdm.dat$allele.1,jdm.dat$allele.2,sep='/'))
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$jdm
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
jdm.dat[,g.class:=align.class]
trait <- 'ind:jdm'
DT <- jdm.dat[!g.class %in% c('impossible','ambig'),]
DT <- DT[CHR!=23,]
DT[,trait:=trait]

ggplot(DT,aes(x=uk10_A2_AF,y=RAF,col=g.class)) + geom_point() + ggtitle("jdm")

#DT.f <- DT[,.(CHROM,POS,ID=jdm,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO=sprintf("CLASS=%s",class),FORMAT='GT',study='1|1')]
DT.f <- DT[,.(CHROM=CHR,POS=position,ID=jdm,REF=uk10_A1,ALT=uk10_A2,QUAL='.',FILTER='.',INFO='.',FORMAT='GT',study='1|1')]
fname <- file.path(out.dir,sprintf("%s.vcf",trait))

write(sprintf(mock_header,trait),file=fname)
write.table(DT.f,file=fname,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
bcfname <- gsub("vcf$","bcf",fname)
cmd <- sprintf("%s convert %s -Ob -o %s",BCFTOOLS,fname,bcfname)
system(cmd)
unlink(fname)
cmd <- sprintf("%s index %s",BCFTOOLS,bcfname)
system(cmd)

jia <- list.files(path=man.DT[trait=='JIA',]$munged,pattern=".RDS$",full.names=TRUE)
jia.dat <- lapply(jia,function(f){
  obj <- readRDS(f)
  obj$snps
}) %>% rbindlist
jia.dat[,jia:=paste(chromosome,position,sep=':')]
alleles <- data.table(al.x = paste(jia.dat$uk10_A1,jia.dat$uk10_A2,sep='/'),al.y=paste(jia.dat$allele.1,jia.dat$allele.2,sep='/'))
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$jia
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
jia.dat[,g.class:=align.class]
trait <- 'ind:jia'
DT <- jia.dat[!g.class %in% c('impossible','ambig'),]
DT <- DT[chromosome!=23,]
DT[,trait:=trait]

ggplot(DT,aes(x=uk10_A2_AF,y=RAF,col=g.class)) + geom_point() + ggtitle("jia")

#DT.f <- DT[,.(CHROM,POS,ID=jia,REF=a1,ALT=a2,QUAL='.',FILTER='.',INFO=sprintf("CLASS=%s",class),FORMAT='GT',study='1|1')]
DT.f <- DT[,.(CHROM=chromosome,POS=position,ID=jia,REF=uk10_A1,ALT=uk10_A2,QUAL='.',FILTER='.',INFO='.',FORMAT='GT',study='1|1')]
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

  ## merge all into a giant bcf file !
  files <- list.files(path='/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/mockVCF',pattern="*.bcf$",full.names=TRUE)
  write(files,file='~/tmp/as_basis_merged.txt')
  ## hand edit to remove generated files !!!
  system('bcftools merge -l ~/tmp/as_basis_merged.txt -Ob -o /home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/mockVCF/all_including_ind_studies_missing.bcf -m snps')
  ## this add tags
  #export BCFTOOLS_PLUGINS=/usr/local/Cluster-Apps/bcftools/1.2/plugins
  #bcftools plugin missing2ref  -Ob -o all_including_ind_studies.bcf all_including_ind_studies_missing.bcf
  #bcftools plugin fill-AN-AC -Ob -o  all_including_ind_studies_tagged.bcf all_including_ind_studies.bcf
  ## to get a list of SNPs that are covered by all studies we can do

}

## how to use

if(FALSE){
  library(data.table)
  library(rtracklayer)
  ## load in manifest
  m_file<-'/home/ob219/git/as_basis/manifest/as_manifest_mar_2018.csv'
  m.DT <- fread(m_file)
  ## get basis traits
  basis.traits <- m.DT[basis_trait==1,]
  #tfile<-tempfile()
  traits <- paste(c(basis.traits$disease,'ind:jia'),sep=',',collapse=',')
  #write(traits,file=tfile)
  tot_traits <- (nrow(basis.traits) + 1)*2
  BCF_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/mockVCF/all_including_ind_studies_tagged.bcf'
  #cmd <- sprintf("bcftools view -c %d -s %s -Ov -o %s %s ",tot_traits,traits,ofile,BCF_FILE)
  cmd <- sprintf("bcftools view -c %d -s %s -Ov -G %s ",tot_traits,traits,BCF_FILE)

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
  write.table(man.DT,file='/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_uk10k_jia.tab',row.names=FALSE,quote=FALSE)
}
