## uk10k data as our reference data set

## process PBC

## first convert to annotSnpStats
library(annotSnpStats)
library(data.table)
library(magrittr)
library(rtracklayer)


uk10 <- readRDS("/home/ob219/rds/hpc-work/DATA/UK10K/UK10K_0.005_MAF.RDS")
setnames(uk10,'CHROM','CHR')
uk10[CHR=='X',CHR:='23']
uk10[,CHR:=as.numeric(CHR)]
uk10 <- uk10[order(CHR,POS),]
uk10m <- uk10[,.(CHR,BP=POS,uk10_A1=REF,uk10_A2=ALT,uk10_A2_AF=AF)]
uk10m[,pid:=paste(CHR,BP,sep=':')]




DATA.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/pbc-anderson/pbc2018'

# code to do allele switch without requiring y to be annot snp stats
## x is annot snp stats
## y is DT of a1 a2 af_wrt2 - that matches snps in X
## shameless rip off of Chris' code on annotSnpStats

switcheroo <- function(x,y,mafdiff=0.1,do.plot=FALSE){
  asw <-  getFromNamespace("asw", "annotSnpStats")
  x.alleles <- apply(x@snps[,alleles(x)],1,paste,collapse="/")
  y.alleles <- apply(y[,.(a1,a2)],1,paste,collapse="/")
  names(y.alleles)<-names(x.alleles)
  message("These are the allele codes as currently defined, before any switching:")
  print(tt <- as.matrix(table(x.alleles, y.alleles)))
  ## genotype classes
  sw.class <- g.class(x.alleles,y.alleles)
  any.comp <- any(sw.class %in% c("comp","revcomp"))
  any.ambig <- any(sw.class=="ambig")
  sw <- sw.class %in% c("rev","revcomp")

  if(length(wh <- which(sw.class=="impossible"))) {
    message(length(wh)," pairwise impossible allele labels found. These genotypes will be set to missing.")
    sw[wh] <- NA
  }
  sw[sw.class=="impossible"] <- NA

  if(any.comp & any.ambig) { # there are reverse complements in the distinguishable cases
    ind <- which(sw.class=="ambig")
    message(length(ind)," SNPs have alleles not completely resolvable without strand information, confirming guess by checking allele freqs.")
    x.cs <- col.summary(x[,ind])
    rdiff <- x.cs[,"RAF"] - y[ind,]$raf
    sw2 <- ifelse(abs(x.cs[,"RAF"] - y[ind,]$raf) < abs(1 - x.cs[,"RAF"] - y[ind,]$raf), FALSE, TRUE)
    too.close <- abs(x.cs[,"MAF"]-0.5)<mafdiff
   # if(any(too.close)) {
   #   can.match <- sw.class %in% c("comp","nochange","rev","revcomp")
   #   xsw <- switch.alleles(x[,-ind],which(sw[-ind]))
   #   ysw <- y[,-ind]
   #   ## step through too.close SNPs checking signed correlation
   #   message("using signed correlation for ",sum(too.close)," SNPs too close to 50% MAF")
   #   ldx <- ld(xsw,x[,ind[too.close],drop=FALSE], stats="R")
   #   ldy <- ld(ysw,y[,ind[too.close],drop=FALSE], stats="R")
   #   ldx[abs(ldx)<0.04] <- NA ## drop uncorrelated - have no information
   #   ldy[abs(ldy)<0.04] <- NA ## drop uncorrelated - have no information
   #   cor.sw <- sapply(1:ncol(ldx), function(j) cor(ldx[,j], ldy[,j], use="pair"))
   #   cor.sw[ abs(cor.sw)<0.8 ] <- NA # NA unless correlation is pretty strong
   #   sw2[too.close] <- cor.sw < 0
   #   too.close <- too.close[is.na(cor.sw)]
   # }
    message(sum(is.na(sw2))," SNPs not resolvable (MAF too close to 0.5).")
    sw[ind] <- sw2
 }

 if(!any.comp & any.ambig) { # there are no reverse complements in distinguishable cases
   ind <- which(sw.class=="ambig")
   message(length(ind)," SNPs have alleles not completely resolvable without strand information,\nbut there is no evidence of strand switches amongst SNPs which are resolvable.\nAssuming fixed strand.")
   ind <- which(sw.class=="ambig")
   sw2 <- x.alleles[ind]==g.rev(y.alleles[ind])
   sw[ind] <- sw2
 }
 ## do switch
 if(any(is.na(sw))){
   x@.Data[,is.na(sw)] <- as.raw("00")
 }


 if(length(wh <- which(sw))) {
     xsw <- asw(x@.Data,wh)
     x@.Data=xsw
     ##x <- switch.alleles(x, wh)
     x@snps[wh, alleles(x) ] <- y[wh, .(a1,a2)]
 }

 if(do.plot) {
   x.cs <- col.summary(x)
   plot(x.cs[,"RAF"], y$raf,main="RAF after switching",xlab="x",ylab="y",pch="+")
   abline(0,1,col="red",lty=1)
 }
 return(x)
}


save_by_chr <- function(X,out_dir='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/pbc-anderson/pbc2018/annotSnpStats',prefix){
  Y <- copy(X)
  snps <- cbind(snps(Y) %>% data.table,col.summary(Y) %>% data.table)
  snps[,posid:=1:.N]
  snps <- merge(snps,uk10m,by.x=c('chromosome','position'),by.y=c('CHR','BP'))
  snps[,pid:=paste(chromosome,position,sep=':')]
  ## if there are duplicates this indicate non-binary alleles which we remove
  dup.pid <- snps[duplicated(pid),]$pid
  snps <- snps[!pid %in% dup.pid,]
  Y <- Y[,snps$posid]
  y<-snps[,.(a1=uk10_A1,a2=uk10_A2,raf=uk10_A2_AF)]
  Y<-switcheroo(Y,y,do.plot=TRUE)
  lapply(split(1:nrow(snps(Y)),snps(Y)$chr),function(i){
    G<-Y[,i]
    fname=sprintf("%s_%d.RData",prefix,unique(snps(G)$chr))
    save(G,file=file.path(out_dir,fname))
  })
  samples <- cbind(samples(Y) %>% data.table,row.summary(Y) %>% data.table)
  snps <- cbind(snps(Y) %>% data.table,col.summary(Y) %>% data.table)
  snps <- merge(snps,uk10m,by.x=c('chromosome','position'),by.y=c('CHR','BP'))
  list(snps=snps,samples=samples)
}

## as this is exome back bone there are a tonne of rare SNPs that are genotyped that have a very low MAF - only take forward the ones in
## the refefence panel (i.e MAF>-.5%)

X.ilmn <- annot.read.plink(file.path(DATA.DIR,'coreex_pbcsmg_20180226.illuminus'))
summ.ilmn <- save_by_chr(X.ilmn,prefix='iluminus')
saveRDS(summ.ilmn,file='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/pbc-anderson/pbc2018/annotSnpStats/ilmn.summary.RDS')

X.gc <- annot.read.plink(file.path(DATA.DIR,'coreex_pbcsmg_20180226.gencall.smajor'))
summ.gc <- save_by_chr(X.gc,prefix='gencall')
saveRDS(summ.gc,file='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/pbc-anderson/pbc2018/annotSnpStats/gc.summary.RDS')

## GTex

## now dealt with by switch_alleles_GTEX_q.R




## JDM

DATA.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/imputed-1kg-apr2018'
gt.files <- list.files(path=DATA.DIR,pattern="*.RData",full.names=TRUE)
f<-gt.files[1]

processJDM <- function(f,out.dir='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/jdm-gwas-data/imputed-1kg-apr2018/as_basis'){
  (load(f))
  ## add chromosome to snps
  snps <- snps(A)
  snps$CHR <- gsub("chr[-]([0-9]+)\\.RData","\\1",basename(f)) %>% as.numeric
  ## alter alleles to more cannonical names
  colnames(snps)[4:5] <- c('allele.1','allele.2')
  snps(A) <- snps
  snps <- snps(A) %>% data.table
  ## add chr to snps object
  #snps[,CHR:=gsub("chr[-]([0-9]+)\\.RData","\\1",basename(f)) %>% as.numeric]
  snps[,posid:=1:.N]
  #setnames(snps,c('a0','a1'),c('allele.1','allele.2'))
  alleles(A) <- c('allele.1','allele.2')

  snps <- merge(snps,uk10m,by.x=c('CHR','position'),by.y=c('CHR','BP'))
  snps[,pid:=paste(CHR,position,sep=':')]
  ## if there are duplicates this indicate non-binary alleles which we remove
  dup.pid <- snps[duplicated(pid),]$pid
  if(length(dup.pid)>0)
    snps <- snps[!pid %in% dup.pid,]
  Y <- copy(A)
  Y <- Y[,snps$posid]
  y<-snps[,.(a1=uk10_A1,a2=uk10_A2,raf=uk10_A2_AF)]
  Y<-switcheroo(Y,y,do.plot=TRUE)
  fname <- file.path(out.dir,gsub("[-]","",basename(f)))
  save(Y,file=fname)
  samples <- cbind(samples(Y) %>% data.table,row.summary(Y) %>% data.table)
  snps <- cbind(snps(Y) %>% data.table,col.summary(Y) %>% data.table)
  snps <- merge(snps,uk10m,by.x=c('CHR','position'),by.y=c('CHR','BP'))
  fname <- file.path(out.dir,sprintf("summ_%s",gsub("[-]","",basename(f))))
  fname <- gsub("RData","RDS",fname)
  saveRDS(list(snps=snps,samples=samples),file=fname)
  fname
}

summ.files<-lapply(gt.files,processJDM)

## JIA can't process all together need external script note if does not work
## try filtering on samples so only cases are aligned (this is currently what we have done)

DATA.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/JIA-2017-data'
gt.files <- list.files(path=DATA.DIR,pattern="*.RData",full.names=TRUE)
if(FALSE){
  DATA.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/JIA-2017-data'
  gt.files <- list.files(path=DATA.DIR,pattern="*.RData",full.names=TRUE)
  RSCRIPT <- '/home/ob219/git/as_basis/R/Individual_projection/switch_alleleles_JIA_q.R'
  cmds <- lapply(gt.files,function(f){
    sprintf("Rscript %s -f %s",RSCRIPT,f)
  }) %>% do.call('c',.)
  write(cmds,file="~/tmp/qsub/jia_align.txt")
  ## run using sh script
}


## for pid we have in VCF
#you need to load plink
## module load plink-1.9-gcc-5.4.0-sm3ojoi


PLINK="/home/cew54/localc/bin/plink" # plink binary
BCFTOOLS="/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools" # bcftools binary
OUTDIR="/home/ob219/rds/rds-who1000-wgs10k/analysis/pid/gt_snps_maf05/pid_cases/"
SAMPLE_FILE="/home/ob219/rds/rds-who1000-wgs10k/analysis/pid/gt_snps_maf05/pid_cases/case_samples_AbDef.txt"

files <- list.files(path='/home/ob219/rds/rds-who1000-wgs10k/analysis/pid/gt_snps_maf05/',pattern="*.vcf.gz$",full.names=TRUE)
cmds <- lapply(files,function(f){
  fname <- gsub("^([^_]+)\\_.*","\\1_pid_cases.bcf",basename(f))
  cmd <- sprintf("%s view -S %s -v snps -Ob -o %s %s",BCFTOOLS,SAMPLE_FILE,file.path(OUTDIR,fname),f)
}) %>% do.call('c',.)
write(cmds,file="~/tmp/qsub/pid_case_split.txt")

OUTDIR='/home/ob219/rds/rds-who1000-wgs10k/analysis/pid/gt_snps_maf05/pid_cases/plink'
PLINK='/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/plink-1.9-sm3ojoi5m2lcorbwetkripy5tomhy7bg/bin/plink'
files <- list.files(path='/home/ob219/rds/rds-who1000-wgs10k/analysis/pid/gt_snps_maf05/pid_cases/bcf/',pattern="*.bcf$",full.names=TRUE)
cmds <- lapply(files,function(f){
  fname <- file.path(OUTDIR,gsub("\\.bcf","",basename(f)))
  cmd <- sprintf("%s --bcf %s  --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --set-missing-var-ids @:#[b37]\\$1,\\$2 --make-bed --out %s",PLINK,f,fname)
}) %>% do.call('c',.)
write(cmds,file="~/tmp/qsub/pid_plink.txt")


## process UK10K in LD blocks !!!

ld <- fread("/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37_ordered.bed")[,.(chr=V1,coords=V4,ld.id=1:.N)]
UK10K_DIR <- '/home/ob219/rds/rds-cew54-wallace-share/olly/UK10K/bcf/'
OUT_DIR <- '/home/ob219/rds/hpc-work/DATA/UK10K/by.ld.block/bcf'
BCFTOOLS="/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools" # bcftools binary
for(i in 1:nrow(ld)){
  ## first spit out as a bcf file
  chr <- paste0('chr',ld$chr[i])
  id <- ld$ld.id[i]
  coords <- ld$coords[i]
  bcf_file <- file.path(UK10K_DIR,sprintf("%s.bcf.gz",chr))
  out_file <- file.path(OUT_DIR,sprintf("%d_%s.bcf",id,coords))
  cmd <- sprintf("%s view -r %s -v snps --min-af 0.005:minor --max-alleles 2 --min-alleles 2 -Ob -o %s %s",BCFTOOLS,paste0('chr',coords),out_file,bcf_file)
}


## next take this and create plink files

cmd <- sprintf("%s --bcf %s  --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --set-missing-var-ids @:#[b37]\\$1,\\$2 --make-bed --out %s",PLINK,"/home/ob219/rds/hpc-work/DATA/UK10K/by.ld.block/bcf/1_1:1-888659.bcf","~/tmp/foo")

## bim file is screwed can we recompute from summary data ?


export i=22
export PREFIX=/home/ob219/rds/rds-cew54-wallace-share/olly/UK10K/bcf/
export INDIR=/home/ob219/rds/rds-cew54-wallace-share/Data/reference/UK10K/
export BCFTOOLS=/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools
for i in `seq 1 22`; do
    fout=~/tmp/qsub/douk_chr${i}.sh
    echo $i
    echo "zcat ${INDIR}_EGAZ00001017893_${i}.UK10K_COHORT.REL-2012-06-02.beagle.anno.csq.shapeit.20140306.legend.gz | sed \"s/chr${i}[^ ]*  /chr${i}:/\" | sed 's/ /_/g' | gzip -c > ${PREFIX}chr${i}.legend.gz" >> $fout
    echo "cp ${INDIR}_EGAZ00001017893_${i}.UK10K_COHORT.REL-2012-06-02.beagle.anno.csq.shapeit.20140306.hap.gz ${PREFIX}chr${i}.hap.gz" >> $fout
    echo "cat ${INDIR}_EGAZ00001017893_${i}.UK10K_COHORT.REL-2012-06-02.beagle.anno.csq.shapeit.20140306.sample | tail -n +2 >  ${PREFIX}chr${i}.samples" >> $fout
    echo "$BCFTOOLS convert -o ${PREFIX}chr${i}.bcf -Ob  --haplegendsample2vcf ${PREFIX}chr${i}" >> $fout
    echo "$BCFTOOLS index ${PREFIX}chr${i}.bcf" >> $fout
    echo "rm ${PREFIX}chr${i}.hap.gz ${PREFIX}chr${i}.samples ${PREFIX}chr${i}.legend.gz" >> $fout
    echo "/usr/bin/sh $fout" >> ~tmp/qsub/run_uk10k.txt
done

## now convert into tomohawk output

export INDIR=/home/ob219/rds/rds-cew54-wallace-share/olly/UK10K/bcf/
export OUTDIR=/home/ob219/rds/rds-cew54-wallace-share/olly/UK10K/tmk/
export TOMAHAWK=/home/ob219/git/tomahawk/tomahawk
export i=1
for i in `seq 1 22`; do
  echo "$TOMAHAWK import -i ${INDIR}chr${i}.bcf.gz -o ${OUTDIR}chr${i}"
  $TOMAHAWK calc -pdfi ${OUTDIR}chr${i}.twk -o ${OUTDIR}chr${i} -a 1 -r 0.5 -P 0.1 -t 28 -c 990 -C 1
done

## wrapper for tomohawk so we can obtain r2 for an LD block



## load in plink files align and save as annotsnpStats

processPID <- function(f,out_dir='/home/ob219/rds/rds-who1000-wgs10k/analysis/pid/gt_snps_maf05/pid_cases/annotSnpStats'){
  Y <- annot.read.plink(gsub("\\.bed","",f))
  snps <- cbind(snps(Y) %>% data.table,col.summary(Y) %>% data.table)
  snps[,posid:=1:.N]
  snps <- merge(snps,uk10m,by.x=c('chromosome','position'),by.y=c('CHR','BP'))
  snps[,pid:=paste(chromosome,position,sep=':')]
  ## if there are duplicates this indicate non-binary alleles which we remove
  dup.pid <- snps[duplicated(pid),]$pid
  snps <- snps[!pid %in% dup.pid,]
  Y <- Y[,snps$posid]
  y<-snps[,.(a1=uk10_A1,a2=uk10_A2,raf=uk10_A2_AF)]
  Y<-switcheroo(Y,y,do.plot=TRUE)
  fname <- file.path(out_dir,paste(gsub("\\.bed","",basename(f)),'RData',sep='.'))
  save(Y,file=fname)
  samples <- cbind(samples(Y) %>% data.table,row.summary(Y) %>% data.table)
  snps <- cbind(snps(Y) %>% data.table,col.summary(Y) %>% data.table)
  snps <- merge(snps,uk10m,by.x=c('chromosome','position'),by.y=c('CHR','BP'))
  fname <- gsub("RData","RDS",fname)
  saveRDS(list(snps=snps,samples=samples),file=fname)
  fname
}



files <- list.files(path='/home/ob219/rds/rds-who1000-wgs10k/analysis/pid/gt_snps_maf05/pid_cases/plink',pattern="*.bed",full.names=TRUE)
out <- lapply(files,processPID)


## need to realign summary stats that are used to build the basis



# if(any(too.close)) {
#   can.match <- sw.class %in% c("comp","nochange","rev","revcomp")
#   xsw <- switch.alleles(x[,-ind],which(sw[-ind]))
#   ysw <- y[,-ind]
#   ## step through too.close SNPs checking signed correlation
#   message("using signed correlation for ",sum(too.close)," SNPs too close to 50% MAF")
#   ldx <- ld(xsw,x[,ind[too.close],drop=FALSE], stats="R")
#   ldy <- ld(ysw,y[,ind[too.close],drop=FALSE], stats="R")
#   ldx[abs(ldx)<0.04] <- NA ## drop uncorrelated - have no information
#   ldy[abs(ldy)<0.04] <- NA ## drop uncorrelated - have no information
#   cor.sw <- sapply(1:ncol(ldx), function(j) cor(ldx[,j], ldy[,j], use="pair"))
#   cor.sw[ abs(cor.sw)<0.8 ] <- NA # NA unless correlation is pretty strong
#   sw2[too.close] <- cor.sw < 0
#   too.close <- too.close[is.na(cor.sw)]
# }


library(parallel)
foo<-mclapply(files,checkallele,mc.cores=8)
