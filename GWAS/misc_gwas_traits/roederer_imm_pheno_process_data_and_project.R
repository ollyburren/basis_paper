## code to process data from Roederer et al [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4393780/]
library(rtracklayer)
library(annotSnpStats)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/roederer/'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'
DATA_DIR <- '/home/ob219/rds/hpc-work/roederer/twinr-ftp.kcl.ac.uk/ImmuneCellScience/2-GWASResults'
N_FILE <- '/home/ob219/rds/hpc-work/roederer/twinr-ftp.kcl.ac.uk/ImmuneCellScience/2-GWASResults/sample_size.txt'


all.files <- list.files(path=DATA_DIR,pattern="txt.gz",full.names=TRUE)

for(file in all.files){
  ## use this for traitname
  tra <- basename(file) %>% gsub(".txt.gz","",.) %>% sprintf("roederer_%s",.)
  outfile <- file.path(OUT_DIR,sprintf("%s.RDS",tra))
  if(file.exists(outfile)){
    sprintf("Already processed %s skipping",file) %>% message
    next
  }
  sprintf("Reading in %s ",file) %>% message
  dat <- fread(sprintf("zcat %s",file),select=c('name','chromosome','position','allele1','allele2','effallele','build','n','beta','sebeta','p'))
  message("Done")
  ## note at least in this example effect allele is a2 !
  ## first remap to b37

  dat[,pid:=paste(chromosome,position,sep=':')]
  dat <- dat[!pid %in% dat[duplicated(pid),],]
  dat[,id:=1:.N]
  setnames(dat,make.names(names(dat)))
  dat.36.gr <- with(dat,GRanges(seqnames=Rle(paste0('chr',chromosome)),ranges=IRanges(start=position,width=1L),id=id))
  c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
  dat.37.gr<-unlist(liftOver(dat.36.gr,c))
  DT.37 <- data.table(id=dat.37.gr$id,position.37=start(dat.37.gr))
  dat <- merge(dat,DT.37,by.x='id',by.y='id',all.x=TRUE)
  ## 279 dat don't match after coord conversion
  dat <- dat[!is.na(position.37),]
  dat[,pid.37:=paste(chromosome,position.37,sep=":")]

  man.DT <- fread(SNP_MANIFEST)
  M <- merge(dat[,.(pid=pid.37,a1=allele1,a2=allele2,beta,p.value=p,n)],man.DT,by='pid')
  ss <- sprintf("%s %.1f %.1f",tra,M$n %>% mean(.,rm.na=TRUE),M$n %>% max(.,rm.na=TRUE))
  write(ss,file=N_FILE,append=TRUE)
  ## still missing SNPs not sure worth imputing as even with that 1900 missing
  alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
  align.class <- rep('match',nrow(alleles))
  idx<-which(alleles$al.x!=alleles$al.y)
  x.alleles <- alleles[idx,]$al.x
  names(x.alleles)<-alleles[idx,]$pid
  y.alleles <-  alleles[idx,]$al.y
  names(y.alleles)<-names(x.alleles)
  align.class[idx] <- g.class(x.alleles,y.alleles)
  print(table(align.class))
  alleles[,g.class:=align.class]
  idx<-which(alleles$g.class=='impossible')
  if(length(idx) >0){
    M <- M[-idx,]
    alleles <- alleles[-idx,]
  }
  ## check direction which is the effect allele ? It appears that a2 is the effect allele
  M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
  M <- M[g.class %in% c('rev','revcomp'),beta:=beta*-1]
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * tmp$beta
  tmp[is.na(metric),metric:=0]
  tmp[,trait:=tra]
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",tra))
  saveRDS(tmp[,.(pid,beta,p.value,ws_emp_shrinkage)],file=pfile)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  saveRDS(all.proj,file=outfile)
}
