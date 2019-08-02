## script to work out whether fdr results contain any novel findings

## need latest version of R for this - this did not work as could not install gwascat
# export R_LIBS_USER=/home/ob219/R/x86_64-pc-linux-gnu-library/3.6
# module load r-3.6.0-gcc-5.4.0-bzuuksv


basis.sig <- fread("/home/ob219/share/as_basis/GWAS/RESULTS/basis-sig-drivers.csv")[p.value>5e-8,.(pid=gsub("^P","",pid),trait,RefSNP_id)] %>% unique
basis.sig <- basis.sig[,list(traits=paste(trait,sep=',',collapse=',')),by=c('RefSNP_id','pid')]
basis.sig[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
basis.sig <- basis.sig[order(chr,pos),]
basis.sig[,id:=1:.N]

## code to generate cut down bcf files for the lookup

INTERVAL <- 1e6
OUTDIR <- '/home/ob219/share/as_basis/GWAS/cfdr/lookup'
options('scipen'=999)
for(chrs in basis.sig$chr){
  #chrname <- paste('chr',chrs,sep='')
  fname <- file.path(OUTDIR,sprintf("chr%s.txt",chrs))
  dt <- basis.sig[chr==chrs,][order(chr,as.numeric(pos)),]
  write.table(dt[,.(chr=chrs,start=pos-INTERVAL,end=pos + INTERVAL)],file=fname,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

}
options('scipen'=0)

#UK10KDIR <- '/home/ob219/share/Data/reference/UK10K/BCF/'
UK10KDIR <- '/home/ob219/rds/hpc-work/DATA/1kgenome/VCF/EUR/by.chr.phase3/bcf'
BCFOUTDIR <- '/home/ob219/share/as_basis/GWAS/cfdr/lookup/1gp'
for(chrname in paste('chr',basis.sig$chr %>% unique,sep='')){
  message(chrname)
  ukf <- file.path(UK10KDIR,sprintf("ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz",chrname))
  rfile <- file.path(OUTDIR,sprintf("%s.txt",chrname))
  ofile <- file.path(BCFOUTDIR,sprintf("%s.bcf",chrname))
  #cmd <- sprintf("module load bcftools-1.6-gcc-5.4.0-6p3lqye;`which bcftools` view -R %s -Ob %s > %s",rfile,ukf,ofile)
  ## could not get the module stuff to work with the queue so used my own compiled version
  cmd <- sprintf("/home/ob219/bin/bcftools-1.4/bcftools view -R %s -Ob %s > %s",rfile,ukf,ofile)
  write(cmd,file="~/tmp/qstuff/filter.txt",append=TRUE)
  #system(cmd)
  ## run these outside of R so as to preserve session from watchdog
}


library(snpStats)
library(rtracklayer)
bcf2snpmatrix <- function(pid,quiet=FALSE,interval=INTERVAL,min.r2=0.1){
  message(pid)
  ## convert pid into a region based on interval
  chr <- sub("^([^:]+):.*","\\1",pid); start <- sub("^[^:]+:(.*)","\\1",pid) %>% as.numeric
  region <- sprintf('%s:%d-%d',chr,pmax(start-interval,1),start+interval)
  bcf.file <- file.path(BCFOUTDIR,sprintf("chr%s.bcf",chr))
  header_cmd <- sprintf("bcftools view -h %s",bcf.file)
  if(!quiet)
    message(header_cmd)
  my.pipe<-pipe(header_cmd)
  header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
  close(my.pipe)
  cnames<-unlist(strsplit(header,"\t"))
  bcftools_cmd<-sprintf("/home/ob219/bin/bcftools-1.4/bcftools view -q 0.01:minor -r %s -O v %s | grep -v '^#'",region,bcf.file)
  if(!quiet)
    message(bcftools_cmd)
  #tryCatch(as.data.frame(fread(cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE)),error=function(e){print(sprintf("Error=%s fread CMD=%s",e,cmd));return(NA)})
  tmp<-tryCatch(fread(bcftools_cmd),error=function(e){print(sprintf("Error=%s fread CMD=%s",e,bcftools_cmd));return(data.table(pid=pid,error=TRUE))})
  #print(tmp)
  #message(grepl("error",names(tmp)))
  if(any(grepl("error",names(tmp))))
    return(tmp)
  setnames(tmp,cnames)
  gt<-tmp[,10:ncol(tmp),with=FALSE]


  if(nrow(gt)==0)
    return(NA)
  info<-tmp[,1:9,with=FALSE]
  setnames(info,'#CHROM','CHROM')
  if(!quiet)
    message("Creating snpMatrix obj")
  sm<-apply(gt,1,function(x) sub("0\\|0","1",x))
  sm<-apply(sm,1,function(x) sub("(0\\|1)|(1\\|0)","2",x))
  sm<-apply(sm,1,function(x) sub("1\\|1","3",x))
  ## set anything else to a missing value
  sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\|[0-9]","0",x))))
  info[,pid:=paste(CHROM,POS,sep=':') %>% gsub("^chr","",.)]
  colnames(sm)<-info$pid
  rownames(sm)<-colnames(gt)
  ## compute r2 for selected pid with everything else
  sm <- new("SnpMatrix", sm)
  idx <- which(info$pid==pid)
  r2 <- ld(sm[,idx],sm,stats="R.squared") %>% t
  dt <- data.table(qpid=pid,pid=rownames(r2),r2=r2[,1])[r2>min.r2,]
  dt
}

library(parallel)
# works although could get watchdogged !
ld.snps <- mclapply(basis.sig$pid,bcf2snpmatrix,mc.cores=8)


all.ld.snps <- rbindlist(ld.snps)
all.ld.snps[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
all.ld.snps[,id:=1:.N]
## move these onto build 38 so we can see how they match GWAS catalog
all.ld.snps.37.gr <- with(all.ld.snps,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=pos,width=1L),id=id))
c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg19ToHg38.over.chain') ## e.g. hg19ToHg18.over.chain
all.ld.snps.38.gr<-unlist(liftOver(all.ld.snps.37.gr ,c))
DT.38 <- data.table(id=all.ld.snps.38.gr$id,position.38=start(all.ld.snps.38.gr))
all.ld.snps <- merge(all.ld.snps,DT.38,by.x='id',by.y='id',all.x=TRUE)
b38.ld.gr <- with(all.ld.snps[!is.na(position.38),],GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position.38,width=1L),qpid=qpid,r2=r2))


#basis.gr <- with(all.ld.snps,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position.38-1e6,end=position.38+1e6),id=id,traits=traits))


gwas.cat <- fread("~/tmp/gwas_catalog_v1.0-associations_e96_r2019-07-12.tsv")
#gwas.cat <- gwas.cat[`P-VALUE`<5e-8,]
setnames(gwas.cat,make.names(names(gwas.cat)))
gwas.cat <- gwas.cat[grepl("^ankylosing spondylitis$|^latent autoimmune diabetes$|psoriatic|iga nephropathy",`DISEASE.TRAIT`,ignore.case=TRUE),]
gwas.cat <- gwas.cat[!is.na(as.numeric(CHR_POS)),]
gwas.cat[,CHR_POS:=as.numeric(CHR_POS)]
gwas.cat[,id:=1:.N]

as.gr <- with(gwas.cat[DISEASE.TRAIT=='Ankylosing spondylitis',],GRanges(seqnames=Rle(paste0('chr',CHR_ID)),ranges=IRanges(start=CHR_POS,width=1L),gwas.id=id,p=P.VALUE))
psa.gr <- with(gwas.cat[DISEASE.TRAIT=='Psoriatic arthritis',],GRanges(seqnames=Rle(paste0('chr',CHR_ID)),ranges=IRanges(start=CHR_POS,width=1L),gwas.id=id,p=P.VALUE))
iga.gr <- with(gwas.cat[DISEASE.TRAIT=='IgA nephropathy',],GRanges(seqnames=Rle(paste0('chr',CHR_ID)),ranges=IRanges(start=CHR_POS,width=1L),gwas.id=id,p=P.VALUE))
lada.gr <- with(gwas.cat[DISEASE.TRAIT=='Latent autoimmune diabetes',],GRanges(seqnames=Rle(paste0('chr',CHR_ID)),ranges=IRanges(start=CHR_POS,width=1L),gwas.id=id,p=P.VALUE))

as.pids <- basis.sig[grep("ank_spond",traits),]$pid
psa.pids <- basis.sig[grep("bowes_psa",traits),]$pid
iga.pids <- basis.sig[grep("IgA_nephropathy",traits),]$pid
lada.pids <- basis.sig[grep("cousminer_lada",traits),]$pid

as <- mergeByOverlaps(b38.ld.gr[b38.ld.gr$qpid %in% as.pids,],as.gr)
novel.as <- basis.sig[pid %in% as.pids[!as.pids %in% as$qpid],]
psa <- mergeByOverlaps(b38.ld.gr[b38.ld.gr$qpid %in% psa.pids,],psa.gr)
novel.psa <- basis.sig[pid %in% psa.pids[!psa.pids %in% psa$qpid],]
iga <- mergeByOverlaps(b38.ld.gr[b38.ld.gr$qpid %in% iga.pids,],iga.gr)
novel.iga <- basis.sig[pid %in% iga.pids[!iga.pids %in% iga$qpid],]
lada <- mergeByOverlaps(b38.ld.gr[b38.ld.gr$qpid %in% lada.pids,],lada.gr)
novel.lada <- basis.sig[pid %in% lada.pids[!lada.pids %in% lada$qpid],]

## looks ok so recreate original table using pids and trait filters

basis.verb <- fread("/home/ob219/share/as_basis/GWAS/RESULTS/basis-sig-drivers.csv")
as.verb <- basis.verb[pid %in% paste0('P',novel.as$pid) & trait=='ank_spond',]
psa.verb <- basis.verb[pid %in% paste0('P',novel.psa$pid) & trait=='bowes_psa',]
iga.verb <- basis.verb[pid %in% paste0('P',novel.iga$pid) & trait=='IgA_nephropathy',]
lada.verb <- basis.verb[pid %in% paste0('P',novel.lada$pid) & trait=='cousminer_lada',]

novel <- rbindlist(list(as.verb,psa.verb,iga.verb,lada.verb))

library(xlsx)
write.xlsx(novel, file="/home/ob219/tmp/basis-sig-drivers-novel.xlsx", sheetName="Sheet1",col.names=TRUE, row.names=FALSE, append=FALSE, showNA=TRUE)
