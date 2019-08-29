## compute variance for projection of ssimp imputed myogen gwas

## this file contains the driver SNP information from Chris it has to be
## paired with the correct shrinkage file and basis !
DRIVER_SNP_FILE <- '~/share/as_basis/sparse-basis/basis-sparse-0.999.RData'
SHRINKAGE_METHOD<-'ws_emp_shrinkage'
## location of genotypes for LD computation filtered by basis SNPs for speed.
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
## basis shrinkage metrics
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
load(DRIVER_SNP_FILE)
## the only SNPs we care about are in use.pca
keep.snps.dt <- strsplit(rownames(rot.pca),':') %>% do.call('rbind',.) %>% data.table
keep.snps.dt <- data.table(pid=rownames(rot.pca))
keep.snps.dt[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]


## here DT is a data table returned from loading ssimp output with fread(.,select=c('chr','pos','z_imp','maf','r2.pred'))
## N0 number of controls
## N1 number of cases
## ref_gt_dir location of genotypes
## quiet whether to show debug information

compute_ssimp_var <- function(DT,N0,N1,ref_gt_dir,quiet=FALSE){
  imp.filt <- merge(DT,keep.snps.dt,by=c('chr','pos'),all.y=TRUE)
  ## standard error of the beta is computed like so
  imp.filt[,seb:=1/sqrt(2 * ((N0*N1)/(N0+N1)) * r2.pred * maf * (1-maf))]
  ## where seb is NA i.e. the SNP is missing from the dataset set its seb to 0 as it cannot contribute to the variance
  imp.filt[is.na(seb),seb:=0]
  rot.pca.dt <- data.table(pid=rownames(rot.pca),rot.pca)
  M <- merge(imp.filt[,.(pid,seb)],data.table(pid=rownames(rot.pca),rot.pca),by='pid')
  ## add in the shrinkage
  shrink.dt <- readRDS(SHRINKAGE_FILE)
  M <- merge(M,shrink.dt[,.(pid,shrink=get(`SHRINKAGE_METHOD`))],by='pid')
  by.chr <- split(M,gsub("^([^:]+):.*","\\1",M$pid))
  all.chr <- lapply(names(by.chr),function(chr){
    if(!quiet)
      message(sprintf("Processing %s",chr))
    ss.file<-file.path(ref_gt_dir,sprintf("%s.RDS",chr))
    #message(ss.file)
    sm <- readRDS(ss.file)
    ## there are sometimes duplicates that we need to remove
    pids <- colnames(sm)
    dup.idx<-which(duplicated(pids))
    if(length(dup.idx)>0){
      if(!quiet)
        message(sprintf("Warning removing %d duplicated SNPs",length(dup.idx)))
      sm <- sm[,-dup.idx]
      pids <- pids[-dup.idx]
    }
    M <- by.chr[[chr]]
    sm.map <- match(M$pid,pids)

    if(any(is.na(sm.map))){
      message("SNPs in manifest that don't have genotypes")
    }
    r <- ld(sm[,sm.map],sm[,sm.map],stats="R")
    pc.cols <- which(grepl("^PC[0-9]+$",names(M)))
    sapply(names(M)[pc.cols],function(pc) (tcrossprod(M[,.(tot=get(`pc`) * shrink * seb)]$tot) * r) %>% sum)
  })
  colSums(do.call('rbind',all.chr))
}

## example

imp.DT <- fread("~/share/Data/GWAS-summary/MYOGEN/imputed_summary_stats/gwas_dmjdmpm_new_ssimp.meta.txt",select=c('chr','pos','z_imp','maf','r2.pred'))
compute_ssimp_var(imp.DT,N0=4724,N1=1711,ref_gt_dir=REF_GT_DIR,quiet=FALSE)
