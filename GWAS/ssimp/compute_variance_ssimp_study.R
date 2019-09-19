library(parallel)

## compute variance for projection of ssimp imputed myogen gwas

## this file contains the driver SNP information from Chris it has to be
## paired with the correct shrinkage file and basis !
#DRIVER_SNP_FILE <- '~/share/as_basis/sparse-basis/basis-sparse-0.999.RData'
## new basis
DRIVER_SNP_FILE <- '~/share/as_basis/sparse-basis/basis-sparse-13-0.999.RData'
SHRINKAGE_METHOD<-'ws_emp_shrinkage'
## location of genotypes for LD computation filtered by basis SNPs for speed.
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
## basis shrinkage metrics
#SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
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
  ## infinite happens if the snp cannot be imputed at all
  imp.filt[is.na(seb) | !is.finite(seb),seb:=0]
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


## run on all myositis subtypes
DATA_DIR <- '~/share/Data/GWAS-summary/MYOGEN/imputed_summary_stats/'
OUT_FILE <- '/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_empirical_variances_13_traits_0919.RDS'
n_controls <- 4724
cases <- list(dm=705,jdm=473,pm=533,dmjdmpm=1711)

variances <- mclapply(names(cases),function(trait){
  in.file <- sprintf("gwas_%s_new_ssimp.meta.txt",trait) %>% file.path(DATA_DIR,.)
  imp.DT <- fread(in.file,select=c('chr','pos','z_imp','maf','r2.pred'))
  compute_ssimp_var(imp.DT,N0=n_controls,N1=cases[[trait]],ref_gt_dir=REF_GT_DIR,quiet=TRUE)
},mc.cores=8)

traits <- names(cases) %>% sprintf("gwas_%s_new_ssimp",.)
variances.dt <- do.call('rbind',variances) %>% data.table(trait=traits,.) %>% melt(.,id.vars='trait')
saveRDS(variances.dt,file=OUT_FILE)

## run on nmo
DATA_DIR <- '~/share/Data/GWAS-summary/estrada_nmo/imputed_summary_stats/'
OUT_FILE <- '/home/ob219/share/as_basis/GWAS/nmo/nmo_empirical_variances_13_traits_0919.RDS'
n_controls <- 1244
cases <- list(IgPos=132,IgNeg=83,combined=215)
variances.nmo <- lapply(names(cases),function(trait){
  message(trait)
  pattern <- sprintf("*.%s_imputed.txt",trait)
  all.files <- list.files(path=DATA_DIR,pattern=pattern,full.names=TRUE)
  imp.DT <- mclapply(all.files,function(f){
    fread(f,select=c('chr','pos','Allele1','Allele2','maf','z_imp','P.imp','bst.imp','N.imp','r2.pred','source'))
  },mc.cores=8) %>% rbindlist
  compute_ssimp_var(imp.DT,N0=n_controls,N1=cases[[trait]],ref_gt_dir=REF_GT_DIR,quiet=TRUE)
})

traits <- names(cases) %>% sprintf("NMO_%s_ssimp",.)
variances.dt <- do.call('rbind',variances.nmo) %>% data.table(trait=traits,.) %>% melt(.,id.vars='trait')
saveRDS(variances.dt,file=OUT_FILE)

## aterido_psa

DATA_DIR <- '/home/ob219/share/Data/GWAS-summary/psa-aterido/ssimp_imputed'
OUT_FILE <- '/home/ob219/share/as_basis/GWAS/psa_aterido/psa_aterido_empirical_variances_13_traits_0919.RDS'
n_controls <- 1454 ## slight cheat as n controls slightly different unlikely to make huge difference
cases <- list(span_psa=744,na_psa=1430)

variances.psa <- mclapply(names(cases),function(trait){
  message(trait)
  in.file <- sprintf("%s_imputed.txt",trait) %>% file.path(DATA_DIR,.)
  imp.DT <- fread(in.file,select=c('chr','pos','z_imp','maf','r2.pred'))
  compute_ssimp_var(imp.DT,N0=n_controls,N1=cases[[trait]],ref_gt_dir=REF_GT_DIR,quiet=TRUE)
},mc.cores=8)

traits <- names(cases) %>% sprintf("%s_ssimp",.)
variances.dt <- do.call('rbind',variances.psa) %>% data.table(trait=traits,.) %>% melt(.,id.vars='trait')
saveRDS(variances.dt,file=OUT_FILE)
