## once we have aligned and processed files
## it is much easier to add traits to the basis as
## the fdr source files contain the harmonised beta/or values
## that are required for projection - all we then need to do
## is assign the weights and hey presto.


library(parallel)
#BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'
DATA_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr/'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr_13_traits_0919/'
pc.emp <- readRDS(BASIS_FILE)
man.DT <- fread(SNP_MANIFEST_FILE)
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]

all.filez <- list.files(path=DATA_DIR,pattern='*.RDS',full.names=TRUE)
all.filez <- all.filez[-grep("snps_only_source",all.filez)]

quant <- grep("^roederer",basename(all.filez))
## strictly the cytokine data is not on the or scale but
## i to the exp so the case control routines should work for it

case_control <- mclapply(all.filez[-quant],function(f){
  message(f)
  trait <- basename(f) %>% gsub("\\_source.RDS","",.)
  dat <- readRDS(f)
  dat <- dat[pid %in% man.DT$pid,.(pid,or,p.value)]
  tmp <- merge(dat,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  tmp$trait <- trait
  pfile <- file.path(OUT_DIR,basename(f))
  saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
  tmp[is.na(metric) | ! is.finite(metric),c('metric','trait'):=list(0,trait)]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
},mc.cores=8) %>% do.call('rbind',.)

library(pheatmap)

case_control[-grep("^CK",rownames(case_control)),] %>% pheatmap(.,cluster_cols=FALSE)

f<-all.filez[quant][1]
qua <- mclapply(all.filez[quant],function(f){
  message(f)
  trait <- basename(f) %>% gsub("\\_source.RDS","",.)
  dat <- readRDS(f)
  dat <- dat[pid %in% man.DT$pid,.(pid,beta,p.value)]
  tmp <- merge(dat,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * tmp$beta
  tmp$trait <- trait
  pfile <- file.path(OUT_DIR,basename(f))
  saveRDS(tmp[,.(pid,beta,p.value,ws_emp_shrinkage)],file=pfile)
  tmp[is.na(metric) | ! is.finite(metric),c('metric','trait'):=list(0,trait)]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
},mc.cores=8) %>% do.call('rbind',.)


all.results <- rbind(qua,case_control)

saveRDS(all.results,'/home/ob219/share/as_basis/GWAS/RESULTS/non_uk_bb_for_fdr_13_traits_0919.RDS')
