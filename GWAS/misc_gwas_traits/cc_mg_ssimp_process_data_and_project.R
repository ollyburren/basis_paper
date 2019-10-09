## NMO  GWAS
library(annotSnpStats)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
DATA.DIR <- '/home/ob219/share/Data/GWAS-summary/NMO_summary_stats_2018/NMO/imputed_summary_stats'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/nmo/'


## NOT DONE YET !!!
n_controls <- 1244
cases <- list(IgPos=132,IgNeg=83,combined=215)
trait <- names(cases)[1]

res <- lapply(names(cases),function(trait){
  message(trait)
  pattern <- sprintf("*.%s_imputed.txt",trait)
  all.files <- list.files(path=DATA.DIR,pattern=pattern,full.names=TRUE)
  nmo.DT <- mclapply(all.files,function(f){
    fread(f,select=c('chr','pos','Allele1','Allele2','maf','z_imp','P.imp','bst.imp','N.imp','r2.pred','source'))
  },mc.cores=8) %>% rbindlist
  #myo.DT <- fread(in.file,select=c('chr','pos','Allele1','Allele2','maf','z_imp','P.imp','bst.imp','N.imp','r2.pred'))
  snps.DT <- fread(SNP_MANIFEST)
  nmo.DT[,pid:=paste(chr,pos,sep=':')]
  nmo.DT <- nmo.DT[!pid %in% nmo.DT[duplicated(pid),],]
  nmo.DT[,id:=1:.N]
  M <- merge(snps.DT,nmo.DT,by='pid',all.x=TRUE)
  n_cases <- cases[[trait]]
  M[,or_imp:=exp(z_imp/sqrt(2 * ((n_cases * n_controls)/(n_cases + n_controls)) * r2.pred * maf * (1-maf)))]
  ## note that the odds ratios from the imputation are wrong and need to be wholesale flipped - sorry !
  M[,or_imp:=1/or_imp]
  old.dat <- readRDS(file.path(SRC_OUT_DIR,sprintf("NMO_%s_source.RDS",gsub("^Ig","IgG",trait))))
  M2 <- merge(M,old.dat[,.(pid,old.or=or,old.p=p.value)],by='pid',all.x=TRUE)
  ## fill in or using original data where available - remember that things are flipped !
  M2[!is.na(old.or),c('or_imp','P.imp'):=list(old.or,old.p)]
  #M[is.na(z_imp),]
  #missing 657
  ## a2 is the effect allele
  #M <- M[!is.na(z_imp),.(pid,ref_a1,ref_a2,a1=Allele1,a2=Allele2,maf,z_imp,P.imp,bst.imp,N.imp,ref_a1.af,r2.pred)]
  #n_cases <- cases[[trait]]
  #M[,or_imp:=exp(z_imp/sqrt(2 * ((n_cases * n_controls)/(n_cases + n_controls)) * r2.pred * maf * (1-maf)))]
  #M <- M[,.(pid,ref_a1,ref_a2,a1,a2,P=P.imp,or=or_imp)]
  #alleles <- data.table(al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
  ## to make quick
  #align.class <- rep('match',nrow(alleles))
  #idx<-which(alleles$al.x!=alleles$al.y)
  #x.alleles <- alleles[idx,]$al.x
  #names(x.alleles)<-alleles[idx,]$pid
  #y.alleles <-  alleles[idx,]$al.y
  #names(y.alleles)<-names(x.alleles)
  #align.class[idx] <- g.class(x.alleles,y.alleles)
  #print(table(align.class))
  #M[,g.class:=align.class]
  #sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M2,pid)
  tmp <- merge(M2,stmp,by='pid',all.y=TRUE)
  tra <- sprintf("NMO_%s_ssimp",trait)
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS", tra))
  saveRDS(tmp[,.(pid,or=or_imp,p.value=P.imp,ws_emp_shrinkage,source)],file=pfile)
  ## think I messed this up so all are aligned to a1 should be a2 to fix flip direction
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or_imp)
  ## note that there are some SNPs where r2.pred is 0 we set to these to zero
  ## as the or is 1 or 0 and the std.err associated is infinite also
  ## where snp is missing make it zero
  #tmp[is.na(metric),or=1]
  tmp[is.na(metric) | !is.finite(metric),metric:=0]
  #tmp[,trait:= 'myositis_myogen']
  tmp[,trait:= tra]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  #saveRDS(all.proj,file='/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis.RDS')
  ofile <- file.path(OUT_DIR,sprintf("NMO_%s_ssimp_0619.RDS", tra))
  saveRDS(all.proj,file=ofile)
  all.proj
})
