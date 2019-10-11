## ARTERIDO PSA  GWAS
library(annotSnpStats)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'
SRC_OUT_DIR <- '~/share/as_basis/GWAS/for_fdr_13_traits_0919/'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/renton_mg/'
DATA.DIR <- '/home/ob219/share/Data/GWAS-summary/MyastheniaGravis_Renton_JAMA_Neurol_2015/ssimp_imputed'



n_controls <- 1977
cases <- list(YoungOnset=235,LateOnset=737,Overall=972)
convert2original <- list(YoungOnset='renton_mg_early_source.RDS',LateOnset='renton_mg_late_source.RDS',Overall='renton_mg_source.RDS')
trait <- names(cases)[1]
snps.DT <- fread(SNP_MANIFEST)
res <- lapply(names(cases),function(trait){
  message(trait)
  mg.DT <- fread(file.path(DATA.DIR,sprintf("MyastheniaGravis_%s_Renton_JAMA_Neurol_2015_ssimp_imputed.tab",trait)),select=c('chr','pos','Allele1','Allele2','maf','z_imp','P.imp','bst.imp','N.imp','r2.pred','source'))
  #myo.DT <- fread(in.file,select=c('chr','pos','Allele1','Allele2','maf','z_imp','P.imp','bst.imp','N.imp','r2.pred'))
  src.file <- convert2original[[trait]] %>% file.path(SRC_OUT_DIR,.)
  orig.DT <- merge(readRDS(src.file),snps.DT,all.x=TRUE)
  missing <- orig.DT[is.na(or),]$pid
  mg.DT[,pid:=paste(chr,pos,sep=':')]
  mg.DT <- mg.DT[!pid %in% mg.DT[duplicated(pid),],]
  mg.DT[,id:=1:.N]
  M <- merge(snps.DT,mg.DT,by='pid',all.x=TRUE)
  M <- M[!is.na(z_imp),.(pid,ref_a1,ref_a2,a1=Allele1,a2=Allele2,maf,z_imp,P.imp,bst.imp,N.imp,ref_a1.af,r2.pred)]
  n_cases <- cases[[trait]]
  M[,or:=exp(z_imp/sqrt(2 * ((n_cases * n_controls)/(n_cases + n_controls)) * r2.pred * maf * (1-maf)))]
  M[,seb:=1/sqrt(2 * ((n_cases * n_controls)/(n_cases + n_controls)) * r2.pred * maf * (1-maf))]
  #M[,seb:=1/sqrt(2 * ((n_cases * n_controls)/(n_cases + n_controls)) * r2.pred * maf * (1-maf))]
  M <- M[pid %in% missing,.(pid,ref_a1,ref_a2,a1,a2,p.value=P.imp,or,seb)]
  #M[,or:=1/or]
  M <- rbind(orig.DT[!pid %in% missing,.(pid,or,p.value)],M[,.(pid,or,seb,p.value)],fill=TRUE)
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  ## need to compute the sebs of the missing data
  tmp[is.na(seb) & !is.na(or),seb:=log(or)/qnorm(p.value/2,lower.tail=FALSE)]
  ## these are impossible to deal with so should be set to missing
  tmp[!is.finite(seb) &!is.na(or),c('or','p.value','seb'):=list(NA,NA,NA)]
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_ssimp_source.RDS", trait))
  saveRDS(tmp[,.(pid,or,seb,p.value,ws_emp_shrinkage)],file=pfile)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  ## note that there are some SNPs where r2.pred is 0 we set to these to zero
  ## as the or is 1 or 0 and the std.err associated is infinite also
  ## where snp is missing make it zero
  #tmp[is.na(metric),or=1]
  tmp[is.na(metric) | !is.finite(metric),metric:=0]
  #tmp[,trait:= 'myositis_myogen']
  tmp[,trait:= trait]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  #saveRDS(all.proj,file='/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis.RDS')
  ofile <- file.path(OUT_DIR,sprintf("%s_ssimp_13_traits_0919.RDS", trait))
  saveRDS(all.proj,file=ofile)
  all.proj
})
