library(annotSnpStats)
library(rtracklayer)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'
DATA_DIR <- '~/share/Data/GWAS-summary/MYOGEN/imputed_summary_stats/'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/myogen_myositis/'
SRC_OUT_DIR <- '~/share/as_basis/GWAS/for_fdr_13_traits_0919/'

n_controls <- 4724
cases <- list(dm=705,jdm=473,pm=533,dmjdmpm=1711)

trait <- names(cases)[4]
snps.DT <- fread(SNP_MANIFEST)
res <- mclapply(names(cases),function(trait){
  in.file <- sprintf("gwas_%s_new_ssimp.meta.txt",trait) %>% file.path(DATA_DIR,.)
  src.file <- sprintf("%s_myogen_source.RDS",trait) %>% file.path(SRC_OUT_DIR,.)
  if(trait=='dmjdmpm')
    src.file <- file.path(SRC_OUT_DIR,"myositis_myogen_source.RDS")
  myo.DT <- fread(in.file,select=c('chr','pos','Allele1','Allele2','maf','z_imp','P.imp','bst.imp','N.imp','r2.pred'))
  orig.DT <- merge(readRDS(src.file),snps.DT,all.x=TRUE)
  missing <- orig.DT[is.na(or),]$pid
  myo.DT[,pid:=paste(chr,pos,sep=':')]
  myo.DT <- myo.DT[!pid %in% myo.DT[duplicated(pid),],]
  myo.DT[,id:=1:.N]
  M <- merge(snps.DT,myo.DT,by='pid',all.x=TRUE)
  #M[is.na(z_imp),]
  #missing 657
  ## a2 is the effect allele
  M <- M[!is.na(z_imp),.(pid,ref_a1,ref_a2,a1=Allele1,a2=Allele2,maf,z_imp,P.imp,bst.imp,N.imp,ref_a1.af,r2.pred)]
  n_cases <- cases[[trait]]
  M[,or:=exp(z_imp/sqrt(2 * ((n_cases * n_controls)/(n_cases + n_controls)) * r2.pred * maf * (1-maf)))]
  M[,seb:=1/sqrt(2 * ((n_cases * n_controls)/(n_cases + n_controls)) * r2.pred * maf * (1-maf))]
  #M[,seb:=1/sqrt(2 * ((n_cases * n_controls)/(n_cases + n_controls)) * r2.pred * maf * (1-maf))]
  M <- M[pid %in% missing,.(pid,ref_a1,ref_a2,a1,a2,p.value=P.imp,or,seb)]
  alleles <- data.table(al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
  ## to make quick
  align.class <- rep('match',nrow(alleles))
  idx<-which(alleles$al.x!=alleles$al.y)
  x.alleles <- alleles[idx,]$al.x
  names(x.alleles)<-alleles[idx,]$pid
  y.alleles <-  alleles[idx,]$al.y
  names(y.alleles)<-names(x.alleles)
  align.class[idx] <- g.class(x.alleles,y.alleles)
  print(table(align.class))
  M[,g.class:=align.class]
  M <- rbind(orig.DT[!pid %in% missing,.(pid,or,p.value)],M[,.(pid,or,seb,p.value)],fill=TRUE)
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  ## need to compute the sebs of the missing data
  tmp[is.na(seb) & !is.na(or),seb:=log(or)/qnorm(p.value/2,lower.tail=FALSE)]
  ## these are impossible to deal with so should be set to missing
  tmp[!is.finite(seb) &!is.na(or),c('or','p.value','seb'):=list(NA,NA,NA)]
  tra <- basename(in.file) %>% gsub(".meta.txt","",.)
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_ssimp_source.RDS", tra))
  saveRDS(tmp[,.(pid,or,seb,p.value,ws_emp_shrinkage)],file=pfile)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
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
  ofile <- file.path(OUT_DIR,sprintf("%s_ssimp_13_traits_0919.RDS", tra))
  saveRDS(all.proj,file=ofile)
  all.proj
},mc.cores=4)
