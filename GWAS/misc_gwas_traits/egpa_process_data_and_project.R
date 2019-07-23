## EGPA GWAS
library(annotSnpStats)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'

aav.dir <- '~/share/Data/GWAS/egpa_aav/summary'
aav.files <- list.files(path=aav.dir,pattern='*.gwas',full.names=TRUE)
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/lyons_egpa/'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'

samples <- list()
for(f in aav.files){
  a.DT <- fread(f)
  trait=gsub("autosomes\\.([^\\.]+).*","\\1",basename(f))
  samples[[trait]] = data.table(trait=trait,n0=a.DT$controls_total[1],n1=a.DT$cases_total[1])
  a.DT <- a.DT[,.(pid=paste(CHR,BP,sep=':'),a1=a0,a2=a1,or=OR,p.value=P)]
  man.DT <- fread(SNP_MANIFEST)
  M <- merge(a.DT,man.DT,by='pid')
  alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
  #alleles <- alleles[!duplicated(pid),]
  #alleles <- M[,list(al.x=paste(uk10_A1,uk10_A2,sep='/'),al.y=paste(a1,a2,sep='/')),by='pid']
  ## to make quick
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

  ## check direction which is the effect allele ?
  M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
  M <- M[!duplicated(pid),]
  M <- M[g.class!='match',or:=1/or]
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  ## where snp is missing make it zero
  tmp[is.na(metric),metric:=0]
  tmp[,trait:= trait]
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",trait))
  saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  ofile <- file.path(OUT_DIR,sprintf("projections/%s_0619.RDS",trait))
  saveRDS(all.proj,file=ofile)
}

samples.DT <- rbindlist(samples)
saveRDS(samples.DT,file=file.path(OUT_DIR,'sample.RDS'))
