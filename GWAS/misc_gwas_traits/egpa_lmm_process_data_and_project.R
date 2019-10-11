## EGPA GWAS
library(annotSnpStats)

SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'
SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'

aav.dir <- '/home/ob219/share/Data/GWAS-summary/lyons_EGPA/egpa-summstats/'
aav.files <- list.files(path=aav.dir,pattern='*.txt',full.names=TRUE)
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/lyons_egpa_lmm/'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr_13_traits_0919/'
n0 <- 6688
sc <- list(all.egpa.vs.controls.txt=534,anca.negative.egpa.vs.controls.txt=358,mpo.anca.positive.egpa.vs.controls.txt=159)
convertORscale <- function(x,cp) x/(cp * (1-cp))
## HERE I ASSUME THAT ALLELE1 is the scored ALLELE and thus OR are wrt to this
for(f in aav.files){
  a.DT <- fread(f)
  trait=gsub("\\.vs.controls.txt","",basename(f))
  message(trait)
  #samples[[trait]] = data.table(trait=trait,n0=a.DT$controls_total[1],n1=a.DT$cases_total[1])
  prop <- sc[[basename(f)]]/n0
  a.DT <- a.DT[,.(pid=paste(CHR,BP,sep=':'),a1=ALLELE1,a2=ALLELE0,or=convertORscale(BETA,prop) %>% exp,seb=convertORscale(SE,prop),p.value=P)]
  a.DT[,p.value:=pnorm(abs(log(or)/seb),lower.tail=FALSE) * 2]
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
  M <- M[g.class=='match',or:=1/or]
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  ## where snp is missing make it zero
  tmp[is.na(metric),metric:=0]
  tmp[,trait:= trait]
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_lmm_source.RDS",trait))
  saveRDS(tmp[,.(pid,or,seb,p.value,ws_emp_shrinkage)],file=pfile)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  ## note the typo here these are actually 13 trait basis projections !!!
  ## leave in so that we don't get confused.
  ofile <- file.path(OUT_DIR,sprintf("projections/%s_0619.RDS",trait))
  saveRDS(all.proj,file=ofile)
}

samples.DT <- rbindlist(samples)
saveRDS(samples.DT,file=file.path(OUT_DIR,'sample.RDS'))
