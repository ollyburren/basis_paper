## bipolar prognosis GWAS
library(annotSnpStats)
library(rtracklayer)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- '/home/ob219/share/as_basis/GWAS/psa_aterido/psa_aterido_0619.RDS'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'


psa.DT <- fread("/home/ob219/share/Data/GWAS-summary/psa-aterido/psa_Aterido.csv")
psa.DT[,pid:=paste(CHR,POS,sep=':')]
psa.DT <- psa.DT[!pid %in% psa.DT[duplicated(pid),],]
psa.DT[,id:=1:.N]
## note OR are with respect to A1
man.DT <- fread(SNP_MANIFEST)
am.DT <- merge(psa.DT[,.(pid,a1=A1,a2=A2,or=ORN,p.value=PN)],man.DT,by='pid')
span.DT <- merge(psa.DT[,.(pid,a1=A1,a2=A2,or=ORS,p.value=PS)],man.DT,by='pid')

##alignment will be the same for both studies
alleles <- data.table(pid=am.DT$pid,al.x = paste(am.DT$ref_a1,am.DT$ref_a2,sep='/'),al.y=paste(am.DT$a1,am.DT$a2,sep='/'))
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

am.DT <- merge(am.DT,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
span.DT <- merge(span.DT,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
## so here alleles match we need to flip as we want wrt to a2
am.DT[g.class=='match',or:=1/or]
span.DT[g.class=='match',or:=1/or]

am.DT[,trait:='na_psa']
span.DT[,trait:='span_psa']

sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]

res <- lapply(list(N=am.DT,S=span.DT),function(M){
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  tra <- tmp[!is.na(trait),]$trait %>% unique
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",tra))
  saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  tmp[,trait:=tra]
  ## where snp is missing make it zero
  tmp[is.na(metric),metric:=0]
  saveRDS(tmp,file=sprintf('%s/%s_source.RDS',SRC_OUT_DIR,tra))
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
}) %>% do.call('rbind',.)

saveRDS(res,file=OUT_FILE)
