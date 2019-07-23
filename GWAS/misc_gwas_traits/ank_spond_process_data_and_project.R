## bipolar prognosis GWAS
library(annotSnpStats)
library(rtracklayer)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- '/home/ob219/share/as_basis/GWAS/ank_spond/ank_spond_0619.RDS'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'


as.DT <- fread("~/share/Data/GWAS-summary/ank-spond-internal-2019.csv")
as.DT[,p.value:=pnorm(abs(beta/sqrt(varbeta)),lower.tail=FALSE) * 2]
## note OR are with respect to A1
man.DT <- fread(SNP_MANIFEST)
M <- merge(as.DT[,.(pid,a1=allele.1,a2=allele.2,or=exp(beta),p.value)],man.DT,by='pid')

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

## check direction which is the effect allele ? It appears that a1 is the effect allele
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
M <- M[!duplicated(pid),]
## flip revcomp and rev aligned as we want wrt to a2
M <- M[g.class %in% c('rev','revcomp'),or:=1/or]
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)
tmp <- merge(M,stmp,by='pid',all.y=TRUE)
pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",'ank_spond'))
saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
## where snp is missing make it zero
tmp[is.na(metric),metric:=0]
#tmp[,trait:= 'ankylosing_spondylitis']
tmp[,trait:= 'ank_spond']
#saveRDS(tmp,file='/home/ob219/share/as_basis/GWAS/ank_spond/ank_spond_source.RDS')
B <- dcast(tmp,pid ~ trait,value.var='metric')
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
saveRDS(all.proj,file=OUT_FILE)
