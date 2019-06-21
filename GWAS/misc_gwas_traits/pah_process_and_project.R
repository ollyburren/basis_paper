library(annotSnpStats)

## 6 GWAS to project one of each of these for PR3 and MPO
#1 gwas1
#2 gwas2
#3 meta (gwas1 and gwas2)

## create a table of sample sizes for each GWAS
#fname label n1  n0
#bolt_gwas1_mpo_bgen.stats.gz  mpo_gwas1  264 5259
#bolt_gwas1_pr3_bgen.stats.gz  pr3_gwas1 478 5259
#bolt_gwas2_mpo_bgen.stats.gz  mpo_gwas2 609 6717
#bolt_gwas2_pr3_bgen.stats.gz  pr3_gwas2 1142  6717
#meta_mpo_lmm.txt  mpo_meta  873 11976
#meta_pr3_lmm.txt  pr3_meta  1620 11976
SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
DATA.DIR <- '/home/ob219/share/Data/GWAS-summary/aav_limy_wong'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- "/home/ob219/share/as_basis/GWAS/liley_pah/projections/pah_0619.RDS"

(load('/home/ob219/share/Data/GWAS-summary/PAH/results_pah.RData'))

pah.DT <- data.table(snp=rownames(px),px)
setnames(pah.DT,make.names(colnames(pah.DT)))
pah.DT[,c('pid','alleles'):=tstrsplit(snp,"\\_")]
pah.DT[,c('a1','a2'):=tstrsplit(alleles,'/')]

man.DT <- fread(SNP_MANIFEST)
M <- merge(pah.DT[,.(pid,a1,a2,or=exp(Estimate),n0=AN.nonPAH/2,n1=AN.PAH,ctrl.af=AF.nonPAH)],man.DT,by='pid')
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

sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)
tmp <- merge(M,stmp,by='pid',all.y=TRUE)
tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
## where snp is missing make it zero
tmp[is.na(metric),metric:=0]
tmp[,trait:= 'PAH']
saveRDS(tmp,file='/home/ob219/share/as_basis/GWAS/psych/adhd_source.RDS')
B <- dcast(tmp,pid ~ trait,value.var='metric')
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
saveRDS(all.proj,file=OUT_FILE)
