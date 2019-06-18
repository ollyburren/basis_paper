## bipolar prognosis GWAS
library(annotSnpStats)

bip.DT <- fread("zcat /home/ob219/share/Data/GWAS-summary/psyc_traits/BDvsCONT.sumstats.gz")
bip.DT[,pid:=paste(CHR,BP,sep=':')]
## note OR are with respect to A1
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
man.DT <- fread(SNP_MANIFEST)
M <- merge(bip.DT[,.(pid,a1=A1,a2=A2,or=OR)],man.DT,by='pid')
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

## check direction which is the effect allele ? It appears that allele_B/a2 is the effect allele
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
M <- M[!duplicated(pid),]
M <- M[g.class=='match',or:=1/or]
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_vit_t2d.RDS'
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)
tmp <- merge(M,stmp,by='pid',all.y=TRUE)
tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
## where snp is missing make it zero
tmp[is.na(metric),metric:=0]
tmp[,trait:= 'BIP']
saveRDS(tmp,file='/home/ob219/share/as_basis/GWAS/psych/BIP_source.RDS')
B <- dcast(tmp,pid ~ trait,value.var='metric')
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_vit_t2d.RDS'
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
saveRDS(all.proj,file='/home/ob219/share/as_basis/GWAS/psych/BIP_vit_t2d.RDS')
