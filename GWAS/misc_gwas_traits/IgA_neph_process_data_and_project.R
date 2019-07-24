## bipolar prognosis GWAS
library(annotSnpStats)
library(rtracklayer)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- '/home/ob219/share/as_basis/GWAS/IgA_nephropathy/IgA_nephropathy_0619.RDS'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'


iga.DT <- fread("zcat ~/share/Data/GWAS-summary/IgA_nephropathy/METAANALYSIS1.4cohorts.ORDERED")
iga.DT[,pid:=paste(CHR,BP,sep=':')]
iga.DT <- iga.DT[!pid %in% iga.DT[duplicated(pid),],]
iga.DT[,id:=1:.N]


iga.36.gr <- with(iga.DT,GRanges(seqnames=Rle(paste0('chr',CHR)),ranges=IRanges(start=BP,width=1L),id=id))
c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
iga.37.gr<-unlist(liftOver(iga.36.gr,c))
DT.37 <- data.table(id=iga.37.gr$id,position.37=start(iga.37.gr))
iga.DT <- merge(iga.DT,DT.37,by.x='id',by.y='id',all.x=TRUE)
## 173 iga.DT don't match after coord conversion
iga.DT <- iga.DT[!is.na(position.37),]
iga.DT[,pid.37:=paste(CHR,position.37,sep=":")]

## note OR are with respect to A1
man.DT <- fread(SNP_MANIFEST)
M <- merge(iga.DT[,.(pid=pid.37,a1=toupper(Allele1),a2=toupper(Allele2),or=exp(Effect),p.value=P)],man.DT,by='pid')

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
## so here alleles match we need to flip as we want wrt to a2
M <- M[g.class=='match',or:=1/or]
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)
tmp <- merge(M,stmp,by='pid',all.y=TRUE)
tra <- 'IgA_nephropathy'
pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",tra))
saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
## where snp is missing make it zero
tmp[is.na(metric),metric:=0]
tmp[,trait:= 'IgA_nephropathy']
saveRDS(tmp,file='/home/ob219/share/as_basis/GWAS/IgA_nephropathy/IgA_nephropathy_source.RDS')
B <- dcast(tmp,pid ~ trait,value.var='metric')
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
saveRDS(all.proj,file=OUT_FILE)
