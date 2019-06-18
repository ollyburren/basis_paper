## bipolar prognosis GWAS
library(annotSnpStats)
library(rtracklayer)

myo.DT <- fread("zcat ~/share/Data/GWAS-summary/MYOGEN/Sep2018_summary_meta.txt.gz")
snps.DT <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab')
myo.DT[,pid:=paste(CHR,BP,sep=':')]
myo.DT <- myo.DT[!pid %in% myo.DT[duplicated(pid),],]
myo.DT[,id:=1:.N]
setnames(myo.DT,make.names(names(myo.DT)))



myo.36.gr <- with(myo.DT,GRanges(seqnames=Rle(paste0('chr',CHR)),ranges=IRanges(start=BP,width=1L),id=id))
c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
myo.37.gr<-unlist(liftOver(myo.36.gr,c))
DT.37 <- data.table(id=myo.37.gr$id,position.37=start(myo.37.gr))
myo.DT <- merge(myo.DT,DT.37,by.x='id',by.y='id',all.x=TRUE)
## 173 myo.DT don't match after coord conversion
myo.DT <- myo.DT[!is.na(position.37),]
myo.DT[,pid.37:=paste(CHR,position.37,sep=":")]

myo.DT[,c('a1','a2'):=list(A1_risk,ifelse(A1_risk!=A1_minor,A1_minor,A2_major))]

out <- myo.DT[,.(pid,a1,a2,or=OR,p.value=P)]


M <- merge(snps.DT,myo.DT,by.y='pid.37',by.x='pid')

M <- M[,.(pid,ref_a1,ref_a2,ref_a1.af,A1_minor,MAF_jdm=MAF_dmjdmpm,P,OR,a1,a2)]
M[is.na(OR),c('A1_minor','MAF_jdm','P','OR','a1','a2','missing'):=list(ifelse(ref_a1.af<0.5,ref_a1,ref_a2),
  ref_a1.af,
  0.99,1,a1=ref_a1,a2=ref_a2,TRUE)]

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
M[g.class=='comp',c('a1','a2'):=list(ref_a1,ref_a2)]
sw <- align.class %in% c("rev","revcomp")
# M[sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,1/OR,1-MAF_jdm)]
# M[!sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,OR,MAF_jdm)]

## note that here effect is wrt to a1 and not allele two - the basis assumes
## that effect is wrt to allele2 therefore revcomp and rev are OK
M[!sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,1/OR,1-MAF_jdm)]
M[sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,OR,MAF_jdm)]
M[align.class=='impossible',c('OR','P'):=list(NA,NA)]


SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_vit_t2d.RDS'
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)
tmp <- merge(M,stmp,by='pid',all.y=TRUE)
tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
## where snp is missing make it zero
tmp[is.na(metric),metric:=0]
#tmp[,trait:= 'myositis_myogen']
tmp[,trait:= 'myositis_myogen']
B <- dcast(tmp,pid ~ trait,value.var='metric')
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_vit_t2d.RDS'
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
#saveRDS(all.proj,file='/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis.RDS')
saveRDS(all.proj,file='/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis_vit_t2d.RDS')

## process subset data

for (trait in c('jdm','pm','dm')){
  message(trait)
  meta <- fread(sprintf("~/share/Data/GWAS-summary/MYOGEN/Dec2018/gwas_%s.meta",trait))
  fq <- fread(sprintf("~/share/Data/GWAS-summary/MYOGEN/Sep2018/%s.frq",trait))
  myo.DT <- merge(meta[,.(CHR,BP,SNP,A1,P,OR)],fq[,.(SNP,fq.A1=A1,fq.A2=A2,MAF)],by='SNP')
  myo.DT[,A2:=ifelse(A1==fq.A1,fq.A2,fq.A1)]
  snps.DT <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab')
  myo.DT[,pid:=paste(CHR,BP,sep=':')]
  myo.DT <- myo.DT[!pid %in% myo.DT[duplicated(pid),],]
  myo.DT[,id:=1:.N]
  setnames(myo.DT,make.names(names(myo.DT)))
  myo.36.gr <- with(myo.DT,GRanges(seqnames=Rle(paste0('chr',CHR)),ranges=IRanges(start=BP,width=1L),id=id))
  c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
  myo.37.gr<-unlist(liftOver(myo.36.gr,c))
  DT.37 <- data.table(id=myo.37.gr$id,position.37=start(myo.37.gr))
  myo.DT <- merge(myo.DT,DT.37,by.x='id',by.y='id',all.x=TRUE)
  ## 173 myo.DT don't match after coord conversion
  myo.DT <- myo.DT[!is.na(position.37),]
  myo.DT[,pid.37:=paste(CHR,position.37,sep=":")]
  myo.DT[,c('a1','a2'):=list(A1,A2),]
  M <- merge(snps.DT,myo.DT,by.y='pid.37',by.x='pid')
  M <- M[,.(pid,ref_a1,ref_a2,ref_a1.af,MAF,P,OR,a1,a2)]
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
  M[g.class=='comp',c('a1','a2'):=list(ref_a1,ref_a2)]
  sw <- align.class %in% c("rev","revcomp")
  M[!sw,c('a1','a2','or'):=list(ref_a1,ref_a2,1/OR)]
  M[sw,c('a1','a2','or'):=list(ref_a1,ref_a2,OR)]
  M[align.class=='impossible',c('OR','P'):=list(NA,NA)]
  SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_vit_t2d.RDS'
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  ## where snp is missing make it zero
  tmp[is.na(metric),metric:=0]
  tmp[,trait:= sprintf("%s_myogen",trait)]
  of <- sprintf('/home/ob219/share/as_basis/GWAS/myogen_myositis/%s_myositis_source.RDS',trait)
  saveRDS(tmp,file=of)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_vit_t2d.RDS'
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  saveRDS(all.proj,file=sprintf('/home/ob219/share/as_basis/GWAS/myogen_myositis/%s_myositis_vit_t2d.RDS',trait))
}
