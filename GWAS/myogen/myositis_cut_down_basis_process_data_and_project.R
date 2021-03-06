library(annotSnpStats)
library(rtracklayer)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_myositis_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_myositis_gwas_0619.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_myositis_0619.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/myositis_ss_av_0619.RDS'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'


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


sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)
tmp <- merge(M,stmp,by='pid',all.y=TRUE)

pfile <- file.path(SRC_OUT_DIR,sprintf("%s_myo_snps_only_source.RDS",'myositis_myogen'))
saveRDS(tmp[,.(pid,or,p.value=P,ws_emp_shrinkage)],file=pfile)


tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
## where snp is missing make it zero
tmp[is.na(metric),metric:=0]
#tmp[,trait:= 'myositis_myogen']
tmp[,trait:= 'myositis_myogen']
B <- dcast(tmp,pid ~ trait,value.var='metric')
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
#saveRDS(all.proj,file='/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis.RDS')
saveRDS(all.proj,file='/home/ob219/share/as_basis/GWAS/myogen_myositis/myo_snps_only_myogen_myositis_0619.RDS')

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
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_myo_snps_only_source.RDS",sprintf("%s_myogen",trait)))
  saveRDS(tmp[,.(pid,or,p.value=P,ws_emp_shrinkage)],file=pfile)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  ## where snp is missing make it zero
  tmp[is.na(metric),metric:=0]
  tmp[,trait:= sprintf("%s_myogen",trait)]
  #of <- sprintf('/home/ob219/share/as_basis/GWAS/myogen_myositis/%s_myositis_source.RDS',trait)
  #saveRDS(tmp,file=of)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  ofile <- sprintf("/home/ob219/share/as_basis/GWAS/myogen_myositis/%s_myo_snps_only_0619.RDS",trait)
  message(ofile)
  saveRDS(all.proj,file=ofile)
}


## next collate all the results and compute correct Z scores etc.

myositis <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myo_snps_only_myogen_myositis_0619.RDS')
myositis <- data.table(trait=rownames(myositis),myositis)
myositis <- melt(myositis,id.var='trait')
myositis[,c('n0','n1'):=list(4724,1711)]
myositis[,category:='myositis']

## jdm

jdm <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/jdm_myo_snps_only_0619.RDS')
jdm <- data.table(trait=rownames(jdm),jdm)
jdm <- melt(jdm,id.var='trait')
jdm[,c('n0','n1'):=list(4724,473)]
jdm[,category:='myogen']

dm <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/pm_myo_snps_only_0619.RDS')
dm <- data.table(trait=rownames(dm),dm)
dm <- melt(dm,id.var='trait')
dm[,c('n0','n1'):=list(4724,705)]
dm[,category:='myogen']

pm <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/dm_myo_snps_only_0619.RDS')
pm <- data.table(trait=rownames(pm),pm)
pm <- melt(pm,id.var='trait')
pm[,c('n0','n1'):=list(4724,533)]
pm[,category:='myogen']

myogen <- rbindlist(list(myositis,jdm,dm,pm))
myogen[,n:=n1+n0]

#VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_0619.RDS'
var.DT <- readRDS(VARIANCE_FILE)
#BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,value:=(value-control.loading)]
basis.DT[,Z:=sign(value)]
basis.DT[,p.adj:=1]
basis.DT[,control.loading:=NULL]

## compute the variance of a projection
myogen.DT <- merge(myogen,var.DT,by.x='variable',by.y='pc')
myogen.DT <- merge(myogen.DT,control.DT,by.x='variable',by.y='PC')
myogen.DT[,variance:=((log(n)-(log(n1) + log(n-n1)))+ log(mfactor)) %>% exp]
#myogen.DT[!is.na(sdy),variance:=(log(sdy^2/n) + log(mfactor)) %>% exp]
## add in control loading
myogen.DT[,Z:=(value-control.loading)/sqrt(variance)]
myogen.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
myogen.DT[,p.adj:=p.adjust(p.value,method="bonferroni"),by='variable']
myogen.DT[,delta:=value-control.loading]
saveRDS(myogen.DT,'/home/ob219/share/as_basis/GWAS/RESULTS/24_07_19_0619_summary_results_myogen_snps_only.RDS')
