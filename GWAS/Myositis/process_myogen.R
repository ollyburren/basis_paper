library(annotSnpStats)
library(rtracklayer)

myo.DT <- fread("zcat ~/share/Data/GWAS-summary/MYOGEN/Sep2018_summary_meta.txt.gz")
snps.DT <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab')
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

## rather than being risk the A1 allele is effect allele as this is plink ! Current software
## by default assumes that the effect allele is A2 so really this should be flipped !
myo.DT[,c('a1','a2'):=list(A1_risk,ifelse(A1_risk!=A1_minor,A1_minor,A2_major))]

out <- myo.DT[,.(pid,a1,a2,or=OR,p.value=P)]


M <- merge(snps.DT,myo.DT,by.y='pid.37',by.x='pid',all.x=TRUE)

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
M[sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,1/OR,1-MAF_jdm)]
M[align.class=='impossible',c('OR','P'):=list(NA,NA)]

## remove those that don't conform

mod <- lm(data=M,ref_a1.af~MAF_jdm)
## use above to get STDERR of residuals to get outliers
M[,ab_maf_dif:=abs(ref_a1.af - MAF_jdm)]
## 0.999% of the data qnorm(0.999,lower.tail=TRUE)
M[,ba:=ifelse(ab_maf_dif>(0.008194 * 3) & ((ref_a1.af>0.5 & MAF_jdm<0.5) | (ref_a1.af<0.5 & MAF_jdm>0.5)),TRUE,FALSE)]
## set the above to missing i.e. set or to 1 and p.value to 0.99
M[ba==TRUE,c('OR','p.value'):=list(1,0.999)]
## next write these out so we can project !!!
out <- M[,.(pid,a1=ref_a1,a2=ref_a2,or=OR,p.value=P)]
write.table(out,file="/home/ob219/share/as_basis/GWAS/sum_stats/myositis_myogen.tab",quote=FALSE,sep="\t",row.names=FALSE)
