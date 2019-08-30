
## can we load in the whole shooting match

# this command used to work but now seems to time out !
#all.1kg <- fread("bcftools view -H -q 0.01 /home/ob219/reference_panels/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",sep="\t")
all.1kg <- fread("/home/ob219/reference_panels/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.0.01MAF.txt",sep="\t")

### this reads in a map file and spits out myogen Z score file for use with ImpG
library(rtracklayer)
library(annotSnpStats)
myo.DT <- fread("/home/ob219/share/Data/GWAS-summary/MYOGEN/Dec2018/gwas_dm.meta")
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

# note that joining by SNP name is worse than joining by position
#merged <- merge(myo.DT,all.1kg[,.(kg.CHR=V1,kg.BP=V2,SNP=V3,REF=V4,ALT=V5)],by='SNP',all.x=TRUE)
#missing 1289
all.1kg[,pid.37:=paste(V1,V2,sep=":")]
merged.pid <- merge(myo.DT,all.1kg[,.(kg.CHR=V1,kg.BP=V2,kg.SNP=V3,REF=V4,ALT=V5,pid.37)],by='pid.37',all.x=TRUE)
## missing 676
## try and get alternative allele a different way using the frq files as we do for subphenotypes
fq <- fread("~/share/Data/GWAS-summary/MYOGEN/Sep2018/dm.frq")
merged.pid<- merge(merged.pid[,.(CHR,BP,SNP,A1,P,OR,kg.SNP,REF,ALT)],fq[,.(SNP,fq.A1=A1,fq.A2=A2,MAF,n=NCHROBS/2)],by='SNP')
merged.pid[,A2:=ifelse(A1==fq.A1,fq.A2,fq.A1)]
## for impss it seems as if effect allele is a2 from the note in the myogen dir it seems as if a1 is the effect allele
## it says risk allele but this makes no sense as sometimes the A1_risk is the minor allele but yet the OR is < 1 so
## assume it is the effect allele
#myo.filt <- myo.DT[CHR==1,.(MarkerName=SNP,Allele1=A2,Allele2=A1,Z=qnorm(P/2,lower.tail=FALSE) * sign(OR),n)]
#write.table(myo.filt[!is.na(Z),],file='/home/ob219/rds/hpc-work/ssimp/myogwas_chr1.txt',row.names=FALSE,sep=" ",quote=FALSE)
## a1 is the effect allele align so that OR 1kg ref allele is a1 and alternative is a2
alleles <- data.table(al.x = paste(merged.pid$REF,merged.pid$ALT,sep='/'),al.y=paste(merged.pid$A1,merged.pid$A2,sep='/'))
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$pid
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
merged.pid[,g.class:=align.class]
merged.pid <- merged.pid[g.class!='impossible',]
## effect allele is a1 if we have to flip it means that effect allele is allele 2 which we want
## thus we flip aligned alleles to get what we want rather than those that don't match
sw <- merged.pid$g.class %in% c("comp","match")
merged.pid[sw,OR:=1/OR]

myo.filt <- merged.pid[,.(MarkerName=SNP,Allele1=REF,Allele2=ALT,Z=qnorm(P/2,lower.tail=FALSE) * sign(log(OR)),n,P,OR)]
write.table(myo.filt[!is.na(Z),],file='/home/ob219/rds/hpc-work/ssimp/myo_dm_gwas_with_or.txt',row.names=FALSE,sep=" ",quote=FALSE)
myo.filt <- merged.pid[,.(MarkerName=SNP,Allele1=REF,Allele2=ALT,Z=qnorm(P/2,lower.tail=FALSE) * sign(log(OR)),n)]
write.table(myo.filt[!is.na(Z),],file='/home/ob219/rds/hpc-work/ssimp/myo_dm_gwas.txt',row.names=FALSE,sep=" ",quote=FALSE)

#NEED ALL THIS FOR IT TO WORK
#export LC_ALL=C; unset LANGUAGE

/home/ob219/git/ssimp_software/compiled/ssimp-linux-0.5.6 --gwas /home/ob219/rds/hpc-work/ssimp/myo_dm_gwas.txt --ref 1KG/EUR --out /home/ob219/rds/hpc-work/ssimp/myo_dm_gwas_imputed.txt --impute.maf 0.01


## check imputed vs input ?

 imp<-fread("/home/ob219/rds/hpc-work/ssimp/myo_dm_gwas_imputed.txt",select=c('SNP','source','r2.pred','z_imp','Allele1','Allele2','maf','N.imp'))
 input <- fread("/home/ob219/rds/hpc-work/ssimp/myo_dm_gwas_with_or.txt")

 M<-merge(input[,.(SNP=MarkerName,orig.z=Z,in.a1=Allele1,in.a2=Allele2,OR)],imp[,.(SNP,source,r2.pred,z_imp,out.a1=Allele1,out.a2=Allele2,maf,N.imp)],by='SNP')
 alleles <- data.table(al.x = paste(M$out.a1,M$out.a2,sep='/'),al.y=paste(M$in.a1,M$in.a2,sep='/'))
 align.class <- rep('match',nrow(alleles))
 idx<-which(alleles$al.x!=alleles$al.y)
 x.alleles <- alleles[idx,]$al.x
 names(x.alleles)<-alleles[idx,]$pid
 y.alleles <-  alleles[idx,]$al.y
 names(y.alleles)<-names(x.alleles)
 align.class[idx] <- g.class(x.alleles,y.alleles)
 print(table(align.class))
 M[,g.class:=align.class]
 sw <- align.class %in% c("rev","revcomp")
 M[sw,orig.z:=orig.z * -1]
 M[sw,OR:=1/OR]
 library(cowplot)
 ggplot(M,aes(x=orig.z,y=z_imp,color=r2.pred)) + geom_point()
 ## next compute the OR derived from the imputed score
M[,imp_beta_linear:=z_imp/sqrt(2 * N.imp * maf * (1-maf))]
M[,imp_beta_cc:=z_imp/sqrt(2 * ((705 * 4724)/(705+4724)) * r2.pred * maf * (1-maf))]
aa <- ggplot(M,aes(x=log(OR),y=imp_beta_linear,color=r2.pred)) + geom_point() + geom_abline()
bb <- ggplot(M,aes(x=log(OR),y=imp_beta_cc,color=r2.pred)) + geom_point() + geom_abline()

## do new calculation with the reimputed data

reimp <- merge(input[,.(SNP=MarkerName,orig.z=Z,in.a1=Allele1,in.a2=Allele2,OR)],imp[!is.na(Z_reimputed),.(SNP,source,r2.pred,z_imp,out.a1=Allele1,out.a2=Allele2,maf,N.imp,Z_reimputed,r2_reimputed)],by='SNP')
reimp[,imp_beta:=Z_reimputed/sqrt(2 * ((1711 * 4724)/(1711+4724)) * r2_reimputed * maf * (1-maf))]
ggplot(reimp,aes(x=log(OR),y=imp_beta,color=r2_reimputed)) + geom_point() + geom_abline()


 ## recode the sanity check from the software providers in data.table converted from tidyverse

dat.merge <- merge(imp[!is.na(Z_reimputed),.(SNP, Allele1, Allele2, z_imp, r2.pred, Z_reimputed, r2_reimputed)],input[,.(SNP = MarkerName, Z.GWAS = Z, Allele1.GWAS = Allele1, Allele2.GWAS=Allele2)],by='SNP')
dat.merge[Allele1.GWAS != Allele1 & Allele2.GWAS != Allele2,Z.GWAS:=-1 * Z.GWAS]
qplot(Z.GWAS, Z_reimputed, data = dat.merge, color = r2_reimputed) + geom_abline(intercept = 0, slope = 1) ## plot the identity line

## note, that z.imp is the imputation of the tag SNP, while leaving the SNP in the set of tag SNPs. hence r2.pred = 1.
