
## can we load in the whole shooting match

# this command used to work but now seems to time out !
#all.1kg <- fread("bcftools view -H -q 0.01 /home/ob219/reference_panels/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",sep="\t")
all.1kg <- fread("/home/ob219/reference_panels/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.0.01MAF.txt",sep="\t")

### this reads in a map file and spits out myogen Z score file for use with ImpG
library(rtracklayer)
library(annotSnpStats)
nmo.DT <- fread("zcat /home/ob219/share/Data/GWAS-summary/NMO_summary_stats_2018/NMO/EBI/NMO.Combined.gz")
nmo.DT[,pid:=paste(Chr,Pos,sep=':')]
nmo.DT <- nmo.DT[!pid %in% nmo.DT[duplicated(pid),],]
nmo.DT[,id:=1:.N]
setnames(nmo.DT,make.names(names(nmo.DT)))
# note that joining by SNP name is worse than joining by position
#merged <- merge(nmo.DT,all.1kg[,.(kg.CHR=V1,kg.BP=V2,SNP=V3,REF=V4,ALT=V5)],by='SNP',all.x=TRUE)
#missing 1289
all.1kg[,pid:=paste(V1,V2,sep=":")]
merged.pid <- merge(nmo.DT,all.1kg[,.(kg.CHR=V1,kg.BP=V2,kg.SNP=V3,REF=V4,ALT=V5,pid)],by='pid')
## missing 598333
merged.pid<- merged.pid[!is.na(kg.SNP),.(CHR=Chr,BP=Pos,SNP=rsID,A1=toupper(Allele1),A2=toupper(Allele2),P,beta=Effect,kg.SNP,REF,ALT,N)]
## for impss it seems as if effect allele is a2
## looks as if a1 is the effect allele align so that OR 1kg ref allele is a1 and alternative is a2
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
## remove impossible and ambig as cannot align these
merged.pid <- merged.pid[!g.class %in% c('impossible','ambig')]
## effect allele is a1 if we have to flip it means that effect allele is allele 2 which we want
## thus we flip aligned alleles to get what we want rather than those that don't match
sw <- merged.pid$g.class %in% c("comp","match")
merged.pid[sw,beta:=beta * -1]

nmo.filt <- merged.pid[,.(MarkerName=kg.SNP,Allele1=REF,Allele2=ALT,Z=qnorm(P/2,lower.tail=FALSE) * sign(beta),n=N,P,OR=exp(beta))]
write.table(nmo.filt[!is.na(Z),],file='/home/ob219/rds/hpc-work/ssimp/nmo_combined_with_or.txt',row.names=FALSE,sep=" ",quote=FALSE)
nmo.filt <- merged.pid[,.(MarkerName=kg.SNP,Allele1=REF,Allele2=ALT,Z=qnorm(P/2,lower.tail=FALSE) * sign(beta),n=N)]
write.table(nmo.filt[!is.na(Z),],file='/home/ob219/rds/hpc-work/ssimp/nmo_combined.txt',row.names=FALSE,sep=" ",quote=FALSE)

#NEED ALL THIS FOR IT TO WORK
#export LC_ALL=C; unset LANGUAGE




## check imputed vs input ?
if(FALSE){
  /home/ob219/git/ssimp_software/compiled/ssimp-linux-0.5.6 --gwas /home/ob219/rds/hpc-work/ssimp/nmo_combined.txt --ref 1KG/EUR --out /home/ob219/rds/hpc-work/ssimp/nmo_combined_imputed.txt --impute.maf 0.01
n0 <-
n1 <-

 imp<-fread("/home/ob219/rds/hpc-work/ssimp/nmo_combined.txt",select=c('SNP','source','r2.pred','z_imp','Allele1','Allele2','maf','N.imp'))
 input <- fread("/home/ob219/rds/hpc-work/ssimp/nmo_combined_with_or.txt")

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
M[,imp_beta_cc:=z_imp/sqrt(2 * ((n1 * n0)/(n1+n0)) * r2.pred * maf * (1-maf))]
aa <- ggplot(M,aes(x=log(OR),y=imp_beta_linear,color=r2.pred)) + geom_point() + geom_abline()
bb <- ggplot(M,aes(x=log(OR),y=imp_beta_cc,color=r2.pred)) + geom_point() + geom_abline()
}
