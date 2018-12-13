## process Ig titre
library(annotSnpStats)
DATA_DIR <- '/home/ob219/share/Data/GWAS-summary/ig_titre_unpublished'
files <- list.files(path=DATA_DIR,pattern="*.txt",full.names=TRUE)
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
man.DT <- fread(SNP_MANIFEST)

f <- files[1]
i.DT <- fread(f)
i.DT[,pid:=paste(chr,pos,sep=':')]
i.DT <- i.DT[,.(pid,a1=allele_A,a2=allele_B,or=exp(beta),p.value=p_value)]
M <- merge(i.DT,man.DT,by='pid')
alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
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
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
M <- M[!duplicated(pid),]
