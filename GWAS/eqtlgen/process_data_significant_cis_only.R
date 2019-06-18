## process both cis and trans signals from Võsa,U. et al. (2018) Unraveling the polygenic architecture of complex traits using blood eQTL meta-analysis. bioRxiv, 447367.

## note that these are large files - I used UNIX to create a list of unique variants in
## order to facilitate filtering.

CIS_FILE <- '/home/ob219/share/Data/expr/eqtlgen/cis-eQTL_significant_20181017.txt.gz'
cis.DT <- fread(sprintf("zcat %s",CIS_FILE))
cis.DT[,pid:=paste(SNPChr,SNPPos,sep=':')]
## read in basis
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
snp.DT <- fread(SNP_MANIFEST)
basis.DT <- cis.DT[pid %in% snp.DT$pid,]
basis.DT <- basis.DT[,.(pid,Zscore,risk.allele=AssessedAllele,other.allele=OtherAllele,ensg=Gene,NrSamples)]
##merge with the manifest to get MAF so can convert to the correct scale
M <- merge(basis.DT,snp.DT,by='pid')

M[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]

# this to get a derived beta that takes into account differing MAF's and sample sizes.
M[,der.beta:=1/(sqrt(2 * NrSamples) * sqrt(maf * (1-maf))) * Zscore]

library(annotSnpStats)
alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$risk.allele,M$other.allele,sep='/'))
alleles <- alleles[!duplicated(pid),]
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
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
## for basis a2 is the risk allele this means flipped alleles are OK but
## matched alleles are the wrong way round
#M[g.class=='rev',der.beta:=der.beta*-1]
M[g.class=='match',der.beta:=der.beta*-1]
M.out <- M[,.(pid,a1=ref_a1,a2=ref_a2,or=der.beta,p.value=2*pnorm(abs(Zscore),lower.tail=FALSE),ensg)]
by.gene<-split(M.out,M.out$ensg)
#OUT_DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/sum_stats/eqtlgen/'
OUT_DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/sum_stats/eqtlgen_significant_only/'
for(n in names(by.gene)){
  message(n)
  out <- by.gene[[n]][,.(pid,a1,a2,or,p.value)]
  fname <- file.path(OUT_DIR,sprintf("%s.tab",n))
  write.table(out,file=fname,quote=FALSE,row.names=FALSE,sep="\t")
}
