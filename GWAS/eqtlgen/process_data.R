## process both cis and trans signals from Võsa,U. et al. (2018) Unraveling the polygenic architecture of complex traits using blood eQTL meta-analysis. bioRxiv, 447367.

## note that these are large files - I used UNIX to create a list of unique variants in
## order to facilitate filtering.

EQTL_GEN_MAN <- '/home/ob219/share/Data/expr/eqtlgen/cis_manifest.tab'
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
e.DT <- fread(EQTL_GEN_MAN)
setnames(e.DT,c('id','chr','position','effect','other'))
e.DT <- e.DT[id!='SNP',]
e.DT[,pid:=paste(chr,position,sep=':')]
snp.DT <- fread(SNP_MANIFEST)
keep.id <- e.DT[pid %in% snp.DT$pid,]$id
ID_FILE <- "/home/ob219/share/as_basis/GWAS/eqtlgen/gwas_june_match.tab"
write(keep.id,file=ID_FILE)
CIS_FILE <- '/home/ob219/share/Data/expr/eqtlgen/cis-eQTLs_full_20180905.txt.gz'
OUT_FILE <- '/home/ob219/share/as_basis/GWAS/eqtlgen/cis-eQTLs_full_20180905.filtered.txt'
cmd <- sprintf("zcat %s | grep -f %s > %s",CIS_FILE,ID_FILE,OUT_FILE)
system(cmd)
#zcat /home/ob219/share/Data/expr/eqtlgen/cis-eQTLs_full_20180905.txt.gz | grep -f /home/ob219/share/as_basis/GWAS/eqtlgen/gwas_june_match.tab > /home/ob219/share/as_basis/GWAS/eqtlgen/cis-eQTLs_full_20180905.filtered.txt
## above does not work as get out of mem error instead us data.table
DT<-fread("zcat /home/ob219/share/Data/expr/eqtlgen/cis-eQTLs_full_20180905.txt.gz")
f.DT<-DT[SNP %in% keep.id,]

saveRDS(f.DT,"/home/ob219/share/as_basis/GWAS/eqtlgen/cis-eQTLs_full_20180905.filtered.RDS")
f.DT<-readRDS("/home/ob219/share/as_basis/GWAS/eqtlgen/cis-eQTLs_full_20180905.filtered.RDS")
## next do the same for trans eqtl
EQTL_GEN_MAN_TRANS <- '/home/ob219/share/Data/expr/eqtlgen/trans_manifest.tab'
DT.t<-fread("/home/ob219/share/Data/expr/eqtlgen/trans-eQTLs_full_20180905.txt")
e.DT <- fread(EQTL_GEN_MAN_TRANS)
setnames(e.DT,c('id','chr','position','effect','other'))
e.DT <- e.DT[id!='SNP',]
e.DT[,pid:=paste(chr,position,sep=':')]
snp.DT <- fread(SNP_MANIFEST)
keep.id <- e.DT[pid %in% snp.DT$pid,]$id
f.DT.t<-DT.t[SNP %in% keep.id,]

## merge into one file and convert Z scores to the beta scale by multiplying through by se of beta due to variance

f.DT.t[,type:='trans']
f.DT[,type:='cis']

basis.DT <- rbind(f.DT.t,f.DT)

basis.DT[,pid:=paste(SNPChr,SNPPos,sep=':')]
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
OUT_DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/sum_stats/eqtlgen/'
for(n in names(by.gene)){
  message(n)
  out <- by.gene[[n]][,.(pid,a1,a2,or,p.value)]
  fname <- file.path(OUT_DIR,sprintf("%s.tab",n))
  write.table(out,file=fname,quote=FALSE,row.names=FALSE,sep="\t")
}
