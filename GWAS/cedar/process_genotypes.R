#!/usr/bin/env Rscript


odir <- "/home/ob219/share/as_basis/GWAS/cedar/gt"
idir <- "/home/cew54/share/Data/expr/Cedar-published"

message("\n! ---------------------------------------")
message("! ",date())
message("! cedar-genotypes.R running with args:")
message("! output to ",odir)
message("! ---------------------------------------\n")

## basis snps
library(magrittr)
library(data.table)
basis <- fread("~/share/as_basis/GWAS/snp_manifest/gwas_june.tab")
basis[,chr:=sub(":.*","",pid)]
basis[,pos:=as.numeric(sub(".*:","",pid))]


## genotype data
library(annotSnpStats)
A <- annot.read.plink(file.path(idir,"CEDAR_GENO","CEDAR_Genotypes"))
summary(sm(A))
summary(A@snps)
A <- A[ ,A@snps$chromosome <= 22 ]
cs <- col.summary(A)
plot(cs$Call.rate,cs$z.HWE)
A <- A[,cs$MAF>0.01 & abs(cs$z.HWE)<5 & cs$Call.rate > 0.95]
rs <- row.summary(A)
plot(rs$Call.rate,rs$Het)
A <- A[rs$Call.rate > 0.99 & rs$Het > 0.31,]
## against build 37, should match - but ~40% drop out!
table(basis$pos %in% A@snps$position)

## now input genotype dosages for imputed snps
X <- fread(paste("zcat", file.path(idir,"CEDAR_GENO","CEDAR_Imputed.dos.gz")))
M <- as.matrix(X[,-1])
dimnames(M) <- list(X$id,names(X)[-1])
X <- X[,.(id)]
X[,c("chr","pos","a1","a2"):=tstrsplit(id,"_")]
X[,pid:=paste(chr,pos,sep=":")]
table(basis$pid %in% X$pid)
11105/nrow(basis) # 4% of snps missing
rownames(M) <- X$pid
o <- intersect(X$pid,basis$pid)
X <- X[pid %in% o,]
M <- M[X$pid,]

## subset to samples in A
M <- t(M)
dim(M)
dim(A)
table(rownames(A) %in% rownames(M))
M <- M[rownames(A),] # match rows

XA <- as.data.table(snps(A))
XA[,pid:=paste(chromosome,position,sep=":")]
table(XA$pid %in% X$pid)
## subset to basis
use <- X$pid %in% basis$pid
X <- X[use,]
M <- M[,use]
use <- XA$pid %in% basis$pid
XA <- XA[use,]
A <- A[,use]
## remove those also in A
drop <- which(X$pid %in% XA$pid)
X <- X[-drop,]
M <- M[,-drop]
MM <- new("SnpMatrix",matrix(mean2g(M),nrow(M)))
dimnames(MM) <- dimnames(M)

MM <- cbind(MM,sm(A))
setnames(XA,c("allele.1","allele.2"),c("a1","a2"))
X <- rbind(X[,.(pid,a1,a2)],XA[,.(pid,a1,a2)])

dim(X)
dim(MM)
colnames(MM) <- X$pid
table(basis$pid %in% X$pid) # 5113 missing
save(MM,X,file=file.path(odir,"MM.RData"))
