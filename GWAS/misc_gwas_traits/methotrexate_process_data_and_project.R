## bipolar prognosis GWAS
library(annotSnpStats)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
DATA.DIR <- '/home/ob219/share/Data/GWAS-summary/aav_limy_wong'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- '/home/ob219/share/as_basis/GWAS/mtx/mtx_matura_0619.RDS'


met.DT <- fread("zcat /home/ob219/share/Data/GWAS-summary/MATURA-MTXresponse-DAS.meta.gz")
met.DT[,pid:=paste(CHR,BP,sep=':')]
met.DT[,z:=qnorm(P/2,lower.tail=FALSE) * sign(BETA)]
## back compute var.beta
met.DT[,var.beta:=(BETA/z)^2]
## use this to compute sdY using n=1424 and MAF estimates from UK10K

## note OR are with respect to A1
man.DT <- fread(SNP_MANIFEST)
M <- merge(met.DT[,.(pid,a1=A1,a2=A2,beta=BETA,var.beta)],man.DT,by='pid')
M[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
## compute Z scores and therefore estimate of var(beta)


## note that these data are based on linear regressions examining CRP level
## changes from baseline with methotrexate treatment. Need to work out rescaling
## such that we can compare on the same scale.

## in fact we just project as is. The trick is then computing the significance
## for case control we multiply by the number of cases. For quant traits we
## multiply by sdY to get the variance of the residual (under the null) - which
## will be the case for a majority of SNPs

## chris has code in coloc to do this

##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##'
##' @title Estimate trait variance, internal function
##' @param vbeta vector of variance of coefficients
##' @param maf vector of MAF (same length as vbeta)
##' @param n sample size
##' @return estimated standard deviation of Y
##'
##' @author Chris Wallace
sdY.est <- function(vbeta, maf, n) {
    oneover <- 1/vbeta
    nvx <- 2 * n * maf * (1-maf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    return(sqrt(cf))
}

sdy <- M[,sdY.est(var.beta,maf,n=1424)]

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
M <- M[g.class!='match',beta:=beta * -1]
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)
tmp <- merge(M,stmp,by='pid',all.y=TRUE)
tmp$metric <- tmp[['ws_emp_shrinkage']] * tmp$beta
## where snp is missing make it zero
tmp[is.na(metric),metric:=0]
tmp[,trait:= 'methotrexate']
saveRDS(tmp,file='/home/ob219/share/as_basis/GWAS/mtx/mtx_source.RDS')
B <- dcast(tmp,pid ~ trait,value.var='metric')
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
saveRDS(all.proj,file=OUT_FILE)
