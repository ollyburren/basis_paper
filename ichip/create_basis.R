## code to create and store basis and project all summary statistics

#library(devtools)
#install_github('ollyburren/cupcake')
load_all('~/git/cupcake')
library(cupcake)

SHRINKAGE_METHOD<-'ws_emp'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/ichip/support/shrinkage_ic.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/ichip/support/basis_ic.RDS'
ICHIP_DATA_DIR <- '/home/ob219/share/as_basis/ichip/sum_stats'
SNP_MANIFEST_FILE <- '/home/ob219/share/as_basis/ichip/snp_manifest/ichip_all.tab'
MANIFEST <- '/home/ob219/share/as_basis/ichip/trait_manifest/as_manifest_ichip.tsv'
VARIANCE_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/analytical_variances_ichip.RDS'
DEFAULT.SNPSTATS.DIR <- '/home/ob219/rds/hpc-work/as_basis/snpStats/basis_ichip'

## load data
basis.DT<-get_gwas_data(MANIFEST,SNP_MANIFEST_FILE,ICHIP_DATA_DIR,filter_snps_by_manifest=TRUE)
## compute various shrinkage methods and store
shrink.DT<-compute_shrinkage_metrics(basis.DT)
saveRDS(shrink.DT,file=SHRINKAGE_FILE)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
saveRDS(pc.emp,file=BASIS_FILE)

if(FALSE){
  library(cowplot)
  library(ggrepel)
  pc.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
  ggplot(pc.DT,aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text_repel()
}


## compute the variances
gw.DT <- basis.DT[trait==head(sample(unique(trait)),n=1),]
## set the standard error to the null
gw.DT[,emp_se:=se_null(n,n1,maf)]
## data table of PC snp loadings
gw.DT[,c('chr','position'):=tstrsplit(pid,':')]
w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)

if(FALSE){
    #load in manifest file and prepare annotSnpStats files for simulations
    sman.DT <- fread(SNP_MANIFEST_FILE)
    ## code to generate ichip control genotype data object for computing variance for a projection given sample size
    library(annotSnpStats)
    load("~/share/Projects/coloccc/aligned9.build37.RData")
    ## downsample to make simulations quicker
    #n<-1000
    #sample.keep<-sample(1:nrow(X),n) %>% sort
    #X.f <- X[sample.keep,]
    snp.DT <- snps(X.f) %>% data.table(pid=rownames(snps(X)),.)
    keep.snps <- which(snp.DT$pid %in% sman.DT$pid)
    X.f <- X[,keep.snps]
    snp.DT <- snp.DT[keep.snps,]
    ## filter object by manifest file
    by.chr<-split(1:nrow(snp.DT),snp.DT$chromosome)
    out.dir <- '/home/ob219/rds/hpc-work/as_basis/snpStats/basis_ichip'
    for(chr in names(by.chr)){
      message(sprintf("Processing %s",chr))
      tmp <- X.f[,by.chr[[chr]]]
      snps.DT <- snps(tmp) %>% data.table(pid=rownames(snps(tmp)),.)
      obj <- list(sm=sm(tmp),info=snps.DT)
      fname <- file.path(out.dir,sprintf("%s.RData",chr))
      save(obj,file=fname)
    }
}



basis.DT<-get_gwas_data(MANIFEST,SNP_MANIFEST_FILE,ICHIP_DATA_DIR,filter_snps_by_manifest=TRUE)
pc.emp <- readRDS(BASIS_FILE)
shrink.DT <- readRDS(SHRINKAGE_FILE)
gw.DT <- basis.DT[trait==head(sample(unique(trait)),n=1),]
## set the standard error to the null
gw.DT[,emp_se:=se_null(n,n1,maf)]
## these as expected are equivalent
gw.DT[,maf_se:=(maf_se_estimate(maf)/2) * sqrt(n/((n-n1)*n1)) * sqrt(1/2)]
ggplot(gw.DT,aes(x=emp_se,y=maf_se)) + geom_point() + geom_abline()


## does this work i.e. removing reference to sample size alltogether ?
gw.DT[,emp_se:=sqrt(1/maf + 1/(1-maf)) * sqrt(1/2)]


## data table of PC snp loadings
gw.DT[,c('chr','position'):=tstrsplit(pid,':')]
w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
analytical.vars <- compute_proj_var(gw.DT,w.DT,shrink.DT,DEFAULT.SNPSTATS.DIR,SHRINKAGE_METHOD,quiet=FALSE)
## we can convert between by normalising by our factor total.num/(ncases.num * (total.num-ncases.num))
# N <- 2000
# N1 <- 1000
# factor <- N/(N1 * (N-N1))

factor <- 1
#adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars/factor)
adj.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars)
#adj.vars <- split(adj.vars,names(adj.vars))
setkey(adj.vars,pc)
saveRDS(adj.vars,file=VARIANCE_FILE)




## plotting code - delete when happy as not required to create the basis
library(cowplot)
library(ggrepel)
pc.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
N <- 50000
N1 <- 1000
factor <- N/(N1 * (N-N1))


computeEllipse <- function(a,b,x,y,n.points=1000){
  rad <- seq(0,2 * pi,length.out=n.points)
  xe <- (a * cos(rad)) + x
  ye <- (b * sin(rad)) + y
  elipse.DT <- data.table(x=xe,y=ye)
}

fac <- factor

ctrl.ellipse <- computeEllipse(sqrt(adj.vars[1,]$mfactor * fac),sqrt(adj.vars[2,]$mfactor * fac),pc.DT[trait=='control',]$PC1,pc.DT[trait=='control',]$PC2,)
ggplot(pc.DT,aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text_repel() +
geom_line(data=ctrl.ellipse,aes(x=x,y=y),alpha=0.1,inherit.aes=FALSE)


cov_beta <- function(sm,se_lor){
  #beta_hat <- lor * lor_shrink
  if(length(se_lor)==1)
    return(se_lor^2)
  # compute R statistic
  #r<-ld(sm,sm,stats="R.squared")
  r<-ld(sm,sm,stats="R.squared")
  # compute closest pos-def covariance matrix
  r<-as.matrix(mvs_sigma(Matrix(r)))
  ## for beta the covariance matrix is estimates by sigma x SE * SE^T
  cov_se<-tcrossprod(se_lor)
  return(cov_se * r)
}

block<-gw.DT[ld.block==1977,]
ss.file<-"/home/ob219/rds/hpc-work/as_basis/snpStats/basis_ichip/10.RData"
sm<-get(load(ss.file))
dup.idx<-which(duplicated(obj$info$pid))
if(length(dup.idx)>0)
  message("Duplicated")
sm$info$order<-1:nrow(sm$info)
info <-sm$info[pid %in% block$pid ,.(pid,order)]
dat <- merge(info,block,by.x='pid',by.y='pid')
sm <- sm$sm[,dat$order]
r<-ld(sm,sm,stats="R")
cov_se<-tcrossprod(dat$emp_se)
#cov_se <- dat$emp_se %*% t(dat$emp_se)
Sigma <- cov_se * r
sum(Sigma)
sapply(1:10000,function(i) (rnorm(nrow(info),0,1) * dat$emp_se) %>% var) %>% mean


lor_se=se_null(2000,1000,0.3)
var(rnorm(10000,0,1) * lor_se)
## which is the same as
var(rnorm(10000,0,1)) * lor_se^2
## variance between
r<-1
cov_se <- tcrossprod(lor_se)
Sigma <- cov_se * r
sum(Sigma)

Z=rmvnorm(10000,rep(0,2),ld.r)
B=Z*lor_se
Bsum=rowSums(Z)
var(Bsum)


## simulate two SNPs
ld.r <- matrix(c(1,0.8,0.8,1),2,2)
lor_se=se_null(rep(2000,2),rep(1000,2),c(0.3,0.4))
apply(rmvnorm(10000,rep(0,2),ld.r),2,var) * lor_se^2
cov_se <- tcrossprod(lor_se)
(cov_se * ld.r) %>% sum
