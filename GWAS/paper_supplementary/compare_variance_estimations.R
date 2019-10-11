## Comparison of variance estimation for PC projection.

## Select a ukbb trait to estimate variance for
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_gwas_13_traits_0919.RDS'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'



## why is egpa so out of whack ?

dat <- readRDS('/home/ob219/share/as_basis/GWAS/for_fdr_13_traits_0919/egpa_source.RDS')
snp.DT <- fread(SNP_MANIFEST_FILE)
M <- merge(snp.DT,dat,by='pid')
## finally add the PC rotations
pc.emp <- readRDS(BASIS_FILE)
rot <- pc.emp$rotation
rot.DT <- data.table(pid=rownames(rot),rot)
M <- merge(M,rot.DT,by='pid')
M <- M[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
n1 <- 542
n0 <- 6717
M <- M[,est.seb:=1/sqrt(2 * ((n0*n1)/(n0+n1)) * 1 * maf * (1-maf))]
M[,an.seb:=abs(log(or)/qnorm(p.value/2,lower.tail=FALSE))]


foo<-fread("autosomes.egpa.high-info.gwas")
test <- foo[CHR==22,]
test[,est.seb:=1/sqrt(2 * ((controls_total*cases_total)/(controls_total+cases_total)) * info * cntr_maf * (1-cntr_maf))]
test[,an.seb:=abs(log(OR)/qnorm(P/2,lower.tail=FALSE))]
plot(log(test$est.seb),log(test$an.seb))
abline(a=0,b=1,col='red')

foo<-fread("autosomes.anca_Neg.high-info.gwas")
test <- foo[CHR==22,]
test[,est.seb:=1/sqrt(2 * ((controls_total*cases_total)/(controls_total+cases_total)) * 1 * cntr_maf * (1-cntr_maf))]
test[,an.seb:=abs(log(OR)/qnorm(P/2,lower.tail=FALSE))]
test[,calc.or:=(case_maf/(1-case_maf)) * ((1-cntr_maf)/cntr_maf)]
plot(log(test$OR),log(test$calc.or))
plot(test$est.seb,test$an.seb)
abline(a=0,b=1,col='red')



bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
pheno <- fread(bb_phenofile)
setnames(pheno,names(pheno) %>% make.names)
med <- pheno[grepl("20001\\_|20002\\_|20003\\_",Phenotype.Code) & Sex=='both_sexes',]
## load in phenotype file
P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/phenotypes.both_sexes.tsv')
P<-P[,.(phenotype,variable_type,non_missing=n_non_missing,cases=n_cases,controls=n_controls,pheno.source=source)]
## load in manifest
## compose a new command
med[,phe:=make.names(Phenotype.Description) %>% gsub("Cancer.code..self.reported..","SRC:",.)]
med[,phe:=gsub("Treatment.medication.code..","SRM:",phe)]
med[,phe:=gsub("Non.cancer.illness.code..self.reported..","SRD:",phe)]
phe.lu <- merge(med[,.(Phenotype.Code,phe)],P[,.(phenotype,cases,controls)],by.x='Phenotype.Code',by.y='phenotype')

SRC_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr_13_traits_0919'
neale.files <- list.files(path=SRC_DIR,pattern="^UKBB_NEALE\\:SRD.*",full.names=TRUE)

## we can randomly select a file
f <- sample(neale.files,1)
f <- neale.files[grep("myositis",neale.files)]
trait <- basename(f) %>% gsub("UKBB_NEALE:(.*)\\_source\\.RDS","\\1",.)
cc <- phe.lu[phe==trait,]

## compute variance estimates the normal way.
n <- cc$cases + cc$controls
n1 <- cc$cases
n0 <- cc$controls

can.DT <- readRDS(VARIANCE_FILE)
can.DT[,variance:=((log(n)-(log(n1) + log(n-n1)))+ log(mfactor)) %>% exp]

## next compute variances using the actual standard error of betas
sum.DT <- readRDS(f)
sum.DT[,seb:=abs(log(or)/qnorm(p.value/2,lower.tail=FALSE))]
if(FALSE){
  bb.snps.DT <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/variants.tsv')
  tmp <- bb.snps.DT[,.(pid=paste(chr,pos,sep=':'),varid,rsid,bb_ref=ref,bb_alt=alt,info)]
  M<-merge(sum.DT,tmp,by='pid',all.x=TRUE)
  saveRDS(M[,.(pid,info,varid,rsid,bb_ref,bb_alt)],file="~/tmp/ukbb_with_info_13_traits.RDS")
}
tmp <- readRDS("~/tmp/ukbb_with_info_13_traits.RDS")
M <- merge(sum.DT,tmp,by='pid',all.x=TRUE)
M <- M[,c('chr','pos'):=tstrsplit(pid,':')]
## need to add ld block
snp.DT <- fread(SNP_MANIFEST_FILE)
M <- merge(snp.DT,M,by='pid')
## finally add the PC rotations
pc.emp <- readRDS(BASIS_FILE)
rot <- pc.emp$rotation
rot.DT <- data.table(pid=rownames(rot),rot)
M <- merge(M,rot.DT,by='pid')
M <- M[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
M <- M[is.na(seb),seb:=1/sqrt(2 * ((n0*n1)/(n0+n1)) * 1 * maf * (1-maf))]


mvs_sigma<-function(r,quiet=TRUE){
  diag(r)<-1
  if(any(is.na(r)))
    if(!quiet)
      message(sprintf("Found %s where R^2 is NA",sum(is.na(r))))
    r[is.na(r)]<-0
  return(as(corpcor::make.positive.definite(r),"Matrix"))
}


compute_seb_var <- function(DT,ref_gt_dir,quiet=FALSE){
  s.DT <- split(DT,DT$chr)
  all.chr <- lapply(names(s.DT),function(chr){
    if(!quiet)
      message(sprintf("Processing %s",chr))
    ss.file<-file.path(REF_GT_DIR,sprintf("%s.RDS",chr))
    message(ss.file)
    sm <- readRDS(ss.file)
    ## there are sometimes duplicates that we need to remove
    pids <- colnames(sm)
    dup.idx<-which(duplicated(pids))
    if(length(dup.idx)>0){
      if(!quiet)
        message(sprintf("Warning removing %d duplicated SNPs",length(dup.idx)))
      sm <- sm[,-dup.idx]
      pids <- pids[-dup.idx]
    }
    # by ld block
    by.ld <- split(s.DT[[chr]],s.DT[[chr]]$ld.block)
    chr.var <- lapply(by.ld,function(block){
      if(!quiet)
        message(sprintf("Processing %s",block$ld.block %>% unique))
      sm.map <- match(block$pid,pids)
      if(any(is.na(sm.map))){
        message("SNPs in manifest that don't have genotypes")
      }
      r <- ld(sm[,sm.map],sm[,sm.map],stats="R")
      if(any(is.na(r)))
        if(!quiet)
          message(sprintf("Found %s where R^2 is NA",sum(is.na(r))))
      r[is.na(r)]<-0
      # compute closest pos-def covariance matrix
      #Sigma <- as.matrix(mvs_sigma(Matrix(r)))
      pc.cols <- which(grepl("^PC[0-9]+$",names(block)))
      sapply(names(block)[pc.cols],function(pc) (tcrossprod(block[,.(tot=get(`pc`) * ws_emp_shrinkage * seb)]$tot) * r) %>% sum)
    })
    s <- colSums(do.call('rbind',chr.var))
    message(s)
    s
  })
  colSums(do.call('rbind',all.chr))
}

seb.DT <- compute_seb_var(M,ref_gt_dir=REF_GT_DIR,quiet=FALSE)
seb.var <- data.table(pc=names(seb.DT),seb.var=seb.DT)
k <- merge(seb.var,can.DT,by='pc')
library(cowplot)
library(ggrepel)
ggplot(k[pc!='PC14',],aes(x=variance,y=seb.var,label=pc)) + geom_point() + geom_text_repel() + geom_abline(col='red',lty=2) +
xlab("Current variance estimation method") + ylab("Using actual SE(beta)")


## here we wish to compute dot product of sigma.n and sigma.maf
test <- M[,.(chr,pid,ld.block,seb,shrink=ws_emp_shrinkage,info,maf)]
## compute components of seb that are maf and sample size derived
test[,maf.se:=sqrt(1/(maf * (1-maf)))]
test[,n.se:=sqrt(1/(((n0*n1)/(n0+n1)) * info))]
test[,a.seb:=sqrt(1/2) * maf.se * n.se]
test[,sparse:=pid %in% keep.snps.dt$pid]

ggplot(test,aes(x=seb,y=a.seb)) + geom_point(size=0.1) + facet_wrap(~sparse) + geom_abline(col='red',lty=2) + ggtitle(trait)


plot(test$seb,test$a.seb)
abline(a=0,b=1,col='red',lty=2) ## looks properly scaled if noisy



## next do the same but for sparse basis


DRIVER_SNP_FILE <- '~/share/as_basis/sparse-basis/basis-sparse-13-0.999.RData'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
load(DRIVER_SNP_FILE)
## the only SNPs we care about are in use.pca
keep.snps.dt <- strsplit(rownames(rot.pca),':') %>% do.call('rbind',.) %>% data.table
keep.snps.dt <- data.table(pid=rownames(rot.pca))
keep.snps.dt[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]


compute_ssimp_var <- function(DT,N0,N1,ref_gt_dir,quiet=FALSE){
  imp.filt <- merge(DT,keep.snps.dt,by='pid',all.y=TRUE)
  ## standard error of the beta is computed like so
  imp.filt[is.na(seb),seb:=1/sqrt(2 * ((N0*N1)/(N0+N1)) * info * maf * (1-maf))]
  ## where seb is NA i.e. the SNP is missing from the dataset set its seb to 0 as it cannot contribute to the variance
  ## infinite happens if the snp cannot be imputed at all
  imp.filt[is.na(seb) | !is.finite(seb),seb:=0]
  rot.pca.dt <- data.table(pid=rownames(rot.pca),rot.pca)
  M <- merge(imp.filt[,.(pid,seb)],data.table(pid=rownames(rot.pca),rot.pca),by='pid')
  ## add in the shrinkage
  shrink.dt <- readRDS(SHRINKAGE_FILE)
  M <- merge(M,shrink.dt[,.(pid,shrink=ws_emp_shrinkage)],by='pid')
  by.chr <- split(M,gsub("^([^:]+):.*","\\1",M$pid))
  all.chr <- lapply(names(by.chr),function(chr){
    if(!quiet)
      message(sprintf("Processing %s",chr))
    ss.file<-file.path(ref_gt_dir,sprintf("%s.RDS",chr))
    #message(ss.file)
    sm <- readRDS(ss.file)
    ## there are sometimes duplicates that we need to remove
    pids <- colnames(sm)
    dup.idx<-which(duplicated(pids))
    if(length(dup.idx)>0){
      if(!quiet)
        message(sprintf("Warning removing %d duplicated SNPs",length(dup.idx)))
      sm <- sm[,-dup.idx]
      pids <- pids[-dup.idx]
    }
    M <- by.chr[[chr]]
    sm.map <- match(M$pid,pids)

    if(any(is.na(sm.map))){
      message("SNPs in manifest that don't have genotypes")
    }
    r <- ld(sm[,sm.map],sm[,sm.map],stats="R")
    if(any(is.na(r)))
      if(!quiet)
        message(sprintf("Found %s where R^2 is NA",sum(is.na(r))))
    r[is.na(r)]<-0
    pc.cols <- which(grepl("^PC[0-9]+$",names(M)))
    sapply(names(M)[pc.cols],function(pc) (tcrossprod(M[,.(tot=get(`pc`) * shrink * seb)]$tot) * r) %>% sum)
  })
  colSums(do.call('rbind',all.chr))
}


imp.DT <- M[,.(pid,seb,info,maf)]
sparse.DT <- compute_ssimp_var(imp.DT,N0=n0,N1=n1,ref_gt_dir=REF_GT_DIR,quiet=TRUE)
sparse.DT <- data.table(pc=names(sparse.DT),sparse.DT)

k <- merge(k,sparse.DT,by='pc')
ggplot(k[pc!='PC14',],aes(x=variance,y=sparse.DT,label=pc)) + geom_point() + geom_text_repel() + geom_abline(col='red',lty=2) +
xlab("Current variance estimation method") + ylab("Using sparse basis and SE(beta)")




library(snpStats)
library(corpcor)
library(Matrix)
library(mvtnorm)

compute_proj_var <- function(man.DT,w.DT,shrink.DT,ref_gt_dir,method='shrinkage_ws_emp',quiet=TRUE){
  M <- merge(man.DT,w.DT,by='pid')
  M <- merge(M,shrink.DT[,.(pid,shrink=get(`method`))],by='pid')
  M <- M[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
  ## create a set of analytical se_beta_maf
  beta_se_maf <- function(f) sqrt(1/f + 1/(1-f)) * sqrt(1/2)
  M <- M[,beta_se_maf:=beta_se_maf(maf)]
  M <- M[,c('chr','pos'):=tstrsplit(pid,':')]
  s.DT <- split(M,M$chr)
  all.chr <- lapply(names(s.DT),function(chr){
    if(!quiet)
      message(sprintf("Processing %s",chr))
    ss.file<-file.path(ref_gt_dir,sprintf("%s.RDS",chr))
    message(ss.file)
    sm <- readRDS(ss.file)
    ## there are sometimes duplicates that we need to remove
    pids <- colnames(sm)
    dup.idx<-which(duplicated(pids))
    if(length(dup.idx)>0){
      if(!quiet)
        message(sprintf("Warning removing %d duplicated SNPs",length(dup.idx)))
      sm <- sm[,-dup.idx]
      pids <- pids[-dup.idx]
    }
    # by ld block
    by.ld <- split(s.DT[[chr]],s.DT[[chr]]$ld.block)
    chr.var <- lapply(by.ld,function(block){
      if(!quiet)
        message(sprintf("Processing %s",block$ld.block %>% unique))
      sm.map <- match(block$pid,pids)
      if(any(is.na(sm.map))){
        message("SNPs in manifest that don't have genotypes")
      }
      r <- ld(sm[,sm.map],sm[,sm.map],stats="R")
      # compute closest pos-def covariance matrix
      Sigma <- as.matrix(mvs_sigma(Matrix(r)))
      pc.cols <- which(grepl("^PC[0-9]+$",names(M)))
      sapply(names(M)[pc.cols],function(pc) (tcrossprod(block[,.(tot=get(`pc`) * shrink * beta_se_maf)]$tot) * Sigma) %>% sum)
    })
    colSums(do.call('rbind',chr.var))
  })
  colSums(do.call('rbind',all.chr))
}




#' Code to sample multivariate norm
#' \code{mvs_perm} sample from a multivariate normal distribution
#'
#' @param m a vector -  mean values
#' @param sigma a matrix  - a positive definite covariance matrix
#' @param n a scalar - number of samples to generate (default 2)
#' @return a matrix of sampled values

mvs_perm<-function(m,sigma,n=2){
  if(!is.matrix(sigma))
    stop("sigma parameter is not a matrix")
  if(!corpcor::is.positive.definite(sigma,,method="chol"))
    stop("sigma is not positive definite")
  rd<-mvtnorm::rmvnorm(n,mean=m,sigma=sigma,method="chol")
  t(rd)
}

#' Code to compute sigma - genotype covariance matrix
#' \code{mvs_sigma} create a geneotype convariance matrix
#'
#' @param r a matrix of pairwise \eqn{r^2} measures generated by snpStats
#' @param quiet boolean whether to return warnings messages default=TRUE
#' @return a positive definite approximation of covariance matrix

mvs_sigma<-function(r,quiet=TRUE){
  diag(r)<-1
  if(any(is.na(r)))
    if(!quiet)
      message(sprintf("Found %s where R^2 is NA",sum(is.na(r))))
    r[is.na(r)]<-0
  return(as(corpcor::make.positive.definite(r),"Matrix"))
}

#' convert a vcf file to snpMatrix object
#' \code{vcf2snpmatrix} convert a vcf file to snpMatrix object
#'
#' @param vcf a scalar - path to vcf file
#' @param bcf_tools a scalar - path to bcftools binary
#' @param region_file a scalar - path to a file satisfying -R bcftools criteria (optional)
#' @param quiet a boolean - if set to false then debug information shown (default = FALSE)
#' @return a list with two slots. sm is the snpMatrix object and info is a data.table describes SNPs
#' @export

vcf2snpmatrix <- function(vcf,bcf_tools,region_file,quiet=TRUE){
  header_cmd <- sprintf("%s view -h %s",bcft,vcf.file)
  if(!quiet)
    message(header_cmd)
  my.pipe<-pipe(header_cmd)
  header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
  close(my.pipe)
  cnames<-unlist(strsplit(header,"\t"))
  if(!missing(region_file)){
    vcftools_cmd<-sprintf("%s view -R %s -O v %s | grep -v '^#'",bcft,region_file ,vcf.file)
  }else{
    vcftools_cmd<-sprintf("%s view -O v %s | grep -v '^#'",bcft,vcf.file)
  }
  if(!quiet)
    message(vcftools_cmd)
  tmp<-fread(vcftools_cmd)
  setnames(tmp,cnames)
  gt<-tmp[,10:ncol(tmp),with=FALSE]
  if(nrow(gt)==0)
    return(NA)
  info<-tmp[,1:9,with=FALSE]
  setnames(info,'#CHROM','CHROM')
  if(!quiet)
    message("Creating snpMatrix obj")
  sm<-apply(gt,1,function(x) sub("0\\|0","1",x))
  sm<-apply(sm,1,function(x) sub("(0\\|1)|(1\\|0)","2",x))
  sm<-apply(sm,1,function(x) sub("1\\|1","3",x))
  ## set anything else to a missing value
  sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\|[0-9]","0",x))))
  info[,pid:=paste(CHROM,POS,sep=':')]
  colnames(sm)<-info$pid
  rownames(sm)<-colnames(gt)
  return(list(sm=new("SnpMatrix", sm),info=info))
}

#' simulate betas for an ld block
#' \code{simulate_beta} use the multivariate normal to simulate realistic betas
#'
#' @param sm a snpMatrix object
#' @param lor a vector - observed betas
#' @param se_lor a vector - standard error of betas
#' @param lor_shrink a vector - shrinkage values to use to adjust betas (default 1 no shrinkage)
#' @param n_sims a scalar - number of simulations to conduct
#' @return a matrix of simulated betas

simulate_beta <- function(sm,lor,se_lor,lor_shrink=1,n_sims){
  beta_hat<-lor * lor_shrink
  if(length(lor)==1)
    return(t(rnorm(n_sims,mean=beta_hat,sd=se_lor)))
  # compute R statistic
  r<-snpStats::ld(sm,sm,stats="R")
  # compute closest pos-def covariance matrix
  r<-as.matrix(mvs_sigma(Matrix(r)))
  ## for beta the covariance matrix is estimates by sigma x SE * SE^T
  cov_se<-tcrossprod(se_lor)
  cov.beta<-cov_se * r
  ## simulate beta
  return(mvs_perm(beta_hat,cov.beta,n=n_sims))
}

#' covariance matrix of betas for an ld block
#' \code{cov_beta}
#'
#' @param sm a snpMatrix object
#' @param se_lor a vector - standard error of betas
#' @return a covariance matrix for beta

cov_beta <- function(sm,se_lor){
  #beta_hat <- lor * lor_shrink
  if(length(se_lor)==1)
    return(se_lor^2)
  # compute R statistic
  #r<-ld(sm,sm,stats="R.squared")
  r<-snpStats::ld(sm,sm,stats="R")
  # compute closest pos-def covariance matrix
  r<-as.matrix(mvs_sigma(Matrix(r)))
  ## for beta the covariance matrix is estimates by sigma x SE * SE^T
  cov_se<-tcrossprod(se_lor)
  return(cov_se * r)
}


#' simulate betas for a study
#' \code{simulate_study} use the multivariate normal to simulate realistic betas for a study
#'
#' @param DT a data.table - as returned by \code{\link{get_gwas_data}}
#' @param ref_gt_dir scalar - path to a dir of R objects named CHR_1kg.RData containing reference GT in snpMatrix format
#' @param shrink_beta scalar - Whether to use Bayesian shrinkage of betas (default = TRUE)
#' @param n_sims a scalar - number of simulations (default 10)
#' @param quiet a scalar - boolean whether to show progress messages
#' @return a DT of n_sims simulated studies for projection
#' @export

simulate_study <- function(DT,ref_gt_dir,shrink_beta=TRUE,n_sims=10,quiet=TRUE){
  s.DT <- split(DT,DT$chr)
  all.chr <- lapply(names(s.DT),function(chr){
    if(!quiet)
      message(sprintf("Processing %s",chr))
    ss.file<-file.path(ref_gt_dir,sprintf("%s.RDS",chr))
     message(file.path(ref_gt_dir,sprintf("%s.RDS",chr)))
    sm <- readRDS(ss.file)
    #sm<-get(load(ss.file))
    ## there are sometimes duplicates that we need to remove
    info <- data.table(pid=colnames(sm),order=1:ncol(sm))
    dup.idx <- which(duplicated(info$pid))
    #dup.idx<-which(duplicated(obj$info$pid))
    if(length(dup.idx)>0){
      if(!quiet)
        message(sprintf("Warning removing %d duplicated SNPs",length(dup.idx)))
      #sm$info<-sm$info[-dup.idx,]
      #sm$sm <- sm$sm[,-dup.idx]
      info <- info[-dup.idx,]
      sm <- sm[,-dup.idx]
    }
    info$order <- 1:nrow(info)
    #sm$info$order<-1:nrow(sm$info)
    # by ld block
    by.ld <- split(s.DT[[chr]],s.DT[[chr]]$ld.block)
    chr.sims <- lapply(names(by.ld),function(block){
      if(!quiet)
        message(sprintf("Processing %s",block))
      dat <- by.ld[[block]]
      setkey(dat,pid)
      #info <-sm$info[pid %in% dat$pid ,.(pid,order)]
      linfo <- info[pid %in% dat$pid ,.(pid,order)]
      setkey(linfo,pid)
      dat <- dat[linfo][order(order),]
      if(shrink_beta){
        ## compute beta shrinkage
        shrink <- with(dat,wakefield_null_pp(p.val,maf,n,n1/n))
        #M <- with(dat,simulate_beta(sm$sm[,order],log(or),emp_se,shrink,n_sims))
        M <- with(dat,simulate_beta(sm[,order],log(or),beta_se,shrink,n_sims))
      }else{
        #M <- with(dat,simulate_beta(sm$sm[,order],log(or),emp_se,1,n_sims))
        M <- with(dat,simulate_beta(sm[,order],log(or),beta_se,1,n_sims))
      }
      sims <- cbind(dat,M)
      sims <- melt(sims,id.vars=names(dat))
      return(sims[,c('or','trait'):=list(exp(value),paste(trait,variable,sep='_'))][,names(dat),with=FALSE])
    })
    chr.sims<-rbindlist(chr.sims)
  })
  all.chr<-rbindlist(all.chr)
  setkey(all.chr,pid)
}


#' analytically compute the variance of a projection given a reference set of genotypes
#' \code{compute_proj_var}
#'
#' @param man.DT a data.table - manifest of variants included in the basis
#' @param w.DT a data.table - a data table of weights/rotations from a  basis first column is 'pid' and subsequent columns for each principal component.
#' @param shrink.DT a data.table - as returned by \code{\link{compute_shrinkage_metrics}}
#' @param ref_gt_dir scalar - path to a dir of R objects named CHR_1kg.RData containing reference GT in snpMatrix format
#' @param method scalar - shrinkage method to use (default ws_emp)
#' @param quiet a scalar - boolean whether to show progress messages
#' @return a scalar of variances for each principal component
#' @export
