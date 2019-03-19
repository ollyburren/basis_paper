#/# This script simulates what happens under the null i.e when the scaled beta's i.e. gamma hats are equal to one another.

## assumes that cupcake has been loaded with (devtools)
library(devtools)
load_all("~/git/cupcake")
#library(cupcake)
library(optparse)

TEST <- FALSE # set to true when debugging

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
DEFAULT.N.SIMS <- 50
DEFAULT.TRAIT <- 'ret_p'






#DEFAULT.SUPPORT.DIR <- '/home/ob219/rds/hpc-work/as_basis/support_tab'
#DEFAULT.SNPSTATS.DIR <- '/home/ob219/rds/hpc-work/as_basis/snpStats/basis_1kg_feb/'

#DEFAULT.GWAS.DATA.DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned/'
#DEFAULT.MANIFEST.FILE <- 'as_basis_manifest_with_jia_cc.tab'
#DEFAULT.MANIFEST.FILE <- '/home/ob219/git/as_basis/manifest/as_manifest_feb_2018_w_ms.csv'
#DEFAULT.REF.AF.FILE <- '/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab';
## does not matter what this is as long as has correct number of SNPs !
#DEFAULT.TRAIT <- 'ret_p'
#DEFAULT_CACHE_BASIS <- '/home/ob219/rds/hpc-work/as_basis/cache/as_basis_snp_support_feb_2018_w_ms.RDS'
#DEFAULT_CACHE_BASIS <- '/home/ob219/rds/hpc-work/as_basis/cache/as_basis_snp_support_feb_2018_w_ms_ws.RDS'
#SHRINKAGE_METHOD <- 'emp'
#SHRINKAGE_METHOD <- 'ws_emp'

## this routine takes a set of vcf files from a reference panel say 1KG and extracts the variants that are used in our basis.
#In this form it won't work (needs basis info for filtering but acts as a reference)

## to prevent unneccsary computation we can compute the basis ahead of time and just read that in

#cacheBasis <- function(){
#  basis.DT<-get_gwas_data(DEFAULT.MANIFEST.FILE,DEFAULT.REF.AF.FILE,DEFAULT.GWAS.DATA.DIR,FALSE)
#  shrink.DT<-compute_shrinkage_metrics(basis.DT)
#  basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
#  ## need to add control where beta is zero
#  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
#  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
#  saveRDS(list(basis=pc.emp,shrink=shrink.DT),file=DEFAULT_CACHE_BASIS)
#}

## perhaps add options for different manifests etc.
option_list = list(
        make_option(c("-n","--total"), type="numeric",default=DEFAULT.SUPPORT.DIR,
                    help="Location of support files", metavar="character"),
        make_option(c("-c", "--ncases"), type="numeric", default=DEFAULT.GWAS.DATA.DIR,
                    help="location of OR aligned GWAS source files", metavar="character"),
        make_option(c("-o", "--out_dir"), type="character", default=NULL,
                    help="output directory", metavar="character")
        )


if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  print(args)
}else{
  message("IN TESTING MODE ======>!")
  trait <- 'ret_p'
  args <- list(
    total = 10000,
    ncases=5000,
    out_dir = '~/tmp/test_distance/'
  )
  print(args)
}

est_SE<-function(f,n.cases,total){
  n.ctrls<-total-n.cases
  a<-1/(2*f*n.ctrls)
  b<-1/(2*(1-f)*n.ctrls)
  c<-1/(2*f*n.cases)
  d<-1/(2*(1-f)*n.cases)
  sqrt(rowSums(cbind(a,b,c,d)))
}

if(FALSE){
  OUTDIR <- '/home/ob219/share/as_basis/GWAS/variance_simulations'
  cmd_template<-"Rscript /home/ob219/git/basis_paper/GWAS/variance_simulation/null_by_sample_size_special.R -n %d -c %d -o %s"
  all.cmds <- lapply(c(1000,2000,3000,15000,5000),function(n){
    sapply(c(250,500,1000,5000,10000,50000),function(n1){
      if(n1<n)
        rep(sprintf(cmd_template,n,n1,OUTDIR),40)
      })
    }) %>% unlist
    write(all.cmds,"~/tmp/qstuff/run_cal.txt")
}


study1.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE,DEFAULT.TRAIT)
shrink.DT <- readRDS(SHRINKAGE_FILE)
M <- merge(study1.DT,shrink.DT,by='pid')
pc.emp <- readRDS(BASIS_FILE)

# check the above code by computing the SE looks ok but will be completely out for high p.values
test <- M[,.(pid,beta=log(or),se.beta=abs(log(or))/qnorm(p.value/2,lower.tail=FALSE),n,n1,maf,p.value)]
test[,est.se.beta:=est_SE(maf,n1,n)]
## all the noise comes from things that are not associated at all
## probably because qnorm pdf has a lot of area under it therefore big error with back calc of the standard error
keep <- test$p.value>1e-1
plot(test[keep,]$se.beta,test[keep,]$est.se.beta)
## do need to be careful of the instances where log(or) = 0 as these are
## where se will be completely wrong.
study1.DT$n <- args$total
study1.DT$n1 <- args$ncases

study1.DT[,beta_se:=est_SE(maf,n1,n)]
#study1.DT[,emp_se_old := 1/sqrt(2) * 1/sqrt(n) * emp_maf_se]
study1.DT[,or:=1]
study1.DT[,c('chr','position'):=tstrsplit(pid,':')]
message(sprintf("Simulate %d studies",args$n_sims))
s1.sim <- simulate_study(study1.DT,REF_GT_DIR,shrink_beta=FALSE,n_sims=DEFAULT.N.SIMS,quiet=FALSE)
idx <- split(1:nrow(s1.sim),s1.sim$trait)
pred.emp<-do.call('rbind',lapply(split(idx, ceiling(seq_along(idx)/10)),function(i){
  tmp.idx <- unlist(i)
  tmp.sim <- s1.sim[tmp.idx,]
  setkey(tmp.sim,pid)
  sim.mat.emp<-create_ds_matrix(tmp.sim,shrink.DT,SHRINKAGE_METHOD)
  predict(pc.emp,newdata=sim.mat.emp)
}))
rand.string <- tolower(do.call(paste0, replicate(6, sample(LETTERS, 1, TRUE), FALSE)))
of <- sprintf("%d_%d_%s.RDS",args$total,args$ncases,rand.string)
ofile <- file.path(args$out_dir,of)
saveRDS(pred.emp,file=ofile)
message(sprintf("Saved projections in %s",ofile))
