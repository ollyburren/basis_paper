#/# This script simulates what happens under the null i.e when the scaled beta's i.e. gamma hats are equal to one another.

## assumes that cupcake has been loaded with (devtools)
library(devtools)
load_all("~/git/cupcake")
#library(cupcake)
library(optparse)

TEST <- FALSE # set to true when debugging

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
#SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
SHRINKAGE_FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/shrink_ws.RDS"
#BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
BASIS_FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/basis_ws.RDS"
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
#SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
SNP_MANIFEST_FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/support/snp_support_bootstrap_USE.tab"
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
DEFAULT.N.SIMS <- 50
DEFAULT.TRAIT <- 'ret_p'

option_list = list(
        make_option(c("-n","--total"), type="numeric",default=1000,
                    help="Location of support files", metavar="character"),
        make_option(c("-c", "--ncases"), type="numeric", default=500,
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


#study1.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE,DEFAULT.TRAIT)
snps.DT <- fread(SNP_MANIFEST_FILE)
shrink.DT <- readRDS(SHRINKAGE_FILE)
study1.DT <- copy(shrink.DT)
study1.DT <- merge(study1.DT,snps.DT[,.(pid,maf=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af))],by='pid')
pc.emp <- readRDS(BASIS_FILE)
study1.DT$n <- args$total
study1.DT$n1 <- args$ncases


study1.DT[,beta_se:=est_SE(maf,n1,n)]
#study1.DT[,emp_se_old := 1/sqrt(2) * 1/sqrt(n) * emp_maf_se]
study1.DT[,or:=1]
study1.DT[,c('chr','position'):=tstrsplit(pid,':')]
study1.DT[,trait:='test']
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


## this code prints out the commands for slurm job submission
if(FALSE){
  OUTDIR <- '/home/ob219/share/as_basis/GWAS/variance_simulations_bootstrap_basis'
  cmd_template<-"Rscript /home/ob219/git/basis_paper/GWAS/variance_simulation/null_using_bootstrap_basis.R -n %d -c %d -o %s"
  all.cmds <- lapply(c(1000,2000,3000,15000,5000),function(n){
    sapply(c(250,500,1000,5000,10000,50000),function(n1){
      if(n1<n)
        rep(sprintf(cmd_template,n,n1,OUTDIR),40)
      })
    }) %>% unlist
    write(all.cmds,"~/tmp/qstuff/run_sim_bs.txt")
}
