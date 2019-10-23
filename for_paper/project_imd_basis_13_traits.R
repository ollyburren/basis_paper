## once we have aligned and processed files
## it is much easier to add traits to the basis as
## the fdr source files contain the harmonised beta/or values
## that are required for projection - all we then need to do
## is assign the weights and hey presto.

library(devtools)
load_all("~/git/cupcake")
library(parallel)
#BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
ROOT.DIR <- '~/share/as_basis/basis-creation'
SNP_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_snp_manifest.RDS')
SHRINKAGE_FILE <- file.path(ROOT.DIR,'support/13_trait_shrinkage.RDS')
BASIS_FILE <- file.path(ROOT.DIR,'support/13_trait_basis.RDS')
NOWEIGHT_BASIS_FILE <- file.path(ROOT.DIR,'support/13_trait_basis-noweight.RDS')
DATA_DIR <- file.path(ROOT.DIR,'../basis-projection/project_13')
OUT_DIR <- file.path(ROOT.DIR,'../basis-projection/project_13_new')

pc.emp <- readRDS(BASIS_FILE)
man.DT <- readRDS(SNP_MANIFEST_FILE)
sDT <- readRDS(SHRINKAGE_FILE)
#stmp<-sDT[,.(pid,ws_emp_shrinkage)]

all.filez <- list.files(path=DATA_DIR,pattern='*.RDS',full.names=TRUE)
proj <- mclapply(all.filez,function(f){
  trait <- basename(f) %>% gsub("\\_source.RDS","",.)
  sprintf("Processing %s",trait) %>% message
  dat <- readRDS(f)
  x<-project_basis(dat,sDT,pc.emp,trait)
  saveRDS(x$data,file=file.path(OUT_DIR,sprintf("%s_source.RDS",trait)))
  saveRDS(x$missing,file=file.path(OUT_DIR,sprintf("%s_missing.RDS",trait)))
  x$proj
},mc.cores=8)

proj <- do.call('rbind',proj)
saveRDS(proj,"~/share/as_basis/basis-projection/results/13_traits.RDS")
delta <- proj -  pc.emp$x["control",]

if(FALSE){
  #check how different they look
  all.results <- readRDS('/home/ob219/share/as_basis/GWAS/RESULTS/non_uk_bb_for_fdr_13_traits_0919_v2.RDS')
  delta.old <- all.results - readRDS('/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS')$x["control",]
  old <- data.table(trait=rownames(delta.old),all.results) %>% melt(.,id.vars='trait')
  new <- data.table(trait=rownames(delta),proj) %>% melt(.,id.vars='trait')
  M <- merge(old[,.(trait,variable,old.delta=value)],new[,.(trait,variable,new.delta=value)],by=c('trait','variable'),all=TRUE)
  ## possible supplementary figure
  par(mfrow=c(4,4))
  for(pc in paste0('PC',1:13)){
    plot(M[variable==pc,]$old.delta,M[variable==pc,]$new.delta,xlab="old delta",ylab="new delta",main=pc)
    abline(0,1,col='red',lty=2)
    abline(0,-1,col='red',lty=2)
  }
}
