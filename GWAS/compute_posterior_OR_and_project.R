library(data.table)
library(annotSnpStats)
library(magrittr)
library(cupcake)
library(optparse)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
CAL_FILE <- '/home/ob219/rds/hpc-work/as_basis/support//por_2500_2.0_0.01.RDS'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
POR_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/individual_data/individual_por/'
PROJ_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/individual_data/individual_proj/'
SAMP_CHUNK_SIZE <- 50

TEST <- FALSE

option_list = list(
  make_option(c("-d", "--data_dir"), type="character",default='',
              help="Location of processed snpStats objects", metavar="character"))

if(TEST){
  args <- list(
    data_dir='/home/ob219/share/as_basis/GWAS/individual_data/filtered_gt/jiasys'
  )
}else{
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}

if(FALSE){
  DATA.DIR='/home/ob219/share/as_basis/GWAS/individual_data/filtered_gt/'
  RSCRIPT <- '/home/ob219/git/basis_paper/GWAS/compute_posterior_OR_and_project.R'
  all.files <- list.files(path=DATA.DIR,pattern="*",full.names=TRUE)
  cmd <- lapply(all.files,function(f){
    sprintf("Rscript %s --data_dir %s",RSCRIPT,f)
  }) %>% do.call('c',.)
  write(cmd,file="~/tmp/qstuff/por_jia.txt")
}

print(args)

fls <- list.files(path=args$data_dir,pattern="*.RDS",full.names=TRUE)
cal <- readRDS(CAL_FILE)[,2:4,with=FALSE] %>% as.matrix
man<-fread(SNP_MANIFEST_FILE)

get_por <- function(f){
  message(sprintf("Processing %s",f))
  X <- readRDS(f)
  snps <- snps(X) %>% data.table
  snps <- snps[,.(snp.name,uid=1:.N)]
  m <- merge(man,snps,by.x='pid',by.y='snp.name')[order(uid),]
  if(sum(1:nrow(m) != m$uid)!=0)
    stop("Manifest does not match genotyped SNPs - please check !")
  sm <- as(X,'SnpMatrix')
  sm<-matrix(sm,nrow=nrow(sm),ncol=ncol(sm))
  gt.lor <- lapply(seq_along(m$ref_a1.af),function(j){
    i <- round((m$ref_a1.af[j])*1000)
    as.vector(pp(sm[,j]) %*% cal[i,])
  }) %>% do.call("rbind",.)
  list(snps=m[,.(pid,ref_a1,ref_a2,ref_a1.af,ld.block)],proj.lor=gt.lor,samples=samples(X))
}

all <- lapply(fls,get_por)
## combine into a single file for the genome
snps <- lapply(all,'[[','snps') %>% rbindlist
proj.lor <- lapply(all,'[[','proj.lor') %>% do.call('rbind',.)
samples <- all[[1]]$samples
obj <- list(snps=snps,proj.lor=proj.lor,samples=samples)
trait <- basename(args$data_dir)
fname <- file.path(POR_OUT_DIR,sprintf("%s_por.RDS",trait))
saveRDS(obj,file=fname)
pc.emp <- readRDS(BASIS_FILE)
shrink.DT <- readRDS(SHRINKAGE_FILE)

# ## create basis for projection
# basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=FALSE)
# ## these SNPs are PITA
# basis.DT[p.value==0,p.value:=0.99]
# shrink.DT<-compute_shrinkage_metrics(basis.DT)
# ## need to add control where beta is zero
# basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
# ## need to add control where beta is zero
# basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
# pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)

## process proj.lors
DT.lor <- data.table(obj$proj.lor)
setnames(DT.lor,rownames(obj$samples))
DT.lor[,pid:=obj$snps$pid]
ind.proj.DT <- melt(DT.lor,id.vars='pid')[,or:=exp(value)]
setnames(ind.proj.DT,'variable','trait')

by.samp <- split(ind.proj.DT,ind.proj.DT$trait)
samp.no <- length(names(by.samp))
group.size <- ceiling(samp.no/SAMP_CHUNK_SIZE)
samp.groups <- split(names(by.samp),cut(1:samp.no,group.size))

all.proj <- lapply(samp.groups,function(ss){
  message(length(ss))
  sample.proj.DT <- by.samp[ss] %>% rbindlist
  sample.proj.DT[,trait:=as.character(trait)]
  setkey(sample.proj.DT,pid)
  mat.emp <- create_ds_matrix(sample.proj.DT,shrink.DT,SHRINKAGE_METHOD)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
    stop("Something wrong basis and projection matrix don't match")
  predict(pc.emp,newdata=mat.emp)
}) %>% do.call('rbind',.)

fname <- file.path(PROJ_OUT_DIR,sprintf("%s_projection.RDS",trait))
saveRDS(all.proj,file=fname)

## example single use
#Rscript   /home/ob219/git/as_basis/R/Individual_projection/compute_posterior_OR_and_project.R -d /home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/individual_gt/gtex

if(FALSE){
  ### do some checking if we have the alleles correct then we should get correlation between mean OR for inds and
  ### those derived from summary statistics
  sum.stats <- fread("/home/ob219/share/as_basis/GWAS/sum_stats/jiasys_unpub.tab")
  as <- data.table(pid=obj$snps$pid,i.a1=obj$snps$ref_a1,i.a2=obj$snps$ref_a2,mean.por=rowMeans(obj$proj.lor),ref_a1.af=obj$snps$ref_a1.af)
  cm <- merge(sum.stats[,.(pid,a1,a2,lor=log(or))],as,by.x='pid',by.y='pid')
  library(ggplot2)
  library(cowplot)
  cm[,af_cat:=cut(ref_a1.af,seq(0,1,length.out=10))]
  ggplot(cm[-grep("^6:",pid),],aes(x=lor,y=mean.por)) + geom_point() +
  facet_wrap(~af_cat) + xlab("Summary Stats log(OR)") +
  ylab("mean(Individual posterior log(OR))") + ggtitle("JIA SYS")

  ## check allele freq
  X <- readRDS(fls[1])
  cs<-col.summary(X) %>% data.table
  cs[,pid:= snps(X)$snp.name]
  M <- merge(man[,.(pid,ref_a1.af)],cs[,.(pid,RAF)],by='pid')

}
