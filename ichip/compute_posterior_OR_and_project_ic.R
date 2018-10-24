library(data.table)
library(annotSnpStats)
library(magrittr)
library(cupcake)
library(optparse)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SNP_MANIFEST_FILE <- '/home/ob219/share/as_basis/ichip/snp_manifest/ichip_all.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/ichip/trait_manifest/as_manifest_ichip.tsv'
#CAL_FILE <- '/home/ob219/share/as_basis/ichip/support/por_2500_2.0_0.01.RDS'
CAL_FILE <- '/home/ob219/share/as_basis/ichip/support/por_2500_2.0_0.01.RDS'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/ichip/support/shrinkage_ic.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/ichip/support/basis_ic.RDS'
POR_OUT_DIR <- '/home/ob219/share/as_basis/ichip/individual_data/individual_por/'
PROJ_OUT_DIR <- '/home/ob219/share/as_basis/ichip/individual_data/individual_proj/'
ICHIP_DATA_DIR <- '/home/ob219/share/as_basis/ichip/sum_stats'


SAMP_CHUNK_SIZE <- 100

TEST <- FALSE

option_list = list(
  make_option(c("-d", "--data_dir"), type="character",default='',
              help="Location of processed snpStats objects", metavar="character"))

if(TEST){
  args <- list(
    data_dir='/home/ob219/share/as_basis/ichip/individual_data/filtered_gt/ERA'
  )
}else{
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}

print(args)

fls <- list.files(path=args$data_dir,pattern="*.RDS",full.names=TRUE)
cal <- readRDS(CAL_FILE)[,2:4,with=FALSE] %>% as.matrix
man<-fread(SNP_MANIFEST_FILE)

get_por <- function(f){
  message(sprintf("Processing %s",f))
  X <- readRDS(f)
  snps <- snps(X) %>% data.table
  snps <- snps[,.(pid,uid=1:.N)]
  m <- merge(man,snps,by.x='pid',by.y='pid')[order(uid),]
  #if(sum(1:nrow(m) != m$uid)!=0)
  #  stop("Manifest does not match genotyped SNPs - please check !")
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

## create basis for projection
#basis.DT<-get_gwas_data(MANIFEST_FILE,SNP_MANIFEST_FILE,ICHIP_DATA_DIR,filter_snps_by_manifest=TRUE)
#shrink.DT<-compute_shrinkage_metrics(basis.DT)
## need to add control where beta is zero
#basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
#basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
#pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)

## process proj.lors
DT.lor <- data.table(obj$proj.lor)
setnames(DT.lor,rownames(obj$samples))
DT.lor[,pid:=obj$snps$pid]
ind.proj.DT <- melt(DT.lor,id.vars='pid')[,or:=exp(value)]
setnames(ind.proj.DT,'variable','trait')
by.samp <- split(ind.proj.DT,ind.proj.DT$trait)
samp.no <- length(names(by.samp))
group.size <- ceiling(samp.no/SAMP_CHUNK_SIZE)
if(group.size==1){
  samp.groups <- list(1:length(names(by.samp)))
}else{
  samp.groups <- split(names(by.samp),cut(1:samp.no,group.size))
}

all.proj <- lapply(samp.groups,function(ss){
  message(length(ss))
  sample.proj.DT <- by.samp[ss] %>% rbindlist
  sample.proj.DT[,trait:=as.character(trait)]
  setkey(sample.proj.DT,pid)
  mat.emp <- create_ds_matrix(sample.proj.DT,shrink.DT,SHRINKAGE_METHOD)
  if(!identical(colnames(mat.emp), rownames(pc.emp$rotation)))
    stop("Something wrong basis and projection matrix don't match")
  predict(pc.emp,newdata=mat.emp)
}) %>% do.call('rbind',.)

fname <- file.path(PROJ_OUT_DIR,sprintf("%s_projection.RDS",trait))
saveRDS(all.proj,file=fname)

if(FALSE){
  R_SCRIPT <- '/home/ob219/git/basis_paper/ichip/compute_posterior_OR_and_project_ic.R'
  dirs <- list.files(path='/home/ob219/share/as_basis/ichip/individual_data/filtered_gt',pattern="*",full.names=TRUE)
  cmds <- lapply(dirs,function(d){
      sprintf("Rscript %s -d  %s",R_SCRIPT,d)
  }) %>% do.call('c',.)
  write(cmds,"~/tmp/qstuff/ind_proj.txt")
}

## example single use
#Rscript   /home/ob219/git/as_basis/R/Individual_projection/compute_posterior_OR_and_project.R -d /home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/individual_gt/gtex

if(FALSE){
  library(ggplot2)
  library(ggrepel)
  ind.DT <- all.proj %>% data.table
  ind.DT[,trait:='']
  basis.DT <- data.table(pc.emp$x,trait=rownames(pc.emp$x))
  ggplot(rbind(ind.DT,basis.DT),aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text_repel()
}
