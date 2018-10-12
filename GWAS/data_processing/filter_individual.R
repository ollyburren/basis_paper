library(annotSnpStats)
library(optparse)
library(magrittr)
library(data.table)

TEST <- FALSE

option_list = list(
  make_option(c("-f", "--file"), type="character",default='',
              help="File to process", metavar="character"),
  make_option(c("-m", "--manifest_file"), type="character",default='',
              help="manifest file to filter on", metavar="character"),
  make_option(c("-o", "--outdir"), type="character",default='',
              help="Output directory", metavar="character"),
  make_option(c("-s", "--sample_filter"), type="character",default='',
              help="sample_filter", metavar="character"),
  make_option(c("-x", "--chr_col"), type="character",default='chr',
              help="chr column", metavar="character"),
  make_option(c("-y", "--position_col"), type="character",default='position',
              help="position column", metavar="character")
)


if(TEST){
  args <- list(
    manifest_file='/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab',
    file = "/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/JIA-2017-data/as_basis/chr6.RData",
    outdir = "/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/individual_gt//jiasys",
    sample_filter='alt_ilar_code:sys',
    chr_col='chromosome',
    position_col='position'
)
}else{
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}

print(args)

X <- get(load(args$file))
sample.idx <- vector()
if(args$sample_filter!=''){
  filters <- strsplit(args$sample_filter,':') %>% unlist
  sample.idx <- which(samples(X)[[filters[1]]]==filters[2])
  if(length(sample.idx)==0)
    stop("Filter drew a blank")
}


chr<-basename(args$file) %>% gsub("chr([0-9]+).*\\.RData","\\1",.)
man <- fread(args$manifest_file)[grep(sprintf("^%s:",chr),pid),]
pids <- paste(snps(X)[[args$chr_col]],snps(X)[[args$position_col]],sep=':')
#pids <- with(snps(X),paste(`args$chr_col`,`args$position_col`,sep=':'))
snps.idx <- which(pids %in% man$pid)


if(length(snps.idx)==0)
  stop("No SNPs found")
if(length(snps.idx)!=nrow(man))
  stop("SNP lists not equal")

if(length(sample.idx)>0){
  Xf <- X[sample.idx,snps.idx]
}else{
  Xf <- X[,snps.idx]
}
ofile <- file.path(args$outdir,basename(args$file)) %>% gsub("RData$","RDS",.)
repids <- paste(snps(Xf)[[args$chr_col]],snps(Xf)[[args$position_col]],sep=':')
snps(Xf)[['snp.name']] <- repids
saveRDS(Xf,file=ofile)

### We wish to have a single unified manifest file that we use for everything
## at the moment the PBC individual data is too sparse (needs imputation) so
## we exclude

if(FALSE){
  ## run JIA
  library(annotSnpStats)
  DATA.DIR <- '/home/ob219/share/as_basis/GWAS/individual_data/gt/jia/'
  OUT_DIR <- '//home/ob219/share/as_basis/GWAS/individual_data/filtered_gt'
  RSCRIPT <- '/home/ob219/git/basis_paper/GWAS/data_processing/filter_individual.R'
  SNP_MANIFEST_FILE <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
  X<-get(load(file.path(DATA.DIR,'chr22.RData')))
  all.files <- list.files(path=DATA.DIR,pattern="*.RData",full.names=TRUE)
  ilar <- unique(samples(X)$alt_ilar_code) %>% paste('alt_ilar_code',.,sep=':')
  # create dirs
  dirs <- unique(samples(X)$alt_ilar_code) %>% paste0('jia',.) %>% file.path(OUT_DIR,.)
  for(d in dirs)
    dir.create(d,showWarnings=FALSE)

  cmds <- lapply(ilar,function(i){
    lapply(all.files,function(f){
      dname <- gsub("alt_ilar_code:","",i) %>% paste0('jia',.) %>% file.path(OUT_DIR,.)
      sprintf("Rscript %s --file %s --manifest_file %s --sample_filter %s --outdir %s --chr_col chromosome",RSCRIPT,f,SNP_MANIFEST_FILE,i,dname)
    }) %>% do.call('c',.)
  }) %>% do.call('c',.)
  write(cmds,file="~/tmp/qsub/filter_jia.txt")

  ## gtex

  DATA.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Data/gtex/eur-Whole_Blood-001/as_basis/'
  OUT_DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/individual_gt/gtex'
  RSCRIPT <- '/home/ob219/git/as_basis/R/Individual_projection/filter_individual.R'
  MANIFEST_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab'
  all.files <- list.files(path=DATA.DIR,pattern="*.RData",full.names=TRUE)
  cmds <- lapply(all.files,function(f){
    sprintf("Rscript %s --file %s --manifest_file %s --outdir %s",RSCRIPT,f,MANIFEST_FILE,OUT_DIR)
  }) %>% do.call('c',.)
  write(cmds,file="~/tmp/qsub/filter_gtex.txt")


}
