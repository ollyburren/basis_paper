library(data.table)
library(optparse)

TEST <- FALSE

option_list = list(
  make_option(c("-f", "--file"), type="character",default='',
              help="File to process", metavar="character"),
  make_option(c("-m", "--manifest_file"), type="character",default='',
              help="manifest file to filter on", metavar="character"),
  make_option(c("-af", "--allele_flip"), action="store_true",type="logical",default=FALSE,
              help="whether to flip OR", metavar="character"),
  make_option(c("-o", "--outdir"), type="character",default='',
              help="Output directory", metavar="character")
)


if(TEST){
  args <- list(
    manifest_file='/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snp_support_uk10k.tab',
    file = "/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/astle_lymph_p.tab",
    outdir = "/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/all_studies_filtered",
    allele_flip=FALSE
)
}else{
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}

print(args)



man.DT <- fread(args$manifest_file)
DT <- fread(args$file)
out.DT <- DT[pid %in% man.DT$pid,.(pid,a1,a2,or,p.value)]
if(args$allele_flip)
  out.DT[,or:=1/or]
fname <- file.path(args$outdir,basename(args$file))
write.table(out.DT,file=fname,quote=FALSE,row.names=FALSE,sep="\t")

if(FALSE){
  library(magrittr)
  library(data.table)
  MANIFEST_FILE <- '/home/ob219/git/as_basis/manifest/as_manifest_july.tsv'
  BASE_DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/'
  RSCRIPT <- '/home/ob219/git/as_basis/R/Individual_projection/filter_summary_stats_q.R'
  OUTDIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/all_studies_filtered'
  SNP_MANIFEST_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab'
  m.DT <- fread(MANIFEST_FILE)[genotype=='N' & include=='Y',]
  cmds <- lapply(m.DT$file,function(f){
    fs <- file.path(BASE_DIR,f)
    effect_allele <- m.DT[file==f,]$effect_allele
    template <- "Rscript %s --file %s --manifest_file %s --outdir %s"
    if(effect_allele=='a1')
      template <- "Rscript %s --file %s --manifest_file %s --outdir %s --allele_flip"
    sprintf(template,RSCRIPT,fs,SNP_MANIFEST_FILE,OUTDIR)
  }) %>% do.call('c',.)
  write(cmds,file="~/tmp/qsub/filter_summ.txt")
}
