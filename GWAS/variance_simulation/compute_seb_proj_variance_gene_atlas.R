## compute the analytical variances associated with projecting onto a basis

library(devtools)
load_all("~/git/cupcake")
library(optparse)

TEST<-FALSE
option_list = list(
        make_option(c("-f", "--files"), type="character", default=NULL,
              help="file containing list of files to process", metavar="character")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$files)){
	   print_help(opt_parser)
	    stop("Supply a file containing a list of files to process", call.=FALSE)
    }
}else{
  args <- list(files = '/home/ob219/tmp/qstuff/avar/file891822fef6df')
}




REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_gwas_13_traits_0919.RDS'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'
SRC_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr_13_traits_0919'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/seb_proj_var_13_traits_0919'
meta.dt <- fread("~/tmp/41588_2018_248_MOESM3_ESM.csv")
meta.dt <- meta.dt[Category=='Binary',.(ID,Description,Cases,Controls=round(Cases/Sample)-Cases,prop=Sample)]
meta.dt <- meta.dt[,phe:=make.names(Description)]
convertORscale <- function(x,cp) x/(cp * (1-cp))


if(FALSE){
  CHUNK_SIZE <- 10
  lout <- '/home/ob219/tmp/qstuff/avar/'
  all.files <- list.files(path=SRC_DIR,pattern="^GA*",full.names=TRUE)
  foo <- lapply(split(all.files,ceiling(seq_along(all.files)/CHUNK_SIZE)),function(fs){
    tfile <- tempfile(tmpdir=lout)
    write(fs,file=tfile)
    tfile
  })
  rscript <- '/home/ob219/git/basis_paper/GWAS/variance_simulation/compute_seb_proj_variance_gene_atlas.R'
  cmds <- sapply(foo,function(x){
    sprintf("Rscript %s -p %s",rscript,x)
  })
  write(cmds,file="~/tmp/qstuff/variance.txt")
}



files <- scan(args$files,"character")

for(f in files){
  trait <- basename(f) %>% gsub("^(.*)\\_source\\.RDS","\\1",.)
  prop <- meta.dt[phe==gsub("GA:","",trait),]$prop
  message(trait)
  ofile <- file.path(OUT_DIR,sprintf("%s.RDS",trait))
  if(file.exists(ofile)){
    sprintf("Already done %s. skipping",trait)
    next
  }
  ## these files contain the shrinkages so no need to add those at this stage
  gwas.DT <- readRDS(f)
  ## fit to the or scale
  gwas.DT[,tmp.or:=convertORscale(log(or),prop) %>% exp]
  gwas.DT[,or:=tmp.or]
  shrink.DT <- readRDS(SHRINKAGE_FILE)
  setnames(shrink.DT,'ws_emp_shrinkage','shrinkage')
  man.DT <- fread(SNP_MANIFEST_FILE)
  pc.emp <- readRDS(BASIS_FILE)
  rot <- pc.emp$rotation
  w.DT <- data.table(pid=rownames(rot),rot)
  trait <- basename(f) %>% gsub("^(.*)\\_source\\.RDS","\\1",.)
  test <- compute_seb_proj_var(gwas.DT,shrink.DT,man.DT,w.DT,ref_gt_dir=REF_GT_DIR,quiet=TRUE)
  saveRDS(test,file=ofile)
}
