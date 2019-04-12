## work out the list of variants reqiured

library(optparse)

TEST<-FALSE
option_list = list(
        make_option(c("-f", "--fname"), type="character", default=NULL,
              help="index of phenotype to process ", metavar="character")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$fname)){
	   print_help(opt_parser)
	    stop("Supply an file containing a list of genes to process", call.=FALSE)
    }
}else{
  args <- list(fname='/home/ob219/tmp/qstuff/eqtlgen/genelists/gene_90.txt')
}

gene.list <- scan(args$fname,"character")

#gene.list <- c('ENSG00000270170','ENSG00000270172','ENSG00000270175','ENSG00000270177','ENSG00000270179','ENSG00000270184')

## to use maf estimate of se remove ss prefix !

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
## just the one shrinkage file
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
DATA_DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/sum_stats/eqtlgen_notrans/'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
OUT_DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/eqtlgen_projections_no_trans/'


## running on the queue
if(FALSE){
  files <- list.files(path=DATA_DIR,pattern="*.tab",full.names=FALSE)
  ensgs <- basename(files) %>% gsub("\\.tab","",.)
  CHUNK_SIZE <- 100
  gchunk <- split(ensgs, ceiling(seq_along(ensgs)/CHUNK_SIZE))
  for(i in seq_along(gchunk)){
    fname <- file.path('/home/ob219/tmp/qstuff/eqtlgen/genelists',sprintf("gene_%d.txt",i))
    write(gchunk[[i]],file=fname)
  }
  files <- list.files(path='/home/ob219/tmp/qstuff/eqtlgen/genelists',pattern="*.txt",full.names=TRUE)
  cmds <- sapply(files,function(f){
  #cmds <- sapply(c(165,159,78,167,116,57,166),function(i){
    sprintf("Rscript /home/ob219/git/basis_paper/GWAS/eqtlgen/process_eqtlgen.R -f %s",f)
  })
  write(cmds,file="~/tmp/qstuff/gwas_eqtlgen_notrans.txt")
}

shrink.DT <- readRDS(SHRINKAGE_FILE)
shrink.DT<-shrink.DT[,c('pid',shrink=SHRINKAGE_METHOD),with=FALSE]
setkey(shrink.DT,'pid')
pc.emp <- readRDS(BASIS_FILE)
snp.DT <- fread(SNP_MANIFEST_FILE)
dummy.DT <- snp.DT[,.(pid,uid='DUMMY:-999',shrunk.beta=0)]

all.DT <- lapply(gene.list,function(ensg){
  lf <- file.path(DATA_DIR,sprintf("%s.tab",ensg))
  DT <- fread(lf)[,uid:=ensg]
  merge(DT,shrink.DT,by.x='pid',by.y='pid')[,shrunk.beta:=or * get(`SHRINKAGE_METHOD`)][,.(pid,uid,shrunk.beta)]
}) %>% rbindlist



all.DT <- rbind(all.DT,dummy.DT)
all.DT <- melt(all.DT,id.vars=c('pid','uid'),measure.vars='shrunk.beta')
setkey(all.DT,'pid')
r.DT <- dcast(all.DT,pid~uid+variable,fill=0)
mat <- as.matrix(r.DT[,-1])
rownames(mat) <- r.DT[[1]]
bc <- predict(pc.emp,newdata=t(mat))
res.DT <- data.table(trait = rownames(bc)  %>% gsub("_shrunk.beta","",.),bc)[trait!='DUMMY:-999',]


fname <- basename(args$fname) %>% gsub("\\.txt","",.)

saveRDS(res.DT,file=sprintf("%s%s.RDS",OUT_DIR,fname))
message(sprintf("Wrote to %s%s.RDS",OUT_DIR,fname))
#file.remove(file.path(odir,med$ofile[i]))



if(FALSE){
  library(cowplot)
  OUT_DIR <- '/home/ob219/share/as_basis/GWAS/eqtlgen_projections'
  BASIS_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/basis_june10k.RDS'
  fs <- list.files(path=OUT_DIR,pattern="*.RDS",full.names=TRUE)
  res.DT <- lapply(fs,readRDS) %>% rbindlist
  ## define a Z score for loading to see if any are significant across pc's
  M <- melt(res.DT,id.vars='trait')
  M[,Z:=(value-mean(value))/sqrt(var(value)),by='variable']
  ## use biomart to get gene information of name, entrez_id and biotype
  library(biomart)
  test <- M[variable=='PC3',]

}
