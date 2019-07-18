## a list of all phenotypes is available from the paper
## supplementary table 1 - I downloaded this and saved as a data.table

library(optparse)

TEST<-FALSE
option_list = list(
        make_option(c("-p", "--phenotypes"), type="character", default=NULL,
              help="file containing list of phenotypes to process", metavar="character")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$phenotypes)){
	   print_help(opt_parser)
	    stop("Supply a file containing a list of phenotypes to process", call.=FALSE)
    }
}else{
  args <- list(phenotypes = '/home/ob219/tmp/qstuff/geneatlas/file48cf2588d1f17')
}


if(FALSE){
  ## create a list of phenotype ids
  BLOCK.SIZE<-50
  meta.dt <- fread("~/tmp/41588_2018_248_MOESM3_ESM.csv")
  meta.dt <- meta.dt[Category=='Binary',.(ID,Description,Cases,Controls=round(Cases/Sample)-Cases,prop=Sample)]
  lout <- '/home/ob219/tmp/qstuff/geneatlas/'
  foo<-lapply(split(meta.dt$ID,ceiling(seq_along(meta.dt$ID)/BLOCK.SIZE)),function(ids){
    tfile <- tempfile(tmpdir=lout)
    write(ids,file=tfile)
    tfile
  })
  rscript <- '/home/ob219/git/basis_paper/GWAS/GENEATLAS/process_gene_atlas.R'
  cmds <- sapply(foo,function(x){
    sprintf("Rscript %s -p %s",rscript,x)
  })
  write(cmds,file="~/tmp/qstuff/geneatlas.txt")


}


phids<-scan(args$phenotypes,"character()")



OUT.DIR <- '/home/ob219/share/as_basis/GWAS/geneatlas/0619'
meta.dt <- fread("~/tmp/41588_2018_248_MOESM3_ESM.csv")
meta.dt <- meta.dt[Category=='Binary',.(ID,Description,Cases,Controls=round(Cases/Sample)-Cases,prop=Sample)]

all.urls <- sprintf("ftp://ftp.igmm.ed.ac.uk/pub/GeneATLAS/%s.v2.tar",meta.dt$ID)


## this code make lookup files to speed things up and not use so much memory
if(FALSE){
  trait <- meta.dt$Description[i] %>% make.names
  res.file <- sprintf("%s/%s.RDS",OUT.DIR,trait)
  trait <- meta.dt$Description[i] %>% make.names
  ftp_url <- all.urls[i]
  out_dir <- '/home/ob219/share/Data/GWAS-summary/tmp'
  ofile <- paste(out_dir,basename(ftp_url),sep='/')
  sprintf("Processing %s",ofile) %>% message
  #ftp_url <- 'ftp://ftp.igmm.ed.ac.uk/pub/GeneATLAS/clinical_c_M05.v2.tar'
  cmd <- sprintf("wget -nv %s -O %s",ftp_url,ofile)
  system(cmd)
  ## next list the files
  cmd <- sprintf("tar -tf %s",ofile)
  chr.files <- system(cmd,intern=TRUE)
  chr.files <- chr.files[grep("\\.gz$",chr.files)]
  chr.files <- chr.files[grep("chr[0-9]+",chr.files)]
  chr.files <- chr.files[grep("imputed",chr.files)]
  #tmp.file <- file.path(OUT.DIR,'tmp.gz')
  tmp.file <- tempfile(pattern = basename(ftp_url) %>% gsub("\\.tar","",.), tmpdir = OUT.DIR, fileext = ".gz")
  all.results <- lapply(chr.files,function(x){
    #cmd3 <- sprintf("tar -Oxf %s %s | zcat",ofile,x)
    # load in the lookup for the chromosome
    cmd3 <- sprintf("tar -Oxf %s %s > %s",ofile,x,tmp.file)
    system(cmd3)
    message(gsub(".*(chr[0-9]+).csv.gz","\\1",x))
    chrom <- gsub(".*(chr[0-9]+).csv.gz","\\1",x)
    fname <- sprintf('/home/ob219/share/as_basis/GWAS/geneatlas/%s_lookup.RDS',chrom)
    lu <- readRDS(fname)
    cmd4 <- sprintf("zcat %s",tmp.file)
    dt <- fread(cmd4)[SNP %in% lu$SNP,]
    dt <- merge(dt,lu,by='SNP')
    #dt[,chr:=chrom]
  }) %>% rbindlist
  unlink(tmp.file)
  setnames(all.results,c('SNP','ALLELE','iscores','beta','seb','pval','pid'))
  uk10 <- readRDS("/home/ob219/rds/hpc-work/DATA/UK10K/UK10K_0.005_MAF.RDS")
  m <- merge(all.results,uk10,by.x='SNP',by.y='ID')
  m[,pid:=paste(CHROM,POS,sep=':')]
  ## how many from current basis do we get ?
  m.DT <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab')
  filt <- m[pid %in% m.DT$pid,]
  wout <- filt[,.(SNP,pid,chr)]
  lapply(wout$chr %>% unique,function(x){
    fname <- sprintf('/home/ob219/share/as_basis/GWAS/geneatlas/%s_lookup.RDS',x)
    saveRDS(wout[chr==x,.(SNP,pid)],fname)
  })
}

mainfunc <- function(id){
  i <- which(meta.dt$ID==id)
  trait <- meta.dt$Description[i] %>% make.names
  res.file <- sprintf("%s/%s.RDS",OUT.DIR,trait)
  if(file.exists(res.file)){
    sprintf("Already done %s skipping",trait) %>% message
    return()
    message("here")
  }else{
    message("Go")
  }
  trait <- meta.dt$Description[i] %>% make.names
  ftp_url <- all.urls[i]
  out_dir <- '/home/ob219/share/Data/GWAS-summary/tmp'
  ofile <- tempfile(pattern = "file", tmpdir = out_dir, fileext = "")
  ofile <- paste(out_dir,basename(ftp_url),sep='/')
  sprintf("Processing %s",ofile) %>% message
  #ftp_url <- 'ftp://ftp.igmm.ed.ac.uk/pub/GeneATLAS/clinical_c_M05.v2.tar'
  cmd <- sprintf("wget -nv %s -O %s",ftp_url,ofile)
  system(cmd)
  ## next list the files
  cmd <- sprintf("tar -tf %s",ofile)
  chr.files <- system(cmd,intern=TRUE)
  chr.files <- chr.files[grep("\\.gz$",chr.files)]
  chr.files <- chr.files[grep("chr[0-9]+",chr.files)]
  chr.files <- chr.files[grep("imputed",chr.files)]
  #tmp.file <- file.path(OUT.DIR,'tmp.gz')
  tmp.file <- tempfile(pattern = basename(ftp_url) %>% gsub("\\.tar","",.), tmpdir = OUT.DIR, fileext = ".gz")
  all.results <- lapply(chr.files,function(x){
    #cmd3 <- sprintf("tar -Oxf %s %s | zcat",ofile,x)
    # load in the lookup for the chromosome
    cmd3 <- sprintf("tar -Oxf %s %s > %s",ofile,x,tmp.file)
    system(cmd3)
    message(gsub(".*(chr[0-9]+).csv.gz","\\1",x))
    chrom <- gsub(".*(chr[0-9]+).csv.gz","\\1",x)
    fname <- sprintf('/home/ob219/share/as_basis/GWAS/geneatlas/%s_lookup.RDS',chrom)
    lu <- readRDS(fname)
    cmd4 <- sprintf("zcat %s",tmp.file)
    dt <- fread(cmd4)[SNP %in% lu$SNP,]
    dt <- merge(dt,lu,by='SNP')
    #dt[,chr:=chrom]
  }) %>% rbindlist
  unlink(tmp.file)
  setnames(all.results,c('SNP','ALLELE','iscores','beta','seb','pval','pid'))
  BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
  pc.emp <- readRDS(BASIS_FILE)
  convertORscale <- function(x,cp) x/(cp * (1-cp))
  all.results[,c('beta.log','se.beta.log'):=list(convertORscale(beta,meta.dt$prop[i]),convertORscale(seb,meta.dt$prop[i]))]
  SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  tmp <- merge(all.results,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * tmp$beta.log
  ## where snp is missing make it zero
  tmp[,trait:=meta.dt$Description[i] %>% make.names]
  tmp[is.na(metric),metric:=0]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
    stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  out.DT <- data.table(trait=rownames(all.proj),all.proj)
  saveRDS(out.DT,sprintf("%s/%s.RDS",OUT.DIR,unique(out.DT$trait)))
  sprintf("Saved %s/%s.RDS",OUT.DIR,unique(out.DT$trait)) %>% message()
  unlink(ofile)
  sprintf("Deleted %s",ofile) %>% message
}

#clinical_c_N95.v2

for(id in phids){
  tryCatch(mainfunc(id),error=function(e){print(sprintf("Error=%s fread phen=%d",e,id));return(NA)})
}
