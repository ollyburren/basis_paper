## a list of all phenotypes is available from the paper
## supplementary table 1 - I downloaded this and saved as a data.table

OUT.DIR <- '/home/ob219/share/as_basis/GWAS/geneatlas/0619'
meta.dt <- fread("~/tmp/41588_2018_248_MOESM3_ESM.csv")
meta.dt <- meta.dt[Category=='Binary',.(ID,Description,Cases,Controls=round(Cases/Sample),prop=Sample)]

all.urls <- sprintf("ftp://ftp.igmm.ed.ac.uk/pub/GeneATLAS/%s.v2.tar",meta.dt$ID)

for(i in 1:nrow(meta.dt)){
  trait <- meta.dt$Description[i] %>% make.names
  res.file <- sprintf("%s/%s.RDS",OUT.DIR,trait)
  if(file.exists(res.file))
    next
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

  tmp.file <- file.path(OUT.DIR,'tmp.gz')
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
setnames(all.results,c('SNP','ALLELE','iscores','beta','seb','pval','pid'))

## can we convert my pids to rs#

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
pc.emp <- readRDS(BASIS_FILE)

## create a filter lookup to make things quicker

if(FALSE){
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

#wout <- readRDS("/home/ob219/share/as_basis/GWAS/geneatlas/lookup.RDS")
#filt <-  merge(all.results,wout,by='SNP')
## no need to align ? all wrt to allele2
convertORscale <- function(x,cp) x/(cp * (1-cp))
all.results[,c('beta.log','se.beta.log'):=list(convertORscale(beta,meta.dt$pro[i]),convertORscale(seb,meta.dt$prop[i]))]
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
