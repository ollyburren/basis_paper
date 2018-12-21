##code to process sun data and work out which variants to take through
##to next stage - not for the feintharted - 933M variants

PQTL_DIR <- '/home/ob219/share/as_basis/sun_pqtl/gwas_basis_june10k_pqtl'
## remove dirs that we have already processed
all.dirs <- list.dirs(path=PQTL_DIR,recursive = FALSE)

af <- list.files(path=PQTL_DIR,pattern="*.tsv.gz",full.names=TRUE,recursive=TRUE)

f <- af[1]
library(parallel)
all.p <- mclapply(af,function(f){
  message(f)
  prot <- dirname(f) %>% basename
  dt <- fread(sprintf("zcat %s",f),select=c(2,3,8),col.names=c('chr','pos','lp'))
  dt[,prot:=prot]
  dt[,.(prot,pid=paste(chr,pos,sep=':'),p=10^lp)]
},mc.cores=8)

files <- data.table(fname=af,begin_code=dirname(af) %>% basename %>% substr(.,1,3))

table(files$begin_code)

byletter <- split(files$fname,files$begin_code)

for(letter in names(byletter)){
ofile <- file.path('/home/ob219/share/as_basis/sun_pqtl/fdr_source',sprintf("%s.RDS",letter))
if(file.exists(ofile)){
  message("Already done")
  next
}
all.p <- mclapply(byletter[[letter]],function(f){
  message(f)
  prot <- dirname(f) %>% basename
  dt <- fread(sprintf("zcat %s",f),select=c(2,3,8),col.names=c('chr','pos','lp'))
  dt[,prot:=prot]
  dt[,.(prot,pid=paste(chr,pos,sep=':'),p=10^lp)]
},mc.cores=8) %>% rbindlist

saveRDS(all.p,file=ofile)
}

library(parallel)
all.files <- list.files(path='/home/ob219/share/as_basis/sun_pqtl/fdr_source',pattern='*.RDS',full.names=TRUE)
all.dat<-lapply(all.files,readRDS)

by.chr<-split(all,all$chr)

filtered <- mclapply(seq_along(by.chr),function(i){
  DT <- by.chr[[i]]
  sprintf("Processing %d",i) %>% message
  DT[p.adjust(p,method='fdr')<0.05,]
},mc.cores=8)

filtered <- rbindlist(filtered)
saveRDS(filtered,file="/home/ob219/share/as_basis/GWAS/sun_pqtl/all_by_chr_fdr0.05.RDS")
