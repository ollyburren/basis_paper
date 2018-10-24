library(annotSnpStats)
DAT.DIR <- '/home/ob219/share/as_basis/ichip/individual_gt/jia/'

files <- list.files(path=DAT.DIR,pattern="*.RDS",full.names=TRUE)
## read in all the genotype

all.dat <- lapply(files,readRDS)
names(all.dat) <- basename(files)



dat.DT <- readRDS("/home/ob219/share/as_basis/ichip/ind_sample_info/jia_samples_info.RDS")

f <- files[1]
tmp <- readRDS(f)
sample.DT <- samples(tmp) %>% data.table
sample.DT[,uid:=1:.N]
sample.DT <- merge(sample.DT,dat.DT[,.(sample_id,ilar_final)],by.x="member",by.y="sample_id")[order(uid),]

out.dir <- '/home/ob219/share/as_basis/ichip/individual_data/filtered_gt/'
by.type <- split(sample.DT$uid,sample.DT$ilar_final)

man.DT <- fread('/home/ob219/share/as_basis/ichip/snp_manifest/ichip_all.tab')
for(n in names(by.type)){
  out.path <- file.path(out.dir,n)
  dir.create(out.path, showWarnings = FALSE)
    idx <- by.type[[n]]
  for(fn in names(all.dat)){
    X <- all.dat[[fn]][idx,]
    snp.idx <- which(snps(X)$pid %in% man.DT$pid)
    X <- X[,snp.idx]
    saveRDS(X,file=file.path(out.path,fn))
  }
}

if(FALSE){
## just check the alleles are aligned
ind.cs <- lapply(all.dat, function(x){
  cs <- col.summary(x)
  data.table(pid=rownames(cs),cs)
} ) %>% rbindlist

M<-merge(man.DT,ind.cs[,.(pid,RAF)],by.x='pid',by.y='pid')
plot(M$ref_a1.af,M$RAF)

## all looks good.
}
