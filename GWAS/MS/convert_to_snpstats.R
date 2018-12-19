## convert MS data to snpStats objects and add in phenotypic data

library(annotSnpStats)

if(FALSE){
  ## here we prepare a dossier sample information
  clin <- fread("/home/ob219/share/Data/GWAS/ms-gwas-data/sample_data/DataSetInfo_V9_Clinical.csv")
  md <- fread("/home/ob219/share/Data/GWAS/ms-gwas-data/sample_data/DataSetInfo_V9_11376.csv")

  clin <- clin[,.(Code,Gender,MOB,YOB,AAO,AGE,EDSS,MSSS,PPMS,Severity)]
  md <- md[,.(Code,Region,Ethnicity,md.Gender=Gender,AffectionStatus)]
  M <- merge(md,clin,by='Code',all.x=TRUE)
  ## check that genders match between tables
  M[Gender!=md.Gender,] ## one sample that is Aus so no problem as we use UK only
  saveRDS(M,"/home/ob219/share/Data/GWAS/ms-gwas-data/sample_data/Clinical.RDS")
}

meta <- readRDS("/home/ob219/share/Data/GWAS/ms-gwas-data/sample_data/Clinical.RDS")

## process a dir at a time

# DATA.DIR <- '/home/ob219/share/Data/GWAS/ms-gwas-data/CCC2_Cases_UKN/'
# DATA.DIR <- '/home/ob219/share/Data/GWAS/ms-gwas-data/CCC2_Cases_UKC/'
# DATA.DIR <- '/home/ob219/share/Data/GWAS/ms-gwas-data/CCC2_Cases_UKW/'
# DATA.DIR <- '/home/ob219/share/Data/GWAS/ms-gwas-data/CCC2_Controls_Illu58C/'
DATA.DIR <- '/home/ob219/share/Data/GWAS/ms-gwas-data/CCC2_Controls_Illu58C/'




OUT.DIR <- gsub("CCC2","annotsnpstats-CCC2",DATA.DIR)


gt.files <- list.files(path=DATA.DIR,pattern="*.gz",full.names=TRUE)
sample.file <- list.files(path=DATA.DIR,pattern="*.sample",full.names=TRUE)

library(parallel)
mclapply(gt.files,function(f){
  chr <- basename(f) %>% gsub(".*([0-9][^\\_]+)\\_.*","\\1",.)
  sprintf("Processing chromosome %s",chr) %>% message
  s <- fread(sample.file)[ID_1!=0,]
  s[,sid:=1:.N]
  sample.DT <- merge(s,meta,by.x='ID_1',by.y='Code',all.x=TRUE)[order(sid),]
  X <- read.impute(f)
  ## need snp information from this file
  snps<-sprintf("zcat %s",f) %>% fread(.,select=1:5)
  setnames(snps,c('snpid','rs','position','allele.1','allele.2'))
  snps[,chr:=chr]
  setcolorder(snps,c('snpid','rs','chr','position','allele.1','allele.2'))
  snps.df <- data.frame(snps)
  rownames(snps.df) <- snps$rs
  rownames(X) <- sample.DT$ID_1
  samples.df <- data.frame(sample.DT)
  rownames(samples.df) <- sample.DT$ID_1
  Xf <- new("aSnpMatrix",
      .Data=X,
      snps=snps.df,
      samples=samples.df,
      alleles=c("allele.1","allele.2"),
      phenotype="AffectionStatus")
  fname <- basename(f) %>% gsub(".gen.gz","",.) %>% sprintf("annotsnpstats-%s.RDS",.) %>% file.path(OUT.DIR,.)
  saveRDS(Xf,file=fname)
},mc.cores=8)
