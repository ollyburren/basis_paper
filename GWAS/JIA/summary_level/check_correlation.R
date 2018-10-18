##helper code to extract 'interesting' self reported disease data from biobank and get beta's

computeOR <- function(n,n1,Sx,Sxy) {
    ## estimated allele freqs in cases and controls
    fe1 <- Sxy/(2*n1)
    fe0 <- (Sx - Sxy)/(2*(n-n1))
    ## estimated odds ratio
    fe1 * (1-fe0) / ( (1-fe1) * fe0 )
}

SElor<-function(n,n1,Sx,Sxy){
    n0<-n-n1
    #fe1 is the af in cases
    c <- Sxy/(2*n1)
    #fe0 is af in the controls
    a <- (Sx - Sxy)/(2*(n-n1))
    b<-1-a
    d<-1-c
    ## normalise
    a<-(a*n0)/n
    b<-(b*n0)/n
    c<-(c*n1)/n
    d<-(d*n1)/n
    ## estimated odds ratio bc/ad
    sqrt(1/2) * sqrt(1/n) * sqrt(1/a + 1/b + 1/c + 1/d)
}


bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
pheno <- fread(bb_phenofile)
setnames(pheno,names(pheno) %>% make.names)
pheno <- pheno[grepl("20002\\_",Phenotype.Code) & Sex=='both_sexes',]
pheno <- pheno[,.(Phenotype.Code,Phenotype.Description,File)]


## load in phenotype file

P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/phenotypes.both_sexes.tsv')
P<-P[,.(phenotype,variable_type,non_missing=n_non_missing,cases=n_cases,controls=n_controls,pheno.source=source)]
P <- merge(pheno,P,by.x='Phenotype.Code',by.y='phenotype')
P[,phe:=make.names(Phenotype.Description) %>% gsub("Non.cancer.illness.code..self.reported..","",.)]

P[,File:=gsub("tsv.bgz$","RDS",File)]


## load in manifest
## compose a new command
med[,phe:=make.names(Phenotype.Description) %>% gsub("Non.cancer.illness.code..self.reported..","",.)]





files <- list.files(path='/home/ob219/share/Data/GWAS-summary/uk_biobank_neale_summary_stats_2018/as_basis_tmp',pattern='*.RDS',full.names=TRUE)
files <- files[grepl('^20002',basename(files))]

for(f in files){
message(f)
DT <-  readRDS(f)
p <- P[File==basename(f),]
DT[,or:=computeOR(p$non_missing,p$cases,AC,ytx)]
DT[,c('theta','se.theta'):=list(log(or),SElor(p$non_missing,p$cases,AC,ytx))]
DT[,c('theta.pval','theta.Z','n0','n1'):=list(2*(pnorm(abs(theta/se.theta),lower.tail = FALSE)),theta/se.theta,p$non_missing-p$cases,p$cases)]
out <- DT[,.(variant,or,p.value=theta.pval)]
out[,c('chr','pos','a1','a2'):=tstrsplit(variant,':')]
out<-out[,pid:=paste(chr,pos,sep=':')][,.(pid,a1,a2,or,p.value)]
ofile <- file.path("/home/ob219/share/as_basis/GWAS/tmp/bb_sum_stats/",paste(p$phe,'tab',sep='.'))
write.table(out,file=ofile,row.names=FALSE,quote=FALSE,sep="\t")
}

## get the list of traits that we want so we can create the matrix for correlation

dat.DT <- readRDS("~/share/as_basis/GWAS/tmp/jia_bb_summary.RDS")
traits <- dat.DT[!grepl('^jia|unclassifiable',trait),]$trait %>% unique
traits <- paste(traits,'tab',sep='.')
tra <- traits[1]
bb_traits <- lapply(traits,function(t){
  DT <- file.path("/home/ob219/share/as_basis/GWAS/tmp/bb_sum_stats/",t) %>% fread
  DT[,trait:=sub(".tab","",t)]
}) %>% rbindlist

## add in those for JIA

jia_tl <- c('jiaeo_unpub.tab','jiaera_unpub.tab','jiapo_unpub.tab','jiapsa_unpub.tab','jiarfneg_unpub.tab','jiarfpos_unpub.tab','jiasys_unpub.tab')

jia_traits <- lapply(jia_tl,function(t){
  DT <- file.path("/home/ob219/share/as_basis/GWAS/sum_stats/",t) %>% fread
  DT[,trait:=sub(".tab","",t)]
}) %>% rbindlist

all <- rbind(jia_traits,bb_traits)

all[,beta:=log(or)]
all<-all[!duplicated(paste(pid,trait)),]

foo <- melt(all,id.vars=c('pid','trait'),measure.vars='beta')
M <- dcast(foo,pid~variable+trait)
Mat <- M[,-1] %>% as.matrix
## remove odd cases where things are missing or have an odd OR
Mat <- Mat[-which(is.na(Mat),arr.ind=TRUE)[,1],]
cfoo<-cor(Mat)
rownames(cfoo) <- gsub("beta\\_|\\_unpub","",rownames(cfoo))
colnames(cfoo) <- gsub("beta\\_|\\_unpub","",colnames(cfoo))
library(pheatmap)
pheatmap(cfoo)
