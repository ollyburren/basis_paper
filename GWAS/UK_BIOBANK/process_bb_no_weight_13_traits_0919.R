## work out the list of variants reqiured


library(parallel)

## to use maf estimate of se remove ss prefix !


## just the one shrinkage file
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
#BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_no_weight_gwas_13_traits_0919.RDS'
BNEALE_DIR <- '/home/ob219/share/Data/GWAS-summary/uk_biobank_neale_summary_stats_2018/'
BASIS_FILT_DIR <- file.path(BNEALE_DIR,'as_basis_tmp')
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'
BB_BASIS_LU_FILE <- '/home/ob219/share/as_basis/GWAS/support/gwas_13_traits_0919_var_man.RDS'
#SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr_13_traits_0919/'SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr_no_weight_13_traits_0919/'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/no_weight_13_traits_0919/'

## load in biobank link file
bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
pheno <- fread(bb_phenofile)
setnames(pheno,names(pheno) %>% make.names)
med <- pheno[grepl("20001\\_|20002\\_|20003\\_",Phenotype.Code) & Sex=='both_sexes',]

## load in phenotype file

P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/phenotypes.both_sexes.tsv')
P<-P[,.(phenotype,variable_type,non_missing=n_non_missing,cases=n_cases,controls=n_controls,pheno.source=source)]


## load in manifest
## compose a new command
med[,c('wget','db','o','ofile'):=tstrsplit(wget.command,' ')]
#odir <- '/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/self_reported_disease/'
med[,new.cmd:=sprintf("wget -nv %s -O %s%s",Dropbox.File,BNEALE_DIR,ofile)]
med[,phe:=make.names(Phenotype.Description) %>% gsub("Cancer.code..self.reported..","SRC:",.)]
med[,phe:=gsub("Treatment.medication.code..","SRM:",phe)]
med[,phe:=gsub("Non.cancer.illness.code..self.reported..","SRD:",phe)]

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


## this assumes that previous process_bb.R has been run so that files
## are downloaded and filtered
pc.emp <- readRDS(BASIS_FILE)
man.DT <- fread(SNP_MANIFEST_FILE)
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,none)]
ids <- grepl("vitiligo|crohns disease|ulcerative colitis|malabsorption/coeliac disease|multiple sclerosis|rheumatoid arthritis|systemic lupus erythematosis|type 1 diabetes|asthma",med$Phenotype.Description) %>% which

all.results <- mclapply(ids,function(i){
  sprintf("Processing %s",med$phe[i]) %>% message
  filt.fname <- gsub("\\.tsv\\.bgz",".RDS",med$ofile[i])
  DT <- readRDS(file.path(BASIS_FILT_DIR,filt.fname))
  p <- P[phenotype==med$Phenotype.Code[i],]
  DT[,or:=computeOR(p$non_missing,p$cases,AC,ytx)]
  DT[,c('theta','se.theta'):=list(log(or),SElor(p$non_missing,p$cases,AC,ytx))]
  DT[,c('theta.pval','theta.Z','n0','n1'):=list(2*(pnorm(abs(theta/se.theta),lower.tail = FALSE)),theta/se.theta,p$non_missing-p$cases,p$cases)]
  out <- DT[,.(variant,or,p.value=theta.pval)]
  out[,c('chr','pos','a1','a2'):=tstrsplit(variant,':')]
  out[,pid:=paste(chr,pos,sep=':')]
  dat <- out[pid %in% man.DT$pid,.(pid,or,p.value)]
  tmp <- merge(dat,stmp,by='pid',all.y=TRUE)
  tmp <- tmp[!duplicated(pid),]
  tmp$metric <- tmp[['none']] * log(tmp$or)
  tmp$trait <- med$phe[i]
  pfile <- file.path(SRC_OUT_DIR,sprintf("UKBB_NEALE:%s_source.RDS",med$phe[i]))
  saveRDS(tmp[,.(pid,or,p.value,none)],file=pfile)
  tmp[is.na(metric) | ! is.finite(metric),c('metric','trait'):=list(0,trait)]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
},mc.cores=8) %>% do.call('rbind',.)


## quick look at hclust

rbind(all.results,pc.emp$x) %>% dist %>% hclust %>% plot
saveRDS(rbind(all.results,pc.emp$x),file='/home/ob219/share/as_basis/GWAS/RESULTS/13_traits_no_weights_0919.RDS')

saveRDS(all.results,'/home/ob219/share/as_basis/GWAS/RESULTS/uk_bb_for_fdr_13_traits_0919.RDS')
