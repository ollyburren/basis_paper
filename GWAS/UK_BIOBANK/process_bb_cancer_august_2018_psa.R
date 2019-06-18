## work out the list of variants reqiured

library(optparse)

TEST<-FALSE
option_list = list(
        make_option(c("-i", "--integer"), type="numeric", default=NULL,
              help="index of phenotype to process ", metavar="numeric")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$integer)){
	   print_help(opt_parser)
	    stop("Supply an integer for phenotype to process", call.=FALSE)
    }
}else{
  args <- list(integer=10)
}

i<-args$integer

## to use maf estimate of se remove ss prefix !

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
#SHRINKAGE_METHOD<-'recip.emp_maf_se'
#SHRINKAGE_METHOD<-'none'
## just the one shrinkage file
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_vit_t2d.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_vit_t2d.RDS'
#BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_noshrink_gwas.RDS'
#BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS'
BNEALE_DIR <- '/home/ob219/share/Data/GWAS-summary/uk_biobank_neale_summary_stats_2018/'
BASIS_FILT_DIR <- file.path(BNEALE_DIR,'as_basis_tmp')
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
BB_BASIS_LU_FILE <- '/home/ob219/share/as_basis/GWAS/support//sept_bb_gwas_var_man.RDS'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/vit_t2d_2019/'
#OUT_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/beta_2018/'


## running on the queue
if(FALSE){
  bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
  pheno <- fread(bb_phenofile)
  setnames(pheno,names(pheno) %>% make.names)
  med <- pheno[grepl("20001\\_",Phenotype.Code) & Sex=='both_sexes',]
  cmds <- sapply(1:nrow(med),function(i){
  #cmds <- sapply(c(165,159,78,167,116,57,166),function(i){
    sprintf("Rscript /home/ob219/git/basis_paper/GWAS/process_bb_cancer_august_2018.R -i %d",i)
  })
  write(cmds,file="~/tmp/qstuff/gwas_bb_cancer_proj_t2d_vit.txt")
}


## generating lookup file
if(FALSE){
  snp.DT <- fread(SNP_MANIFEST_FILE)
  bb.snps.DT <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/variants.tsv')
  tmp <- bb.snps.DT[,.(pid=paste(chr,pos,sep=':'),varid,rsid,bb_ref=ref,bb_alt=alt)]
  tmp[,lookup:=paste(pid,bb_ref,bb_alt,sep=':')]
  M<-merge(snp.DT,tmp,by.x='pid',by.y='pid',all.x=TRUE)

  ## note that there are 25 SNPs missing that were included before so will need to update basis
  bb_man <- M[!is.na(varid),]
  bb_man[ref_a1==bb_ref & ref_a2==bb_alt,flip:=FALSE]
  bb_man[ref_a2==bb_ref & ref_a1==bb_alt,flip:=TRUE]
  ## everything matches so no flipping required.
  saveRDS(bb_man$lookup,BB_BASIS_LU_FILE)
}
keep <- readRDS(BB_BASIS_LU_FILE) %>% gsub("\\_",':',.)


## load in biobank link file
bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
pheno <- fread(bb_phenofile)
setnames(pheno,names(pheno) %>% make.names)
med <- pheno[grepl("20001\\_",Phenotype.Code) & Sex=='both_sexes',]

## load in phenotype file

P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/phenotypes.both_sexes.tsv')
P<-P[,.(phenotype,variable_type,non_missing=n_non_missing,cases=n_cases,controls=n_controls,pheno.source=source)]


## load in manifest
## compose a new command
med[,c('wget','db','o','ofile'):=tstrsplit(wget.command,' ')]
#odir <- '/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/self_reported_disease/'
med[,new.cmd:=sprintf("wget -nv %s -O %s%s",Dropbox.File,BNEALE_DIR,ofile)]
med[,phe:=make.names(Phenotype.Description) %>% gsub("Cancer.code..self.reported..","",.)]

## download data if required
if(!file.exists(file.path(BNEALE_DIR,med$ofile[i]))){
  message(med$new.cmd[i])
  system(med$new.cmd[i],ignore.stdout=TRUE,ignore.stderr=FALSE)
}
filt.fname <- gsub("\\.tsv\\.bgz",".RDS",med$ofile[i])
if(!file.exists(file.path(BASIS_FILT_DIR,filt.fname))){
  message(sprintf("zcat %s",file.path(BNEALE_DIR,med$ofile[i])))
  DT<-fread(sprintf("zcat %s",file.path(BNEALE_DIR,med$ofile[i])),showProgress=FALSE)[variant %in% keep, ]
  saveRDS(DT,file=file.path(BASIS_FILT_DIR,filt.fname))
  q(save="no")
}else{
  DT <- readRDS(file.path(BASIS_FILT_DIR,filt.fname))
}




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

p <- P[phenotype==med$Phenotype.Code[i],]

DT[,or:=computeOR(p$non_missing,p$cases,AC,ytx)]
DT[,c('theta','se.theta'):=list(log(or),SElor(p$non_missing,p$cases,AC,ytx))]
DT[,c('theta.pval','theta.Z','n0','n1'):=list(2*(pnorm(abs(theta/se.theta),lower.tail = FALSE)),theta/se.theta,p$non_missing-p$cases,p$cases)]


out <- DT[,.(variant,or,p.value=theta.pval)]
out[,c('chr','pos','a1','a2'):=tstrsplit(variant,':')]
out[,pid:=paste(chr,pos,sep=':')]

snp.DT <- fread(SNP_MANIFEST_FILE)
out<-merge(snp.DT,out,by.x='pid',by.y='pid')

## if any of the OR are 0 or inf set these to 1
out[or==0 | is.infinite(or),or:=1]




shrink.DT <- readRDS(SHRINKAGE_FILE)
shrink.DT<-shrink.DT[,c('pid',shrink=SHRINKAGE_METHOD),with=FALSE]
setkey(shrink.DT,'pid')
pc.emp <- readRDS(BASIS_FILE)

all.DT <- out[!duplicated(pid),.(pid,uid=med$phe[i],beta=log(or))]


## add shrinkage here
all.DT <- merge(all.DT,shrink.DT,by.x='pid',by.y='pid')[,shrunk.beta:=beta * get(`SHRINKAGE_METHOD`)][,.(pid,uid,shrunk.beta)]


## next add in a dummy gene that has data for all variants
dummy.DT <- snp.DT[,.(pid,uid='DUMMY:-999',shrunk.beta=0)]
all.DT <- rbind(all.DT,dummy.DT)
all.DT <- melt(all.DT,id.vars=c('pid','uid'),measure.vars='shrunk.beta')
setkey(all.DT,'pid')
r.DT <- dcast(all.DT,pid~uid+variable,fill=0)
mat <- as.matrix(r.DT[,-1])
rownames(mat) <- r.DT[[1]]
bc <- predict(pc.emp,newdata=t(mat))
res.DT <- data.table(trait = med$phe[i]  %>% gsub("_shrunk.beta","",.),bc)[2,]


saveRDS(res.DT,file=sprintf("%s%s.RDS",OUT_DIR,med$phe[i]))
message(sprintf("Wrote to %s%s.RDS",OUT_DIR,med$phe[i]))
#file.remove(file.path(odir,med$ofile[i]))
