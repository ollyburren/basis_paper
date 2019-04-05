## work out the list of variants reqiured

library(optparse)

TEST<-TRUE
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
  args <- list(integer=136)
}

i<-args$integer

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
#SHRINKAGE_METHOD<-'recip.emp_maf_se'
## just the one shrinkage file
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/basis_gwas.RDS'
#BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/basis_noshrink_gwas.RDS'
BNEALE_DIR <- '/home/ob219/share/Data/GWAS-summary/uk_biobank_neale_summary_stats_2017/'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
BB_BASIS_LU_FILE <- '/home/ob219/share/as_basis/GWAS/support//2017_bb_gwas_var_man.RDS'
#OUT_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/no_shrink_2017/'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/shrink_2017/'


## running on the queue
if(FALSE){
  bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_linklist.20170915.csv'
  pheno <- fread(bb_phenofile)
  setnames(pheno,names(pheno) %>% make.names)
  med <- pheno[grepl("20002\\_",Phenotype.code),]
  cmds <- sapply(1:nrow(med),function(i){
  #cmds <- sapply( c(187,181,86,189,136,64,188),function(i){
      sprintf("Rscript /home/ob219/git/basis_paper/GWAS/process_bb_self_reported_disease_september_2017.R -i %d",i)
  })
  write(cmds,file="~/tmp/qstuff/gwas_bb_disease_proj_2017.txt")
}


## generating lookup file
if(FALSE){
  snp.DT <- fread(SNP_MANIFEST_FILE)
  bb.snps.DT <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20170915/variants.tsv')
  bb.snps.DT[,c('chr','pos','a1','a2'):=tstrsplit(variant,':')]
  tmp <- bb.snps.DT[,.(pid=paste(chr,pos,sep=':'),variant,rsid,bb_ref=a1,bb_alt=a2)]
  tmp[,lookup:=paste(pid,bb_ref,bb_alt,sep=':')]
  M<-merge(snp.DT,tmp,by.x='pid',by.y='pid',all.x=TRUE)

  ## note that there are 25 SNPs missing that were included before so will need to update basis
  bb_man <- M[!is.na(variant),]
  bb_man[ref_a1==bb_ref & ref_a2==bb_alt,flip:=FALSE]
  bb_man[ref_a2==bb_ref & ref_a1==bb_alt,flip:=TRUE]
  ## everything matches so no flipping required.
  saveRDS(bb_man$lookup,BB_BASIS_LU_FILE)
}
keep <- readRDS(BB_BASIS_LU_FILE) %>% gsub("\\_",':',.)


## load in biobank link file
bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_linklist.20170915.csv'
pheno <- fread(bb_phenofile)
setnames(pheno,names(pheno) %>% make.names)
med <- pheno[grepl("20002\\_",Phenotype.code),]

## load in phenotype file

P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20170915/phenosummary_final_11898_18597.tsv')
P<-P[,.(Field.code,non_missing=N.non.missing,cases=N.cases,controls=N.controls)]


## load in manifest
## compose a new command
med[,c('wget','db','o','ofile'):=tstrsplit(wget.command,' ')]
#odir <- '/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/self_reported_disease/'
med[,new.cmd:=sprintf("wget %s -O %s%s",db,BNEALE_DIR,ofile)]
med[,phe:=make.names(Description) %>% gsub("Non.cancer.illness.code..self.reported..","",.)]

## download data if required
if(!file.exists(file.path(BNEALE_DIR,med$ofile[i])))
  system(med$new.cmd[i])
DT<-fread(sprintf("zcat %s",file.path(BNEALE_DIR,med$ofile[i])))[variant %in% keep, ]




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

p <- P[Field.code==med$Phenotype.code[i],]

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


if(FALSE){
  old.out <- fread('/home/ob219/share/as_basis/GWAS/sum_stats/bb:20002_1222:self_reported_type_1_diabetes.tab')
  setnames(old.out,paste('old',names(old.out),sep='.'))
  M <- merge(old.out,out,by.x='old.pid',by.y='pid')
  M[,old.or:=as.numeric(old.or)]


  ## compare with the ones that give great results
  old.out <- fread('/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/bb:20002_1222:self_reported_type_1_diabetes.tab')
  setnames(old.out,paste('old',names(old.out),sep='.'))
  M <- merge(old.out,out,by.x='old.pid',by.y='pid')

}


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

if(FALSE){
  ## some had numerical errors - code to identify and rerun
  OUT_DIR <- '/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/self_reported_disease/june_10k/'
  fs <- list.files(path=OUT_DIR,pattern="*.RDS",full.names=TRUE)
  all.res <- lapply(fs,readRDS) %>% rbindlist
  all.res[is.nan(PC1),]
  missing<-all.res[is.nan(PC1),]$trait
  cmds <- sapply(which(med$phe %in% missing),function(i){
    sprintf("Rscript /home/ob219/git/as_basis/R/Individual_projection/process_bb_self_reported_disease.R -i %d",i)
  })
  write(cmds,file="~/tmp/qstuff/bb_med_proj.txt")
  ## what about infinite ones looks as if there was an issue - that I can't repeat !
  missing<-all.res[is.infinite(PC1),]$trait
  cmds <- sapply(which(med$phe %in% missing),function(i){
    sprintf("Rscript /home/ob219/git/as_basis/R/Individual_projection/process_bb_self_reported_disease.R -i %d",i)
  })
  write(cmds,file="~/tmp/qstuff/bb_med_proj.txt")
}


if(FALSE){
  library(cowplot)
  OUT_DIR <- '/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/self_reported_disease/june_10k/'
  BASIS_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/basis_june10k.RDS'
  fs <- list.files(path=OUT_DIR,pattern="*.RDS",full.names=TRUE)
  res.DT <- lapply(fs,readRDS) %>% rbindlist
  pc.emp <- readRDS(BASIS_FILE)
  basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x,cat='basis')
  res.DT[,cat:='bb']
  plot.DT <- rbind(res.DT,basis.DT)
  plot.DT[,label:=gsub("\\.[0-9]*mg.*","",trait)]
  ggplot(plot.DT,aes(x=PC1,y=PC4,label=label,col=cat)) + geom_point() + geom_text()
  ## compute Z scores for chris
  mdt<-melt(res.DT,id.vars='trait',measure.vars=sprintf("PC%d",1:11))
  mdt[,c('mean','sd'):=list(mean(value),sd(value)),by='variable']
  mdt[,Z:=(value-mean)/sd]
  mdt[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
  mdt[,p.adj:=p.adjust(p.value,method="fdr"),by=variable]
  save(mdt,file="/home/ob219/rds/rds-cew54-wallace-share/as_basis/bb/basis_june10k_mself_reported_disease.RDS")
}
