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
  args <- list(integer=2)
}

i<-args$integer


if(FALSE){
  bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
  pheno <- fread(bb_phenofile)
  setnames(pheno,names(pheno) %>% make.names)
  med <- pheno[grepl("20001\\_",Phenotype.Code) & Sex=='both_sexes',]
  cmds <- sapply(1:nrow(med),function(i){
    sprintf("Rscript /home/ob219/git/as_basis/R/Individual_projection/process_bb_cancer.R -i %d",i)
  })
  write(cmds,file="~/tmp/qstuff/bb_cancer_proj.txt")
}



if(FALSE){
  snp.DT <- fread('/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab')
  bb.snps.DT <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/variants.tsv')
  tmp <- bb.snps.DT[,.(pid=paste(chr,pos,sep=':'),varid,rsid,bb_ref=ref,bb_alt=alt)]
  tmp[,lookup:=paste(pid,bb_ref,bb_alt,sep=':')]
  M<-merge(snp.DT,tmp,by.x='pid',by.y='pid',all.x=TRUE)

  ## note that there are 25 SNPs missing that were included before so will need to update basis
  bb_man <- M[!is.na(varid),]
  bb_man[ref_a1==bb_ref & ref_a2==bb_alt,flip:=FALSE]
  bb_man[ref_a2==bb_ref & ref_a1==bb_alt,flip:=TRUE]
  ## everything matches so no flipping required.
  saveRDS(bb_man$lookup,"/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/june_10k_bb_varid.RDS")
}
keep <- readRDS("/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/june_10k_bb_varid.RDS") %>% gsub("\\_",':',.)


## load in biobank link file
bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
pheno <- fread(bb_phenofile)
setnames(pheno,names(pheno) %>% make.names)
med <- pheno[grepl("20001\\_",Phenotype.Code) & Sex=='both_sexes',]

## load in phenotype file

P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/phenotypes.both_sexes.tsv')
P<-P[,.(phenotype,variable_type,non_missing=n_non_missing,cases=n_cases,controls=n_controls)]


## load in manifest
## compose a new command
med[,c('wget','db','o','ofile'):=tstrsplit(wget.command,' ')]
odir <- '/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/cancer/'
med[,new.cmd:=sprintf("wget %s -O %s%s",Dropbox.File,odir,ofile)]
med[,phe:=make.names(Phenotype.Description) %>% gsub("Cancer.code..self.reported..","",.)]

## reactivate this if we ever want to redownload everthing
if(TRUE)
  system(med$new.cmd[i])


DT<-fread(sprintf("zcat %s",file.path(odir,med$ofile[i])))[variant %in% keep, ]




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

snp.DT <- fread('/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab')
out<-merge(snp.DT,out,by.x='pid',by.y='pid')

## if any of the OR are 0 or inf set these to 1
out[or==0 | is.infinite(or),or:=1]


SHRINKAGE_METHOD<-'ws_emp'
SHRINKAGE_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/shrinkage_june10k.RDS'
BASIS_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/basis_june10k.RDS'

shrink.DT <- readRDS(SHRINKAGE_FILE)
shrink.DT<-shrink.DT[,c('pid',shrink=sprintf("%s_shrinkage",SHRINKAGE_METHOD)),with=FALSE]
setkey(shrink.DT,'pid')
pc.emp <- readRDS(BASIS_FILE)

all.DT <- out[!duplicated(pid),.(pid,uid=med$phe[i],beta=log(or))]


## add shrinkage here
all.DT <- merge(all.DT,shrink.DT,by.x='pid',by.y='pid')[,shrunk.beta:=beta * ws_emp_shrinkage][,.(pid,uid,shrunk.beta)]


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

OUT.DIR <- '/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/cancer/june_10k/'
saveRDS(res.DT,file=sprintf("%s%s.RDS",OUT.DIR,med$phe[i]))
message(sprintf("Wrote to %s%s.RDS",OUT.DIR,med$phe[i]))
file.remove(file.path(odir,med$ofile[i]))

if(FALSE){
  ## some had numerical errors - code to identify and rerun
  OUT.DIR <- '/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/cancer/june_10k/'
  fs <- list.files(path=OUT.DIR,pattern="*.RDS",full.names=TRUE)
  all.res <- lapply(fs,readRDS) %>% rbindlist
  all.res[is.nan(PC1),]
  missing<-all.res[is.nan(PC1),]$trait
  cmds <- sapply(which(med$phe %in% missing),function(i){
    sprintf("Rscript /home/ob219/git/as_basis/R/Individual_projection/process_bb_cancer.R -i %d",i)
  })
  write(cmds,file="~/tmp/qstuff/bb_med_proj.txt")
  ## what about infinite ones looks as if there was an issue - that I can't repeat !
  missing<-all.res[is.infinite(PC1),]$trait
  cmds <- sapply(which(med$phe %in% missing),function(i){
    sprintf("Rscript /home/ob219/git/as_basis/R/Individual_projection/process_bb_cancer.R -i %d",i)
  })
  write(cmds,file="~/tmp/qstuff/bb_med_proj.txt")
}


if(FALSE){
  library(cowplot)
  OUT.DIR <- '/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/cancer/june_10k/'
  BASIS_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/basis_june10k.RDS'
  fs <- list.files(path=OUT.DIR,pattern="*.RDS",full.names=TRUE)
  res.DT <- lapply(fs,readRDS) %>% rbindlist
  pc.emp <- readRDS(BASIS_FILE)
  basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x,cat='basis')
  res.DT[,cat:='bb']
  plot.DT <- rbind(res.DT,basis.DT)
  plot.DT[,label:=gsub("\\.[0-9]*mg.*","",trait)]
  ggplot(plot.DT,aes(x=PC1,y=PC2,label=label,col=cat)) + geom_point() + geom_text()
  ## compute Z scores for chris
  mdt<-melt(res.DT,id.vars='trait',measure.vars=sprintf("PC%d",1:11))
  mdt[,c('mean','sd'):=list(mean(value),sd(value)),by='variable']
  mdt[,Z:=(value-mean)/sd]
  mdt[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
  mdt[,p.adj:=p.adjust(p.value,method="fdr"),by=variable]
  mdt[p.adj<0.05,.(trait,pc=variable,Z,p.value,p.adj)]
  saveRDS(mdt,file="/home/ob219/rds/rds-cew54-wallace-share/as_basis/bb/basis_june10k_cancer.RDS")
}
