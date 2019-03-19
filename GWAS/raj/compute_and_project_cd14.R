library(annotSnpStats)
library(parallel)
library(cupcake)
library(optparse)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
CAL_FILE <- '/home/ob219/rds/hpc-work/as_basis/support//por_2500_2.0_0.01.RDS'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'

TEST<-FALSE
option_list = list(
        make_option(c("-i", "--integer"), type="numeric", default=NULL,
              help="index of probesets to process", metavar="numeric")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$integer)){
	   print_help(opt_parser)
	    stop("Supply an integer for phenotype to process", call.=FALSE)
    }
}else{
  args <- list(integer=39)
}

if(FALSE){
  SCRIPT <- '/home/ob219/git/basis_paper/GWAS/raj/compute_and_project_cd14.R'
  cmds <- sapply(1:387,function(i){
    sprintf("Rscript %s -i %s",SCRIPT,i)
  })
  write(cmds,file="~/tmp/qstuff/raj_summ.txt")
}


fs <- list.files(path="/home/ob219/share/as_basis/GWAS/individual_data/filtered_gt/raj/cd14",pattern="^chr[0-9]+.RDS",full.names=TRUE)
all.gt<-lapply(fs,readRDS)

gt <- lapply(all.gt,function(x) x@.Data) %>% do.call('cbind',.) %>% sm
snps <- lapply(all.gt,snps) %>% do.call('rbind',.)
samples <- samples(all.gt[[1]])

merged.gt <- new("aSnpMatrix",
   .Data=gt,
   snps=snps,
   samples=samples,
   alleles=c("allele.1","allele.2"),
   phenotype="affected")

(load("/home/ob219/share/Projects/twas/raj-cd14-expression.RData"))
texpr <- t(expr)
translate <- readRDS('/home/ob219/share/Projects/twas/model_output/trans_raj-cd14.rds')
idxy<-match(rownames(texpr),names(translate))
rownames(texpr) <- translate[idxy]

samp <- samples(merged.gt)
idxz <- match(rownames(samp),rownames(texpr))
colnames(texpr) <- colnames(texpr) %>% paste("P",.,sep='_')
samp <- cbind(samp,texpr[idxz,])

all.probes <- colnames(texpr)

## do in batches of 500

split.probe <- split(all.probes, ceiling(seq_along(all.probes)/50))
all.probes <- split.probe[[args$integer]]

all.lm <- lapply(all.probes,function(prob){
  message(prob)
  res <- snp.rhs.estimates(sprintf("%s~sex",prob) %>% formula,family="gaussian",data=samp,snp.data=sm(merged.gt))
  all.beta <- sapply(res,'[[','beta')
  all.vbeta <- sapply(res,'[[','Var.beta')
  res.DT <- data.table(beta=all.beta,vbeta=all.vbeta,probe=prob,variant=names(all.beta),Z=all.beta/sqrt(all.vbeta))
  res.DT[,variant:=gsub("([^\\.]+)\\..*","\\1",variant)]
  res.DT[,p:=pnorm(abs(Z),lower.tail=FALSE) * 2]
  res.DT
})

pc.emp <- readRDS(BASIS_FILE)
shrink.DT <- readRDS(SHRINKAGE_FILE)

## which genotypes are missing

missing <- shrink.DT[!pid %in% snps(merged.gt)$snp.name,]
missing.DT <- data.table(pid=missing$pid,trait='DUMMY',or=1)

## want data structure with pid,trait(probe),value,or
sample.proj.DT <- lapply(all.lm,function(x){
  probe <- x$probe[1]
  missing.DT[,trait:=probe]
  DT <- rbind(missing.DT,x[,.(pid=snps(merged.gt)$snp.name,trait=probe,or=exp(beta))])
}) %>% rbindlist
setkey(sample.proj.DT,pid)
mat.emp <- create_ds_matrix(sample.proj.DT,shrink.DT,SHRINKAGE_METHOD)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)

all.proj.DT <- data.table(probe=rownames(all.proj),all.proj)
all.proj.m <- melt(all.proj.DT,id.var='probe')
ofile <- file.path("/home/ob219/share/as_basis/GWAS/raj/cd14/summary_projections",sprintf("cd4%s.RDS",args$integer))
saveRDS(all.proj.m,file=ofile)
message(sprintf("Wrote %s",ofile))

if(FALSE){
  dat.dir <- '/home/ob219/share/as_basis/GWAS/raj/cd14/summary_projections'
  dat <- lapply(list.files(path=dat.dir,pattern="*.RDS",full.names=TRUE),readRDS) %>% rbindlist
  dat[,Z:=(value-mean(value))/sd(value),by='variable']

  ## load in individual datasets
  all.lms <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd4/regression_models.RDS")
  names(all.lms)<-paste('PC',1:11,sep="")
  for.reg <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd4/cd4.RDS")
  pc.t <- lapply(seq_along(names(all.lms)),function(i){
    pc<-names(all.lms)[i]
    all.p <- sapply(all.lms[[i]],function(x){
      summary(x)$coefficient["value",3]
    })
    DT <- data.table(PC=pc,probe=names(for.reg)[grep('^P\\_',names(for.reg))],t.stat=all.p)
  })

  pc.t <- rbindlist(pc.t)

  M <- merge(pc.t,dat,by.x=c('PC','probe'),by.y=c('variable','probe'))
  pap <- fread("~/tmp/tableS4_eu_cd4T_cis_fdr05.tsv")
  M[,supp.table:=probe %in% paste('P',pap$PROBESET_ID,sep='_')]
  library(cowplot)
  M[,PC:=paste('PC',1:11)]
  ggplot(M,aes(x=Z,y=t.stat,col=supp.table)) + geom_point(alpha=0.5) +
  facet_wrap(~PC) + xlab("Individual T.stat") + ylab("Summary Empirical Z score")


  ## obvious issue could be that we have the alignment wrong

  cosum <- col.summary(merged.gt)
  cosum.DT <- data.table(SNP=rownames(cosum),cosum)
  ## add the pid
  cosum.DT[,pid:=snps(merged.gt)$snp.name]
  ## load in the snp manifest
  snp.man <- fread(SNP_MANIFEST_FILE)
  Mmaf <- merge(snp.man[,.(pid,ref_a1.af)],cosum.DT[,.(pid,RAF)],by='pid')
}
