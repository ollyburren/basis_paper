library(annotSnpStats)
library(parallel)
library(cupcake)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
CAL_FILE <- '/home/ob219/rds/hpc-work/as_basis/support//por_2500_2.0_0.01.RDS'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'

library(optparse)

TEST<-TRUE
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


fs <- list.files(path="/home/ob219/share/as_basis/GWAS/individual_data/filtered_gt/raj/cd4",pattern="^chr[0-9]+.RDS",full.names=TRUE)
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

(load("/home/ob219/share/Projects/twas/raj-cd4-expression.RData"))
texpr <- t(expr)
translate <- readRDS('/home/ob219/share/Projects/twas/model_output/trans_raj-cd4.rds')
idxy<-match(rownames(texpr),names(translate))
rownames(texpr) <- translate[idxy]

samp <- samples(merged.gt)
idxz <- match(rownames(samp),rownames(texpr))
colnames(texpr) <- colnames(texpr) %>% paste("P",.,sep='_')
samp <- cbind(samp,texpr[idxz,])

all.probes <- colnames(texpr)

## do in batches of 500

split.probe <- split(all.probes, ceiling(seq_along(all.probes)/500))
all.probes <- split.probe[[args$integer]]

all.lm <- mclapply(all.probes,function(prob){
  message(prob)
  res <- snp.rhs.estimates(sprintf("%s~sex",prob) %>% formula,family="gaussian",data=samp,snp.data=sm(merged.gt))
  all.beta <- sapply(res,'[[','beta')
  all.vbeta <- sapply(res,'[[','Var.beta')
  res.DT <- data.table(beta=all.beta,vbeta=all.vbeta,probe=prob,variant=names(all.beta),Z=all.beta/sqrt(all.vbeta))
  res.DT[,variant:=gsub("([^\\.]+)\\..*","\\1",variant)]
  res.DT[,p:=pnorm(abs(Z),lower.tail=FALSE) * 2]
  res.DT
},mc.cores=8)

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
saveRDS(all.proj.m,file=file.path("/home/ob219/share/as_basis/GWAS/raj/cd4/summary_projections",paste("cd4%s.RDS",args$integer))
