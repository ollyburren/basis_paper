library(annotSnpStats)

## 6 GWAS to project one of each of these for PR3 and MPO
#1 gwas1
#2 gwas2
#3 meta (gwas1 and gwas2)

## create a table of sample sizes for each GWAS
#fname label n1  n0
#bolt_gwas1_mpo_bgen.stats.gz  mpo_gwas1  264 5259
#bolt_gwas1_pr3_bgen.stats.gz  pr3_gwas1 478 5259
#bolt_gwas2_mpo_bgen.stats.gz  mpo_gwas2 609 6717
#bolt_gwas2_pr3_bgen.stats.gz  pr3_gwas2 1142  6717
#meta_mpo_lmm.txt  mpo_meta  873 11976
#meta_pr3_lmm.txt  pr3_meta  1620 11976
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
DATA.DIR <- '/home/ob219/share/Data/GWAS-summary/aav_limy_wong'
samples.DT <- fread("/home/ob219/share/Data/GWAS-summary/aav_limy_wong/aav_sample_size.txt")
samples.DT[,n:=n1+n0]
samples.DT[,prop.case:=n1/n]

## meta analysis values are not converted to the OR scale as far as I can see
man.DT <- fread(SNP_MANIFEST)
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
library(parallel)
all.vasc <- mclapply(seq_along(samples.DT),function(i){
  sprintf("Processing %s",samples.DT$fname[i]) %>% message
  fname <- file.path(DATA.DIR,samples.DT$fname[i])
  stat.DT <- fread(sprintf("zcat %s",fname))
  ## in order to convert lmm to the or scale do the following
  ## taken from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/#x1-5000010 (section 10.2)
  convertORscale <- function(x,cp) x/(cp * (1-cp))
  if(grepl("meta",samples.DT$label[i])){
    stat.DT[,c('beta.log','se.beta.log'):=list(convertORscale(beta,samples.DT$prop.case[i]),convertORscale(se,samples.DT$prop.case[i]))]
    ## note the allele recoding
    setnames(stat.DT,c('rsid','allele_A','allele_B'),c('SNP','ALLELE0','ALLELE1'))
  }else{
    stat.DT[,c('beta.log','se.beta.log'):=list(convertORscale(BETA,samples.DT$prop.case[i]),convertORscale(SE,samples.DT$prop.case[i]))]
  }
  ## note here that the counted allele is ALLELE1 hence we flip things as in the basis the counted allele is a2
  M <- merge(stat.DT[SNP %in% man.DT$pid,.(trait=samples.DT$label[i],pid=SNP,a1=ALLELE0,a2=ALLELE1,or=exp(beta.log))],man.DT,by='pid',all.y=TRUE)
  idx <- which(is.na(M$or))
  sprintf("%d missing",length(idx)) %>% message
  if(length(idx)!=0)
    M[idx,c('trait','a1','a2','or'):=list(samples.DT$label[i],ref_a1,ref_a2,or=1)]
  ## align alleles
  alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
  align.class <- rep('match',nrow(alleles))
  idx<-which(alleles$al.x!=alleles$al.y)
  x.alleles <- alleles[idx,]$al.x
  names(x.alleles)<-alleles[idx,]$pid
  y.alleles <-  alleles[idx,]$al.y
  names(y.alleles)<-names(x.alleles)
  align.class[idx] <- g.class(x.alleles,y.alleles)
  print(table(align.class))
  alleles[,g.class:=align.class]
  idx<-which(alleles$g.class=='impossible')
  if(length(idx) >0){
    M <- M[-idx,]
    alleles <- alleles[-idx,]
  }
  ## things are flipped I think but need to check
  flip <- which(alleles$g.class=='rev')
  if(length(flip)>0)
    M[flip,or:=1/or]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  ## where snp is missing make it zero
  tra <- unique(M$trait)
  tmp[is.na(metric),c('metric','trait'):=list(0,tra)]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
},mc.cores=8) %>% do.call('rbind',.)
vasc.DT <- data.table(trait=rownames(all.vasc),all.vasc)
saveRDS(vasc.DT,file="/home/ob219/share/as_basis/GWAS/aav_projections/summary/aav_2019.RDS")


## vasc.DT <- readRDS("/home/ob219/share/as_basis/GWAS/wong_aav/projections/aav_2019.RDS")
vasc.DT <- melt(vasc.DT,id.vars='trait')
control.DT <- data.table(rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='V1')
control.DT <- control.DT[V1=='control',.(pc=variable,control.score=value)]
vasc.DT <- merge(vasc.DT,control.DT,by.x='variable',by.y='pc')
vasc.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
library(cowplot)
ggplot(vasc.DT,aes(x=variable,y=value-control.score,color=trait,group=trait)) + geom_point() + geom_line() +
geom_hline(yintercept=0,lty=2)

##

ggplot(vasc.DT[grep("meta",trait),],aes(x=variable,y=value-control.score,color=trait,group=trait)) + geom_point() + geom_line() +
geom_hline(yintercept=0,lty=2)
