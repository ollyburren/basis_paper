library(annotSnpStats)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- "/home/ob219/share/as_basis/GWAS/renton_mg/projections/renton_mg_0619.RDS"
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'

## process main one first
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
pc.emp <- readRDS(BASIS_FILE)

files <- c('renton_mg' = 'MyastheniaGravis_Renton_JAMA_Neurol_2015.txt.gz',
'renton_mg_late' = 'MyastheniaGravis_LateOnset_Renton_JAMA_Neurol_2015.txt.gz',
'renton_mg_early' = 'MyastheniaGravis_YoungOnset_Renton_JAMA_Neurol_2015.txt.gz')
library(parallel)
all.mg <- mclapply(names(files),function(f){
  message(f)
  fname <- file.path('/home/ob219/share/Data/GWAS-summary',files[f])
  mg.DT <- sprintf("zcat %s",fname) %>% fread
  #mg.DT <- fread("zcat /home/ob219/share/Data/GWAS-summary/MyastheniaGravis_Renton_JAMA_Neurol_2015.txt.gz")
  ## EFFECT ALLELE is ALLELE1
  mg.DT <- mg.DT[,.(trait=f,pid=paste(CHR,BP,sep=":"),a1=ALLELE1,a2=ALLELE2,or=OR,p.value=PVALUE,r2=RSQR,n1=NCASES,n0=NCONTROLS)]
  man.DT <- fread(SNP_MANIFEST)
  M <- merge(man.DT,mg.DT,by='pid')
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
  ## current effect allele is a1 we want this to be a2 so flip matches !
  flip <- which(alleles$g.class=='match')
  if(length(flip)>0)
    M[flip,or:=1/or]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  tra <- unique(M$trait)
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",tra))
  saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
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
},mc.cores=8)

saveRDS(do.call('rbind',all.mg),file=OUT_FILE)
all.mg <- do.call('rbind',all.mg)
mg.DT <- data.table(trait=rownames(all.mg),all.mg)
mg.DT <- melt(mg.DT,id.vars='trait')
control.DT <- data.table(rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='V1')
control.DT <- control.DT[V1=='control',.(pc=variable,control.score=value)]
mg.DT <- merge(mg.DT,control.DT,by.x='variable',by.y='pc')
mg.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
library(cowplot)
ggplot(mg.DT,aes(x=variable,y=value-control.score,color=trait,group=trait)) + geom_point() + geom_line() +
geom_hline(yintercept=0,lty=2)

## howmany missing SNPs ?
lapply(names(files),function(f){
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",f))
  dat <- readRDS(pfile)
  sprintf("%s missing %d",f,is.na(dat$or) %>% sum)
}) %>% do.call('c',.)
