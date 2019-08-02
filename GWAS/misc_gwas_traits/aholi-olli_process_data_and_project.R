library(annotSnpStats)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- "/home/ob219/share/as_basis/GWAS/ahola-olli_cytokine/projections/ck_0619.RDS"
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'
DATA_DIR <- '/home/ob219/share/Data/GWAS-summary/cytokine_gwas/'

man.DT <- fread(SNP_MANIFEST)
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
pc.emp <- readRDS(BASIS_FILE)

files <- list.files(path=DATA_DIR,pattern="*.gz",full.names=TRUE)
library(parallel)
res <- mclapply(files,function(f){
  trait <- basename(f) %>% gsub(".data.gz","",.) %>% sprintf("CK:%s",.)
  sprintf("Processing %s",trait) %>% message
  DT <- sprintf("zcat %s",f) %>% fread
  DT[,pid:=paste(Chromosome,Position,sep=':')]
  M <- merge(DT[,.(trait=trait,pid,a1=toupper(OtherAllele),a2=toupper(EffectAllele),or=exp(Effect),p.value=P.value)],man.DT,by='pid',all.y=TRUE)
  M <- M[!duplicated(pid),]
  idx <- which(is.na(M$or))
  sprintf("%d missing",length(idx)) %>% message
  if(length(idx)!=0)
    M[idx,c('trait','a1','a2','or'):=list(trait,ref_a1,ref_a2,or=1)]
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
  tmp$trait <- trait
  pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",trait))
  saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
  ## where snp is missing make it zero
  tmp[is.na(metric),metric:=0]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
    stop("Something wrong basis and projection matrix don't match")
  predict(pc.emp,newdata=mat.emp)
},mc.cores=8)

res <- do.call('rbind',res)
ck.DT <- data.table(trait=rownames(res),res)
saveRDS(ck.DT,file=OUT_FILE)
