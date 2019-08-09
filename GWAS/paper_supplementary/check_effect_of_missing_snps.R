## what about missing SNPs ?
library(cupcake)
SOURCE_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'
files <- list.files(path=SOURCE_DIR,pattern='*.RDS',full.names=TRUE)

SHRINKAGE_METHOD<-'ws_emp_shrinkage'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/missing_snp_effect'
man.DT <- fread(SNP_MANIFEST_FILE)

#f<-files[12]
library(parallel)
foo<-mclapply(files,function(f){
#for(f in files){
  trait.label <- basename(f) %>% gsub("\\_source.RDS","",.)
  sprintf("Processing %s",trait.label) %>% message
  ofile <- sprintf("%s.RDS",trait.label) %>% file.path(OUT_DIR,.)
  if(file.exists(ofile)){
    sprintf("File %s exists skipping.",ofile) %>% message
    return()
  }
  DT <- readRDS(f)
  if(!any(names(DT)=='or')){
    sprintf("Can't find or column for %s, skipping.",f) %>% message
    return()
  }
  keep.pids <- DT[!is.na(or),]$pid
  shared.basis.pids <- intersect(keep.pids,man.DT$pid)
  missing <- nrow(man.DT) - length(shared.basis.pids)
  if(length(shared.basis.pids)==nrow(man.DT)){
    message("No missing SNPs")
    return()
  }
  stitle <-  length(shared.basis.pids) %>% sprintf("%s:Retained %d/%d, missing %d (%.2f%%)",trait.label,.,nrow(man.DT),nrow(man.DT)-.,((nrow(man.DT)-.)/nrow(man.DT)*100))
  ## create a temporary manifest file
  tfile <- tempfile(pattern=sprintf("%s_",trait.label),fileext = ".tab")
  write.table(man.DT[pid %in% keep.pids,],file=tfile,row.names=FALSE,quote=FALSE,sep="\t")
  ## load data
  basis.DT<-get_gwas_data(TRAIT_MANIFEST,tfile,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
  ## compute various shrinkage methods and store
  shrink.DT<-compute_shrinkage_metrics(basis.DT)
  basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
  ## project on the data and store PC scores.
  stmp<-shrink.DT[,.(pid,ws_emp_shrinkage)]
  tmp <- merge(DT[!is.na(or),.(trait=trait.label,pid,or)],stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  basis.DT <- data.table(ptrait=trait.label,trait=rownames(pc.emp$x),source='basis',pc.emp$x)
  out.DT <- data.table(ptrait=trait.label,trait=trait.label,source='projection',all.proj)
  out.DT <- rbindlist(list(out.DT,basis.DT),fill=TRUE)
  out.DT[,missing:=nrow(man.DT)-length(shared.basis.pids)]
  ## need to keep basis pc scores also
  saveRDS(out.DT,file=ofile)
  unlink(tfile)
},mc.cores=8)


library(cowplot)
library(ggrepel)
## quick check is to look at correlation between PC scores for full basis zero added for missing data
## and cut down basis.
files <- list.files(path=OUT_DIR,pattern='*.RDS',full.names=TRUE)
tailored.basis.res.DT <- lapply(files,readRDS) %>% rbindlist
tailored.basis.res.DT <- melt(tailored.basis.res.DT,id.vars=c('ptrait','trait','source','missing'))
## compute delta
controls <- tailored.basis.res.DT[trait=='control',.(trait=ptrait,cv=value,pc=variable)]
M <- merge(tailored.basis.res.DT[,.(trait,missing,pc=variable,pv=value)],controls,by=c('trait','pc'))
M <- M[,.(trait,pc,missing,tdelta=pv-cv)]

RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/02_08_19_0619_primary_results.RDS'
full.basis.res.DT <- readRDS(RESULTS.FILE)

am <- merge(full.basis.res.DT[,.(trait,fdelta=delta,pc=variable)],M,by=c('trait','pc'))
am[,pc:=factor(pc,levels=paste0('PC',1:12))]
## quick and dirty labelling of outliers
am[,residuals:=abs(fdelta)-abs(tdelta),by=pc]
am[,Z.res:=(residuals)/sd(residuals),by=pc]
am[,label:='']
am[abs(Z.res)>1.96,label:=trait]

ggplot(am,aes(x=abs(fdelta),y=abs(tdelta),color=log10(missing),label=label)) + geom_point() + facet_wrap(~pc) +
geom_abline(col='red',lty=2) + geom_text_repel() + xlab("Full basis delta") + ylab("Tailored basis delta") + labs(color = "log10(#missing SNPs)")
