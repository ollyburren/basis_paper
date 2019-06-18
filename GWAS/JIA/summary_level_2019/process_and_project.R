library(annotSnpStats)
## code to align the latest summary results from JIA GWAS
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
ss.DT <- fread("~/share/Data/GWAS/jia-mar-2019/summary-stats-mar2019.csv")
snp.DT <- fread("~/share/Data/GWAS/jia-mar-2019/summary-stats-snpinfo-mar2019.csv")
sc.DT <- fread("~/share/Data/GWAS/jia-mar-2019/summary-stats-samplecount-mar2019.csv")
setnames(sc.DT,'n','n1')


## split into different subtypes

sub.DT <- data.table(idx=0:9,subtype=c('case','sys','PO','EO','RFneg','RFpos','ERA','PsA','undiff','missing'))
sub.DT <- merge(sub.DT,sc.DT[,.(ilar_pheno,n1)],by.x='idx',by.y='ilar_pheno')
sub.DT[,subtype:=sprintf("jia_%s_19",subtype)]
M <- merge(ss.DT,sub.DT,by.x='ilar',by.y='idx')
M[,p.value:=pnorm(abs(b)/v,lower.tail=FALSE) * 2]
snp.DT[,pid:=paste(chromosome,position,sep=':')]
jia.DT <- merge(M,snp.DT[,.(pid,a1=allele.1,a2=allele.2,snp.name)],by.x='snp',by.y='snp.name')

man.DT <- fread(SNP_MANIFEST)
M <- merge(jia.DT[,.(pid,a1,a2,or=exp(b),trait=subtype)],man.DT,by='pid')
alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
#alleles <- alleles[!duplicated(pid),]
#alleles <- M[,list(al.x=paste(uk10_A1,uk10_A2,sep='/'),al.y=paste(a1,a2,sep='/')),by='pid']
## to make quick
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
## it appears as if everything is reversed
M[,or:=1/or]
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_vit_t2d.RDS'
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)

## do this for each subtype so we can catch all missing SNPS

subtypes <- split(M,M$trait)

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_vit_t2d.RDS'
pc.emp <- readRDS(BASIS_FILE)

proj <- lapply(subtypes,function(M){
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
}) %>% do.call('rbind',.)

jia.DT <- data.table(trait=rownames(proj),proj)
saveRDS(jia.DT,file="/home/ob219/share/as_basis/GWAS/jia_projections/summary/jia_vit_t2d_2019.RDS")
jia.DT <- melt(jia.DT,id.vars='trait')

control.DT <- data.table(rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='V1')
control.DT <- control.DT[V1=='control',.(pc=variable,control.score=value)]
jia.DT <- merge(jia.DT,control.DT,by.x='variable',by.y='pc')
jia.DT[,variable:=factor(variable,levels=paste0('PC',1:13))]
library(cowplot)
ggplot(jia.DT,aes(x=variable,y=value-control.score,color=trait,group=trait)) + geom_point() + geom_line() +
geom_hline(yintercept=0,lty=2)
