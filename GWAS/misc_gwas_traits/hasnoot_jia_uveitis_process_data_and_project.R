## bipolar prognosis GWAS
library(annotSnpStats)
library(rtracklayer)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- '/home/ob219/share/as_basis/GWAS/hasnoot_uveitis/hasnoot_uveitis_0619.RDS'
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'
man.DT <- fread(SNP_MANIFEST)
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
pc.emp <- readRDS(BASIS_FILE)

## first perform the meta analysis


p1.dt <- fread("zcat /home/ob219/share/Data/GWAS-summary/hasnoot_jia_uveitis/jia.uveitis.hg19.phase1.QC.gwas.results.imputed.jia.uviitis.filtered.v2.txt.gz")
setnames(p1.dt,names(p1.dt) %>% paste('p1',.,sep='.'))
## check for infinite effect sizes
p1.inf <- which(!is.finite(p1.dt$p1.EFFECT))
## attempt to fill these in using p.val and se
p1.dt[p1.inf,p1.EFFECT:=2 * p1.SE * qnorm(p1.P,lower.tail=FALSE)]
p2.dt <- fread("zcat /home/ob219/share/Data/GWAS-summary/hasnoot_jia_uveitis/jia.uveitis.hg19.phase2.QC.gwas.results.imputed.jia.uviitis.filtered.v2.txt.gz")
setnames(p2.dt,names(p2.dt) %>% paste('p2',.,sep='.'))
p2.inf <- which(!is.finite(p2.dt$p2.EFFECT))
p2.dt[p2.inf,p2.EFFECT:=2 * p2.SE * qnorm(p2.P,lower.tail=FALSE)]
## combine OR using inverse variance

pall.dt <- merge(p1.dt[,.(pid=paste(p1.CHR,p1.BP,sep=':'),p1.SNP,p1.A1,p1.A2,p1.EFFECT,p1.SE)],p2.dt[,.(p2.SNP,p2.A1,p2.A2,p2.EFFECT,p2.SE)],by.x='p1.SNP',by.y='p2.SNP')
pall.dt <- pall.dt[pid %in% man.DT$pid,]
pall.dt[,c('w1','w2'):=list(1/p1.SE^2,1/p2.SE^2)]
pall.dt[,c('meta.beta','meta.se'):=list(((p1.EFFECT * w1)+(p2.EFFECT * w2))/(w1+w2),sqrt(1/(w1+w2)))]
pall.dt[,c('meta.Z','meta.p.value'):=list(meta.beta/meta.se,2 * pnorm(abs(meta.beta/meta.se),lower.tail=FALSE))]


M <- merge(pall.dt[,.(pid,a1=p1.A1,a2=p1.A2,or=exp(meta.beta),p.value=meta.p.value)],man.DT,by='pid')

##alignment will be the same for both studies
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
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
## as this is plink the test allele is a1 therefore we need to flip all OR
M[g.class=='match',or:=1/or]
tmp <- merge(M,stmp,by='pid',all.y=TRUE)
tra <- 'hasnoot_uveitis_jia'
pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",tra))
saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
tmp[,trait:=tra]
## where snp is missing make it zero
tmp[is.na(metric),metric:=0]
B <- dcast(tmp,pid ~ trait,value.var='metric')
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
uv.DT <- data.table(trait=rownames(all.proj),all.proj)
saveRDS(uv.DT,file=OUT_FILE)
