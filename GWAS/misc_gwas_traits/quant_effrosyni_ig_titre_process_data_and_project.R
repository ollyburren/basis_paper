## process Ig titre
library(annotSnpStats)
library(cupcake)
library(parallel)
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/effrosyni_ig'
DATA_DIR <- '/home/ob219/share/Data/GWAS-summary/ig_titre_unpublished'
files <- list.files(path=DATA_DIR,pattern="*.txt",full.names=TRUE)
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
man.DT <- fread(SNP_MANIFEST)

ret <- mclapply(files,function(f){
i.DT <- fread(f)
i.DT[,pid:=paste(chr,pos,sep=':')]
i.DTf <- i.DT[pid %in% man.DT$pid,]
res.DT <- with(i.DTf,convertBetaToOR(N=N,b=beta,seb=se,m=maf) %>% do.call('cbind',.) %>% data.table)
i.DTf <- cbind(i.DTf[,.(pid,a1=allele_A,a2=allele_B,maf)],res.DT)
#i.DT <- i.DT[,.(pid,a1=allele_A,a2=allele_B,or=exp(beta),p.value=p_value)]
M <- merge(i.DTf,man.DT,by='pid')
M[,maf_ref:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
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
M <- M[!duplicated(pid),]
M[,trait:=gsub("LN\\_([^\\_]+).*","\\1",basename(f))]

SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_vit_t2d.RDS'
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)
tmp<-M[stmp]
tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$OR)
B <- dcast(tmp,pid ~ trait,value.var='metric',fill=0)
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_vit_t2d.RDS'
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
},mc.cores=3)
ret <- do.call('rbind',ret)
res.DT <- data.table(trait=rownames(ret),ret)
saveRDS(res.DT,file.path(OUT_DIR,'effrosyni_ig_vit_t2d.RDS'))
