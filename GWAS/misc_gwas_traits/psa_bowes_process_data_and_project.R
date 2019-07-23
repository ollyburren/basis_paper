## Crohns prognosis GWAS
library(annotSnpStats)

SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
DATA.DIR <- '/home/ob219/share/Data/GWAS-summary/aav_limy_wong'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_FILE <- "/home/ob219/share/as_basis/GWAS/psa_projections/summary/bowes_psa_0619.RDS"
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'

PsA.dir <- '/home/ob219/share/Data/GWAS-summary/psa-2019-unpublished'
psa.files <- list.files(path=PsA.dir,pattern="*.out",full.names=TRUE)
DT <- lapply(psa.files,fread) %>% rbindlist
DT[,ctrl_a1_maf:=(controls_AA*2 + controls_AB)/(2*controls_total)]
DT[,list(ncases=(cases_AA + cases_AB + cases_BB)/2,nctrl=(controls_AA + controls_AB + controls_BB)/2)]
DT.f <- DT[,.(rsid,chr=chromosome,pos=position,a1=alleleA,a2=alleleB,a1_maf=ctrl_a1_maf,p.value=frequentist_add_pvalue,beta=frequentist_add_beta_1,se.beta=frequentist_add_se_1)]
DT.f[,pid:=paste(chr,pos,sep=':')]
man.DT <- fread(SNP_MANIFEST)
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
pc.emp <- readRDS(BASIS_FILE)
## counted allele appears to be a2 so should be aligned but check
M <- merge(DT.f[pid %in% man.DT$pid,.(trait='bowes_psa',pid=pid,a1,a2,or=exp(beta),p.value)],man.DT,by='pid',all.y=TRUE)
idx <- which(is.na(M$or))
sprintf("%d missing",length(idx)) %>% message
if(length(idx)!=0)
  M[idx,c('trait','a1','a2','or'):=list('bowes_psa',ref_a1,ref_a2,or=1)]
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
pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",tmp$trait %>% unique))
saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
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
psa.proj <- predict(pc.emp,newdata=mat.emp)
psa.DT <- data.table(trait=rownames(psa.proj),psa.proj)
saveRDS(psa.DT,file=OUT_FILE)
