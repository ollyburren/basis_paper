library(cupcake)

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'

man.DT <- fread(SNP_MANIFEST_FILE)

## load data
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
## compute various shrinkage methods and store

## get rotations

pc.emp <- readRDS(BASIS_FILE)
rot <- data.table(rownames(pc.emp$rotation),pc.emp$rotation)
pidf <- rot[abs(PC3)>quantile(abs(PC3),0.998),]$V1

fb<-basis.DT[pid %in% pidf,]
fb[,mlps:=-log10(p.value) * sign(log(or))]
library(ggjoy)

ggplot(fb,aes(x=mlps,y=trait)) + geom_histogram_ridges() + coord_cartesian(xlim=c(-10,10))

## alternative approach is to get t1d and ms biobank data and perform rotations before
## summing to see what is happening

basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE,trait=c('bb_MS','bb_asthma','bb_T1D'))
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
shrink.DT <- readRDS(SHRINKAGE_FILE)
SHRINKAGE_METHOD<-'ws_emp_shrinkage'
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## apply loadings manually
predict(pc.emp,basis.mat.exp,retx=TRqUE)


foo <- data.table(pid=colnames(basis.mat.emp),t(basis.mat.emp))
M <- merge(foo,rot[,.(pid=V1,PC3)],by='pid')
M[,c('bb_ms_score','bb_asthma_score','bb_t1d_score'):=list(bb_MS * PC3,bb_asthma * PC3,bb_T1D * PC3)]
M[,use:=abs(PC3)>quantile(PC3,0.80)]
tM <- melt(M[,.(pid,bb_ms_score,bb_asthma_score,bb_t1d_score)],id.vars=c('pid'))

tM[value<0,list(total=sum(value)),by=variable]

ggplot(tM,aes(x=value,y=variable)) + geom_density_ridges() + coord_cartesian(xlim=c(-0.001,0.001))
