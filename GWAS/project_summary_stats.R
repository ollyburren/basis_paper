library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)

SHRINKAGE_METHOD<-'ws_emp'
SHRINKAGE_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/shrinkage_june10k.RDS'
BASIS_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/basis_june10k.RDS'
GWAS_DATA_DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/all_studies_filtered'
SNP_MANIFEST_FILE <-'/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab'
MANIFEST <- '/home/ob219/git/as_basis/manifest/as_manifest_july.tsv'
VARIANCE_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/analytical_variances_june10k.RDS'
SUMMARY_STATS_PROJ_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/summary_stat_proj_june10k.RDS'


shrink.DT <- readRDS(SHRINKAGE_FILE)
pc.emp <- readRDS(BASIS_FILE)

## project case control only
manifest.DT <- fread(MANIFEST)[include=='Y' & basis_trait==0 & pmid != '27863252' & genotype=='N',]
## something wrong with mpo_Pos
traits <- manifest.DT[trait!='mpo_Pos',]$trait
proj.DT<-get_gwas_data(MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,TRUE,traits)
proj.mat.emp<-create_ds_matrix(proj.DT,shrink.DT,SHRINKAGE_METHOD)
pred.emp <- predict(pc.emp,newdata=bb.mat.emp)
pred.DT <- data.table(trait=rownames(pred.emp),individual=NA,pred.emp)
saveRDS(pred.DT,file=SUMMARY_STATS_PROJ_FILE)
