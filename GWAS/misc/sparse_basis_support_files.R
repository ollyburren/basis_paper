## create a data set that we can transport with the package for building the sparse basis
## this will constitute a vignette on how to build the basis

## note that this should not be distributed with the package and is for reference only

## first load in chris' sparse basis

load("~/share/as_basis/sparse-basis/basis-sparse-13-0.999.RData")
keep.pids <- rownames(use.pca)
## first create a snp manifest file
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas_13_traits_0919.tab'
## filter down all the files above

snp.DT <- fread(SNP_MANIFEST_FILE)[pid %in% keep.pids,]
saveRDS(snp.DT,file="~/tmp/sparse_basis/13_traits_sparse_snp_manifest.RDS")
shrink.DT <- readRDS(SHRINKAGE_FILE)[pid %in% keep.pids,]
## save new format
saveRDS(shrink.DT[,.(pid,shrinkage=ws_emp_shrinkage)],file="~/tmp/13_traits_sparse_shrinkage.RDS")

## for the trait manifest import the correct files and then filter so just have the subset of snps in the manifest.

man.DT <- fread(TRAIT_MANIFEST)
man.DT[basis_trait==1,.(trait,disease,cases,controls,pmid,file)]
man.DT[,ofile:=file]

man.DT[file=='cousminer_lada.tab',ofile:='lada_cousminer.tab']
man.DT[file=='cousminer_lada.tab',trait:='LADA']
man.DT[file=='cousminer_lada.tab',pmid:='30254083']

man.DT[file=='kiryluk_neph.tab',ofile:='iga_neph_kiryluk.tab']
man.DT[file=='kiryluk_neph.tab',trait:='IgA_NEPH']
man.DT[file=='kiryluk_neph.tab',pmid:='25305756']

write.table(man.DT[,.(trait,disease,cases,controls,pmid,file=ofile)],file="~/tmp/13_traits_manifest.tab",sep="\t",quote=FALSE,row.names=FALSE)

for(i in 1:nrow(man.DT)){
  message(sprintf("Processing %s",man.DT[i,]$trait))
  tDT<-fread(sprintf("%s/%s",GWAS_DATA_DIR,man.DT[i,]$file))[pid %in% keep.pids,]
  tDT[,or:=signif(or,digits=3)]
  tDT[,p.value:=signif(p.value,digits=3)]
  write.table(tDT,file=sprintf("~/tmp/sparse_basis/%s",man.DT[i,]$ofile),quote=FALSE,sep="\t",row.names=FALSE)
}

## add a ukbb trait for projection

ukbb <- readRDS("/home/ob219/share/as_basis/GWAS/for_fdr_13_traits_0919/UKBB_NEALE:SRD:systemic.lupus.erythematosis.sle_source.RDS")
ukbb <- ukbb[pid %in% keep.pids,]
ukbb[,or:=signif(or,digits=3)]
ukbb[,p.value:=signif(p.value,digits=3)]
write.table(ukbb,file="~/tmp/sparse_basis/UKBB_NEALE:SRD:SLE.tab",quote=FALSE,sep="\t",row.names=FALSE)
