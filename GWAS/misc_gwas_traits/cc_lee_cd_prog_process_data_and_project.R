## Crohns prognosis GWAS


c.DT <- fread("zcat /home/ob219/share/Data/GWAS-summary/CD_prognosis_GWA_results.csv.gz")
c.DT[,pid:=paste(chr,pos,paste=':')]
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
man.DT <- fread(SNP_MANIFEST)
