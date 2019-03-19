
## this is the set of variants for which there is a SNP in the basis I think !
f.DT<-readRDS("/home/ob219/share/as_basis/GWAS/eqtlgen/cis-eQTLs_full_20180905.filtered.RDS")
fi.DT <- f.DT[FDR<0.01,]
## add in trans eqtl
DT.t<-fread("/home/ob219/share/Data/expr/eqtlgen/trans-eQTLs_full_20180905.txt")
DT.t <- DT.t[FDR<0.01,]

DT.all <- rbind(DT.t,fi.DT)
DT.all <- DT.all[,pid:=paste(SNPChr,SNPPos,sep=':')]
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
snp.DT <- fread(SNP_MANIFEST)
DT.all <- DT.all[pid %in% snp.DT$pid,]
## how do we call whether a variant has a significant loading on a given PC.
## are the normally distributed - if so could use Z score method ?
## if not could use absolute rank method
DT.all[,raw.p:=pnorm(abs(Zscore),lower.tail=FALSE)*2,]

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
rot <- pc.emp$rotation
rot  <-  data.table(pid=rownames(rot),rot)
rot <- melt(rot,id.vars='pid')

## for cfdr we used

rot[,unq:=abs(value)]
rot[,use:=unq>quantile(unq,0.998) ,by=variable]

## turn these into genes !

M <- merge(DT.all[,.(pid,Gene,GeneSymbol,raw.p,FDR)],rot,by='pid',allow.cartesian=TRUE)
## foreach component compute fdr
M[use==TRUE,fdr.basis:=p.adjust(raw.p,method="BH"),by=variable]

## this on the gene basis so ditch pids,
by.gene <- M[,list(scomp=any(use)),by=c('Gene','GeneSymbol','variable')]
by.gene$Gene %>% unique

## do pathway analysis on this
library(biomaRt)
ensembl = useEnsembl(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','entrezgene','chromosome_name'),
filters = 'ensembl_gene_id', values = by.gene$Gene %>% unique, mart =ensembl) %>% data.table
genedesc <- genedesc[chromosome_name %in% c(1:22,'X'),]
genedesc <- genedesc[genedesc[ , .I[which.min(entrezgene)], by = 'ensembl_gene_id']$V1]

M <- merge(by.gene,genedesc,by.x='Gene',by.y='ensembl_gene_id')
Mf <- M[gene_biotype=='protein_coding' & !chromosome_name %in% c('X','Y'),]




universe <- Mf$entrezgene %>% unique %>% as.integer

## molsigdb hallmark analysis
f <- scan("~/tmp/h.all.v6.1.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
hm <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(hm) <- sapply(f,'[[',1)
## remove genesets with no members !
hm <- hm[sapply(hm,length)!=0]

## for a given pc (try pc1 to start do fisher test on two by two)

mPC <- Mf[variable=='PC1',]
ft.res.hm <- lapply(seq_along(hm),function(i){
  ft <- table(Mf$entrezgene %in% hm[[i]],Mf$scomp) %>% fisher.test
  data.table(pathway=names(hm)[i],p=ft$p.value,or=ft$estimate)
}) %>% rbindlist

## do the same again for Reactome

f <- scan("~/tmp/c2.cp.reactome.v6.2.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
react <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(react) <- sapply(f,'[[',1)
## remove genesets with no members !
react <- react[sapply(react,length) > 5]

library(parallel)
all.pc <- lapply(paste0('PC',1:10),function(pc){
  message(pc)
  mPC <- Mf[variable==pc,]
  ft.res.react <- mclapply(seq_along(react),function(i){
    ft <- table(mPC$entrezgene %in% react[[i]],mPC$scomp) %>% fisher.test
    data.table(pc=pc,pathway=names(react)[i],p=ft$p.value,or=ft$estimate)
  },mc.cores=8) %>% rbindlist
}) %>% rbindlist

all.pc[,p.adj:=p.adjust(p,method="bonferroni"),by='pc']
all.sig.pc <- all.pc[pathway %in% all.pc[p.adj<0.01,]$pathway,]
library(stringr)
all.sig.pc[,pathway_lab:=gsub("^REACTOME\\_[ ]*","",pathway) %>% str_trunc(.,30)]
all.sig.pc[,pc:=factor(pc,levels=paste0('PC',1:10))]
library(cowplot)
ggplot(all.sig.pc,aes(x=pc,y=pathway_lab,fill=-log10(p.adj))) + geom_tile()


fre.pc1 <- ft.res.react[p.adjust(p,method="bonferroni")<0.05 & or>0,]
