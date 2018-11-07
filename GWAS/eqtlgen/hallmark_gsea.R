library(biomaRt)

OUT_DIR <- '/home/ob219/share/as_basis/GWAS/eqtlgen_projections'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
fs <- list.files(path=OUT_DIR,pattern="*.RDS",full.names=TRUE)
res.DT <- lapply(fs,readRDS) %>% rbindlist
## define a Z score for loading to see if any are significant across pc's
M <- melt(res.DT,id.vars='trait')
M[,Z:=(value-mean(value))/sqrt(var(value)),by='variable']

ensembl = useEnsembl(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','entrezgene'),
filters = 'ensembl_gene_id', values = M$trait, mart =ensembl) %>% data.table

M <- merge(M,genedesc,by.x='trait',by.y='ensembl_gene_id')

dupent <-  M[which(duplicated(M[,.(variable,entrezgene)])),]$entrezgene %>% unique
Mf <- M[!entrezgene %in% dupent & gene_biotype=='protein_coding',]
universe <- M$entrezgene %>% as.integer

#PATHWAY_FILE <- "~/tmp/h.all.v6.1.entrez.gmt"
PATHWAY_FILE <- "~/tmp/c1.all.v6.2.entrez.gmt"


f <- scan(PATHWAY_FILE,"character()",sep="\n") %>% strsplit(.,"\t")
hm <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(hm) <- sapply(f,'[[',1)
## remove genesets with no members !
hm <- hm[sapply(hm,length)!=0]

library(parallel)
# gs.res <- mclapply(paste0('PC',c(1:10)),function(pcq){
#   bpc <- Mf[variable==pcq,]
#   pcr <- bpc$Z
#   names(pcr) <- bpc$entrezgene
#   pall <- lapply(hm,function(gs){
#     message(length(gs))
#     inset <- which(names(pcr) %in% gs)
#     W <- wilcox.test(pcr[inset], pcr[-inset],alternative="greater")
#     data.table(stat=W$statistic,p=W$p.value)
#   }) %>% rbindlist
#   data.table(pc=pcq,gs=names(hm),p=pall$p,stat=pall$stat)
# },mc.cores=6) %>% rbindlist
#
# gs.res[,adj.p:=p.adjust(p,method="fdr")]
# mat <- melt(gs.res,id.vars=c('pc','gs'),measure.vars='adj.p') %>% dcast(pc~gs+variable)
# mat <- as.matrix(mat[,-1],rownames=mat$pc)
# mat <- -log10(mat)
# library(pheatmap)
# cnames <- gsub("\\_adj\\.p|HALLMARK_","",colnames(mat)) %>% gsub("\\_"," ",.)
# colnames(mat) <- cnames
# pheatmap(mat,colnames=cnames)
hm <- hm[sapply(hm,length)>5]
gs.res.t <- mclapply(paste0('PC',c(1:10)),function(pcq){
  bpc <- Mf[variable==pcq,]
  pcr <- bpc$Z
  names(pcr) <- bpc$entrezgene
  pall <- lapply(hm,function(gs){
    inset <- which(names(pcr) %in% gs)
    K <- t.test(pcr[inset], pcr[-inset])
    data.table(stat=K$statistic,p=K$p.value)
  }) %>% rbindlist
  data.table(pc=pcq,gs=names(hm),p=pall$p,stat=pall$stat)
},mc.cores=6) %>% rbindlist
gs.res.t[,p.adj:=p.adjust(gs.res.t$p,method="fdr")]
## or
gs.res.t <- gs.res.t[p.adj<0.05,]
#gs.res.t<-gs.res.t[gs %in% (gs.res.t[p.adj<0.05,]$gs %>% unique),]

## if we look by chr region then we get enrichment FDR<0.05 at only chr19q13 containing FUT2 for PC10 - just p.adj 0.03 !


## remove all classes that don't have at least one significant result

# Z scores
gs.res.t[,as:=qnorm(p.adj/2,lower.tail=FALSE) * sign(stat)]
gs.res.t[,as:=-log10(p.adj) * sign(stat)]


mat <- melt(gs.res.t,id.vars=c('pc','gs'),measure.vars='as') %>% dcast(pc~gs+variable,fill=0)
mat <- as.matrix(mat[,-1],rownames=mat$pc)
library(pheatmap)
cnames <- gsub("\\_as|HALLMARK_","",colnames(mat)) %>% gsub("\\_"," ",.)
colnames(mat) <- cnames
pheatmap(mat,colnames=cnames)





M[variable=='PC3' & gene_biotype=='protein_coding',][order(abs(Z),decreasing=TRUE)]
