## do GTEX regression analysis
library(data.table)
library(magrittr)
library(biomaRt)
library(mygene)


GTEX_WB_COEFF_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/gtex_whole_blood_coeff.RDS'
DATA_FILE <- '/home/ob219/rds/rds-cew54-wallace-share/olly/GTEX/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed'

if(!file.exist(GTEX_WB_COEFF_FILE)){
  # load in GTEX individual projections
  ind.proj <- readRDS("/home/ob219/rds/hpc-work/as_basis/support/ind_proj_june10k.RDS")[trait=='gtex',]

  ## work on whole blood to begin with

  ## these need to be extracted first from a tar file
  #tar -xf  /home/ob219/rds/rds-cew54-wallace-share/Data/gtex/GTEx_Analysis_v8_eQTL_expression_matrices.tar GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz
  #tar -xf  /home/ob219/rds/rds-cew54-wallace-share/Data/gtex/GTEx_Analysis_v8_eQTL_expression_matrices.tar GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz.tbi



  dat <- fread(DATA_FILE)
  genes <- dat$gene
  col.keep <- ind.proj$individual
  E <- dat[,col.keep,with=FALSE] %>% t
  colnames(E) <- genes

  ## next get the covariate matrx

  c <- fread('/home/ob219/rds/rds-cew54-wallace-share/Data/gtex/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt')
  C <- c[,col.keep,with=FALSE] %>% t
  colnames(C) <- c$ID

  #check everything is OK
  identical(rownames(C),ind.proj$individual)
  identical(rownames(E),ind.proj$individual)


  library(parallel)
  all.res <- lapply(paste0('PC',1:11),function(pc){
    message(sprintf("Processing %s",pc))
    DT <- cbind(ind.proj[[pc]],C) %>% data.table
    setnames(DT,c('as_basis',colnames(C)))
    ## next we do the regression of a PC onto gene expression
    mclapply(colnames(E),function(g){
      coef <- summary(lm(E[,g]~.,data=DT))$coef["as_basis",]
      data.table(pc=pc,gene=g,p.val=coef[4],t=coef[3])
    },mc.cores=8) %>% rbindlist
  }) %>% rbindlist

  saveRDS(all.res,file=GTEX_WB_COEFF_FILE)
}else{
  all.res <- readRDS(GTEX_WB_COEFF_FILE)
}

## plot qq
pd <- lapply(split(all.res[,.(pc,p.val)],all.res$pc),function(p){
  po <- -log10(sort(p$p.val))
  n <- length(po)
  q <- -log10((n:1)/(n+1)) %>% rev
  data.table(pc=p$pc,p=po,ep=q)
}) %>% rbindlist

pd[,pc:=factor(pc,levels=paste0('PC',1:11))]

library(cowplot)

ggplot(pd,aes(x=ep,y=p)) + geom_point(size=0.5) + geom_abline(slope=1,col='red',lty=2) + facet_wrap(~pc) + xlab("Expected -log10(P)") + ylab("Observed -log10(P)")

## get genes that have a 5% FDR

all.adj <- all.res[,list(gene=gene,p.val=p.val,adj.pval=p.adjust(p.val,method="fdr"),t.rank=frankv(abs(t),order=-1),t=t),by='pc']
genes.for.chris <- all.adj[t.rank<=10 & pc %in% paste0('PC',c(1:3,5:8)),][,ensg:=gsub("\\.[0-9]+","",gene)]

## ok let us annotate these genes
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','external_gene_name','entrezgene'), filters = 'ensembl_gene_id', values = genes.for.chris$ensg, mart =ensembl) %>% data.table
## use mygene to get entrez stuff
geneentrez <- getGenes(genedesc[!is.na(entrezgene),]$entrezgene, fields=c("summary")) %>% as.data.frame %>% data.table
geneentrez[,query:=as.integer(query)]
mdesc <- merge(genedesc,geneentrez[,.(query,summary)],by.x='entrezgene',by.y='query',all.x=TRUE)
M <- merge(genes.for.chris,mdesc,by.x='ensg',by.y='ensembl_gene_id',all.x=TRUE)
M <- M[order(pc,t.rank),]
library(xlsx)
write.xlsx(M, file="/home/ob219/tmp/gtex_whole_blood.xls", sheetName="all")
