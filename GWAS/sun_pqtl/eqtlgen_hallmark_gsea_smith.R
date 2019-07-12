library(biomaRt)
library(pheatmap)
library(cowplot)
library(parallel)
library(grid)
library(gridExtra)

#OUT_DIR <- '~/share/as_basis/GWAS/eqtlgen_projections_0619_softthresh'
OUT_DIR <- '~/share/as_basis/GWAS/eqtlgen_projections_significant_only_0619'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/noweight_basis_gwas_0619.RDS'
fs <- list.files(path=OUT_DIR,pattern="*.RDS",full.names=TRUE)
res <- lapply(fs,readRDS) %>% do.call('rbind',.)
res.DT <- data.table(trait=rownames(res),res)
res.DT[,hgnc:=gsub("([^\\.]+)\\..*","\\1",trait)]
res.DT$trait <- NULL
## define a Z score for loading to see if any are significant across pc's

M <- melt(res.DT[,-1],id.vars='hgnc')
M[,Z:=(value-mean(value))/sqrt(var(value)),by='variable']
M[,delta:=(value-mean(value)),by='variable']

#control <- readRDS(BASIS_FILE)$x["control",]
#control <- data.table(pc=names(control),cscore=control)
#M <- merge(M,control,by.x='variable',by.y='pc')
#M[,delta:=value-cscore]
#M[,Z:=(delta-mean(delta))/sqrt(var(delta)),by='variable']
M <- merge(unique(res.DT[,.(ensg=trait,hgnc)]),M,by='hgnc')
ensembl = useEnsembl(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','entrezgene_id','chromosome_name'),
filters = 'ensembl_gene_id', values = M$ensg, mart =ensembl) %>% data.table
genedesc <- genedesc[chromosome_name %in% c(1:22,'X'),]
genedesc <- genedesc[genedesc[ , .I[which.min(entrezgene_id)], by = 'ensembl_gene_id']$V1]

M <- merge(M,genedesc,by.x='ensg',by.y='ensembl_gene_id')

Mf <- M[gene_biotype=='protein_coding' & !chromosome_name %in% c('X','Y'),]
Mf[,raw.p:=pnorm(abs(Z),lower.tail=FALSE) * 2]
Mf[,p.adj:=p.adjust(raw.p,method="fdr"),by='variable']
eqtlgen <- copy(Mf)
eqtlgen[,source:='eqtlgen']
pc1 <- eqtlgen[source=='eqtlgen' & p.adj<0.01 & variable=='PC1',]
pc1[,mlp:=-log10(p.adj)]
pc1[,gene_name:=factor(external_gene_name,levels=pc1[order(mlp),]$external_gene_name)]
pc1 <- pc1[order(mlp),]

plot.me <- tail(pc1,n=20)

#ggplot(pc1[delta>0.005 | delta< -0.005],aes(x=gene_name,y=delta)) + geom_bar(stat="identity") + coord_flip()

pp1 <- ggplot(plot.me,aes(x=gene_name,y=mlp)) + geom_bar(stat="identity") + coord_flip() +
xlab("Gene") + ylab("-log10(p)") + ggtitle("Vosa et al. pQTL")

## for talk plot significant genes on pc1


### do pQTL stuff
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/sun_pqtl/fdr_0.05_by_chr_filtered'
fs <- list.files(path=OUT_DIR,pattern="*.RDS",full.names=TRUE)
res <- lapply(fs,readRDS) %>% do.call('rbind',.)
res.DT <- data.table(trait=rownames(res),res)
## define a Z score for loading to see if any are significant across pc's
M <- melt(res.DT,id.vars='trait')
M[,Z:=(value-mean(value))/sqrt(var(value)),by='variable']
M[,hgnc:=gsub("([^\\.]+)\\..*","\\1",trait)]
## compile a list of proteins with significant loadings above fdr threshold
ensembl = useEnsembl(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','entrezgene_id','chromosome_name'),
filters = 'external_gene_name', values = M$hgnc, mart =ensembl) %>% data.table
genedesc <- genedesc[chromosome_name %in% c(1:22,'X'),]
genedesc <- genedesc[genedesc[ , .I[which.min(entrezgene_id)], by = 'ensembl_gene_id']$V1]

M <- merge(M,genedesc,by.x='hgnc',by.y='external_gene_name')
Mf <- M[gene_biotype=='protein_coding' & !chromosome_name %in% c('X','Y'),]
Mf[,raw.p:=pnorm(abs(Z),lower.tail=FALSE) * 2]
Mf[,p.adj:=p.adjust(raw.p,method="fdr"),by='variable']
pqtl <- Mf[,source:='pqtl']


pc1p <- pqtl[source=='pqtl' & p.adj<0.05 & variable=='PC1',]
pc1p[,mlp:=-log10(p.adj)]
pc1p[,gene_name:=factor(hgnc,levels=pc1p[order(mlp),]$hgnc)]
pc1p <- pc1p[order(mlp),]

plot.mep <- tail(pc1p,n=20)

#ggplot(pc1[delta>0.005 | delta< -0.005],aes(x=gene_name,y=delta)) + geom_bar(stat="identity") + coord_flip()

pp2 <- ggplot(plot.mep,aes(x=gene_name,y=mlp)) + geom_bar(stat="identity") + coord_flip()  +
xlab("Gene") + ylab("-log10(p)") + ggtitle("Sun et al. pQTL")


plot_grid(pp1,pp2) %>% save_plot(file="~/tmp/smith_qtl_analysis.pdf",.,base_width=11)


pqtl <- pqtl[pqtl[ , .I[which.max(abs(Z))], by = c('ensembl_gene_id','variable')]$V1]


#foo<-pqtl[,list(count=.N),by=c('ensembl_gene_id','variable')]
## get rid of duplicated genes across pcs
#pqtl <- pqtl[ensembl_gene_id %in% foo[count==1,]$ensembl_gene_id ,]


all <- rbind(pqtl[,.(ensg=ensembl_gene_id,external_gene_name=hgnc,p.adj,pc=variable,Z,source)],eqtlgen[,.(ensg,external_gene_name,p.adj,pc=variable,Z,source)])
melted <- melt(all[p.adj<0.1,],id.vars=c('ensg','external_gene_name','pc','source'),measure.var='Z')
## get rid of those in pQTL twice

foo<-melted[,list(count=.N),by=c('ensg','pc')]
melted <- melted[ensg %in% foo[count==2,]$ensg,]

dat <- dcast(melted,external_gene_name~pc+source+variable)



keep.genes <- Mf[variable %in% c('PC1','PC4','PC7','PC11') & p.adj<0.001,]$external_gene_name

melted <- melt(Mf[variable %in% c('PC1','PC4','PC7','PC11') & external_gene_name %in% keep.genes,.(external_gene_name,pc=variable,Z)],id.vars=c('external_gene_name','pc'))
mat1 <- dcast(melted,external_gene_name~pc)
mat <- mat1[,-1] %>% as.matrix
rownames(mat) <- mat1$external_gene_name
pcs <- 'PC1'
keep.genes <- lapply(c('PC1','PC4','PC7','PC11'),function(pcs){
  tmp <- melted[pc==pcs,][order(value),]
  rbind(head(tmp,n=10),tail(tmp,n=10))$external_gene_name
}) %>% do.call('c',.) %>% unique

melted <- melt(Mf[variable %in% c('PC1','PC4','PC7','PC11') & external_gene_name %in% keep.genes,.(external_gene_name,pc=variable,Z)],id.vars=c('external_gene_name','pc'))
mat1 <- dcast(melted,external_gene_name~pc)
mat <- mat1[,-1] %>% as.matrix
rownames(mat) <- mat1$external_gene_name
#mat[(pnorm(abs(mat),lower.tail=FALSE) * 2 %>% p.adjust)>0.01]<-NA

pheatmap(mat[,1] %>% order,cluster_rows=FALSE,cluster_cols=FALSE)
mat<-mat[,c('PC1','PC4','PC7','PC11')]
pheatmap(mat,cluster_rows=FALSE,cluster_cols=FALSE)


#Mf[,Zd:=(delta-mean(delta))/sqrt(var(delta)),by='variable']
#Mf[,Z:=(value-mean(value))/sqrt(var(value)),by='variable']

#setnames(Mf,'hgnc','trait')

## there are duplicates so get rid of these
 idx <- Mf[,list(idx=which(duplicated(Mf,by=c("variable","hgnc"))))]$idx
 ## save duplicates here but dont use there are 265
 dup <- Mf[idx,]
 Mf <- Mf[-idx,]




gp <- melt(Mf[hgnc %in% Mf[p.adj<0.05,]$hgnc,.(pc=variable,delta,gene=hgnc)],measure.vars='delta') %>% dcast(.,gene~pc)
mat <- as.matrix(gp[,-1])
rownames(mat) <- gp$gene
hc <- dist(mat) %>% hclust(.)

library(cowplot)
Mf <- Mf[variable !='PC11',]
tp <- Mf[hgnc %in% Mf[p.adj<0.05,]$hgnc,.(pc=variable,delta,gene=hgnc,p.adj)]
tp[,is.sig:='']
tp[p.adj<0.05,is.sig:='*']
tp[,pc:=factor(pc,levels=paste('PC',1:10,sep=''))]
tp[,gene:=factor(gene,levels=hc$labels[hc$order])]
bbplot <- ggplot(tp[pc!='PC11'],aes(x=pc,y=gene,fill=delta,label=is.sig))  +
geom_tile() + geom_text() + scale_fill_gradient2("Difference\nfrom control") + xlab("Principal Component") +
ylab("Gene") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
save_plot(bbplot,file="~/tmp/pqtl_genes.pdf",base_height=8,base_aspect_ratio=0.9)

## bogstandard analysis



universe <- Mf$entrezgene %>% as.integer

## molsigdb hallmark analysis
f <- scan("~/tmp/h.all.v6.1.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
hm <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(hm) <- sapply(f,'[[',1)
## remove genesets with no members !
hm <- hm[sapply(hm,length)>20]
## molsigdb location
f <- scan("~/tmp/c1.all.v6.2.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
lom <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(lom) <- sapply(f,'[[',1)
## remove genesets with no members !
lom <- lom[sapply(hm,length) > 20]
## reactome
f <- scan("~/tmp/c2.cp.reactome.v6.2.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
react <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(react) <- sapply(f,'[[',1)
## remove genesets with no members !
react <- react[sapply(react,length) > 20]

f <- scan("~/tmp/c2.cp.kegg.v6.2.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
kegg <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(kegg) <- sapply(f,'[[',1)
## remove genesets with no members !
kegg <- kegg[sapply(kegg,length) > 20]

## plot for talk
newgsea <- function(gl){
  gl <- gl[sapply(gl,length)>20]
  gs.res.t <- mclapply(paste0('PC',c(1:10)),function(pcq){
    bpc <- Mf[variable==pcq,]
    pcr <- bpc$Z
    #pcr <- bpc$Z^2 # use chi.squared distribution
    names(pcr) <- bpc[['entrezgene']]
    pall <- lapply(gl,function(gs){
      inset <- which(names(pcr) %in% gs)
      K <- var.test(pcr[inset], pcr[-inset])
      data.table(stat=K$statistic,p=K$p.value,sign=log(K$statistic) %>% sign,lo=K$conf.int[1],hi=K$conf.int[2])
    }) %>% rbindlist
    data.table(pc=pcq,gs=names(gl),p=pall$p,stat=log(pall$stat),lo=log(pall$lo),hi=log(pall$hi),sign=pall$sign)
  },mc.cores=6) %>% rbindlist
  gs.res.t
}


#gs.res.t[,p.adj:=p.adjust(gs.res.t$p,method='fdr')]
#gs.res.t <- gs.res.t[p.adj<fdr.thresh,]
#gs.res.t[,as:=-log10(p.adj) * sign]

gsea.react.DT <- newgsea(react)
gsea.hm.DT <- newgsea(hm)
gsea.kegg.DT <- newgsea(kegg)
gsea.DT <- rbindlist(list(gsea.react.DT,gsea.hm.DT,gsea.kegg.DT))
#gsea.DT <- gsea.react.DT
## test to see what happens if we just look at one PC
gsea.DT[,p.adj:=p.adjust(p,method='BH')]
## keep KEGG, REACTOME and

#gsea.DT <- gsea.hm.DT

## for thesis show all components

filt.DT <- gsea.DT[gs %in% gsea.DT[p.adj<0.01 & sign >0,]$gs,]
library(stringr)

#filt.DT <- gsea.DT[pc %in% c('PC1','PC3')  & p.adj<0.05 & sign>0,]
filt.DT[,source:='REACTOME']
filt.DT[grep('^HALLMARK_',gs),source:='HALLMARK']
filt.DT[grep('^KEGG_',gs),source:='KEGG']
filt.DT[,gs:=gsub("(HALLMARK|REACTOME|KEGG)\\_","",gs) %>% gsub("\\_"," ",.) %>% str_trunc(., 35, "right")]
filt.DT[,gs:=factor(gs,levels=filt.DT[order(stat,decreasing=FALSE),]$gs %>% unique)]
filt.DT[,pc:=factor(pc,levels=paste0('PC',1:10,sep=''))]
filt.DT[,is.sig:='']
filt.DT[p.adj<0.01,is.sig:='*']
#filt.DT <- filt.DT[source=='REACTOME']


mat <- melt(filt.DT,id.vars=c('pc','gs'),measure.vars='stat') %>% dcast(gs~pc)
m <- as.matrix(mat[,-1])
rownames(m) <- mat$gs
hco <- dist(m) %>% hclust
filt.DT[,gs:=factor(gs,level=hco$labels[hco$order])]
pp <- ggplot(filt.DT,aes(y=gs,x=pc,fill=stat,label=is.sig)) + geom_tile(color='grey') + geom_text() +
xlab("Component") + ylab("Pathway") + #ggtitle(sprintf("V\u00F5sa et al. eQTLGEN PC3 FDR<0.01")) +
theme(axis.text.y = element_text(angle = 0, vjust = 1, hjust=1,size=10)) + scale_fill_gradient2("log(F-statistic)") +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + facet_grid(source~.,scales="free_y",space="free_y")
save_plot(pp,file="~/tmp/pqtl_pw.pdf",base_height=8,base_aspect_ratio=1.1)
