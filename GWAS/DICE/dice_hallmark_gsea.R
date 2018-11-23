library(biomaRt)
library(pheatmap)
library(cowplot)
library(parallel)
library(grid)
library(gridExtra)

OUT_DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/DICE_projections/CD4_STIM'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
fs <- list.files(path=OUT_DIR,pattern="*.RDS",full.names=TRUE)
res.DT <- lapply(fs,readRDS) %>% rbindlist
## define a Z score for loading to see if any are significant across pc's
M <- melt(res.DT,id.vars='trait')
M[,Z:=(value-mean(value))/sqrt(var(value)),by='variable']
ensembl = useEnsembl(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','entrezgene','chromosome_name'),
filters = 'ensembl_gene_id', values = M$trait, mart =ensembl) %>% data.table
M <- merge(M,genedesc,by.x='trait',by.y='ensembl_gene_id')

dupent <-  M[which(duplicated(M[,.(variable,entrezgene)])),]$entrezgene %>% unique
Mf <- M[!entrezgene %in% dupent & gene_biotype=='protein_coding' & !chromosome_name %in% c('X','Y'),]
universe <- Mf$entrezgene %>% as.integer

if(FALSE){
  basis <- readRDS(BASIS_FILE)
  basis.DT <- data.table(trait=rownames(basis$x),basis$x)
  basis.DT <- melt(basis.DT,id.vars="trait")
  Mf <- merge(Mf,basis.DT[trait=='control',][,.(variable,ctrl.value=value)],by='variable')
  ggplot(Mf,aes(x=variable,y=value-ctrl.value,group=trait)) + geom_line(alpha=0.1) +
  xlab("Principal Component") + ylab(expression(Delta~"Control Loading")) +
  ggtitle(expression("eQTLGen consortium cis and trans"~beta*"s")) +
  geom_hline(yintercept=0,colour="white")
  dev.print(pdf,"~/tmp/eqtlgen.pdf")
  #geom_line(data=basis.DT,aes(colour=trait),size=1,alpha=0.5)
}

gsea <- function(gl,geneid='entrezgene',test.type=c('t.test','f.test','wilcox.test'),fdr.method=c('fdr','bonferroni'),fdr.thresh=0.05){
  gl <- gl[sapply(gl,length)>5]
  gs.res.t <- mclapply(paste0('PC',c(1:10)),function(pcq){
    bpc <- Mf[variable==pcq,]
    pcr <- bpc$Z
    #pcr <- bpc$Z^2 # use chi.squared distribution
    names(pcr) <- bpc[[geneid]]
    pall <- lapply(gl,function(gs){
      inset <- which(names(pcr) %in% gs)
      if(test.type=='t.test'){
        K <- t.test(pcr[inset], pcr[-inset])
        ret.DT <- data.table(stat=K$statistic,p=K$p.value,sign=sign(K$statistic))
      }else if(test.type=='wilcox.test'){
        ## test of abs deviation so chi.squared
        K <- wilcox.test(pcr[inset]^2, pcr[-inset]^2)
        ret.DT <- data.table(stat=K$statistic,p=K$p.value,sign=sign(K$statistic))
      }else if(test.type=='f.test'){
        ## test of variance
        K <- var.test(pcr[inset], pcr[-inset])
        ret.DT <- data.table(stat=K$statistic,p=K$p.value,sign=log(K$statistic) %>% sign)
      }else{
        message("not recognised")
      }
      ret.DT
    }) %>% rbindlist
    data.table(pc=pcq,gs=names(gl),p=pall$p,stat=pall$stat,sign=pall$sign)
  },mc.cores=6) %>% rbindlist
  min.p <- gs.res.t[p!=0,]$p %>% min
  gs.res.t[p==0,p:=min.p]
  gs.res.t[,p.adj:=p.adjust(gs.res.t$p,method=fdr.method)]
  gs.res.t <- gs.res.t[p.adj<fdr.thresh,]
  gs.res.t[,as:=-log10(p.adj) * sign]
  gs.res.t
}

phwrap <- function(DT,main){
  mat <- melt(DT,id.vars=c('pc','gs'),measure.vars='as') %>% dcast(pc~gs+variable,fill=0)
  mat <- as.matrix(mat[,-1],rownames=mat$pc)
  cnames <- gsub("\\_as|HALLMARK_|chr","",colnames(mat)) %>% gsub("\\_"," ",.)
  colnames(mat) <- cnames
  pheatmap(mat,colnames=cnames,main=main,cluster_rows=FALSE,fontsize=9)[[4]]

}

## molsigdb hallmark analysis
f <- scan("~/tmp/h.all.v6.1.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
hm <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(hm) <- sapply(f,'[[',1)
## remove genesets with no members !
hm <- hm[sapply(hm,length)!=0]
## molsigdb location
f <- scan("~/tmp/c1.all.v6.2.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
lom <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(lom) <- sapply(f,'[[',1)
## remove genesets with no members !
lom <- lom[sapply(hm,length) > 5]

pl <- list()
pl[[1]] <- gsea(hm,test.type='f.test',fdr.method="fdr") %>% phwrap(.,main="F Test Hallmark")
pl[[2]] <- gsea(lom,test.type='f.test',fdr.method="fdr",fdr.thresh=0.01) %>% phwrap(.,main="F Test Location")
#pl[[2]] <- gsea(hm,test.type='t.test',fdr.method="fdr") %>% phwrap(.,main="T Test Hallmark")
g<-grid.arrange(grobs=pl,ncol=2)
dev.print(pdf,"~/tmp/hallmark_f_eqtlgen.pdf")

## do mean
gsea(hm,test.type='t.test',fdr.method="fdr") %>% phwrap(.,main="t Test Hallmark")
dev.print(pdf,"~/tmp/hallmark_t_eqtlgen.pdf")
library(xtable)
gsea(lom,test.type='t.test',fdr.method="fdr") %>% xtable


## next perform the same analysis using the t-cell modules from burren et al.

mod.DT <- fread("/home/ob219/share/as_basis/tcell_modules_burren_et_al/13059_2017_1285_MOESM1_ESM")
mod.DT.f <- mod.DT[id %in% Mf$trait,]

Mf <- Mf[trait %in% mod.DT.f$id,]
wgcna <- split(mod.DT.f$id,mod.DT.f$all_12)
pl <- list()
pl[[1]] <- gsea(wgcna,geneid='trait',test.type='f.test',fdr.method="bonferroni") %>% phwrap(.,main="F Test Burren et al")
pl[[2]] <- gsea(wgcna,geneid='trait',test.type='t.test',fdr.method="bonferroni") %>% phwrap(.,main="t Test Burren et al")
g<-grid.arrange(grobs=pl,ncol=2)
dev.print(pdf,"~/tmp/burren_eqtlgen.pdf")

gs.res.t <- mclapply(paste0('PC',c(1:10)),function(pcq){
  bpc <- Mf[variable==pcq,]
  pcr <- bpc$Z
  #pcr <- bpc$Z^2 # use chi.squared distribution
  names(pcr) <- bpc$trait
  pall <- lapply(hm,function(gs){
    inset <- which(names(pcr) %in% gs)
    K <- t.test(pcr[inset], pcr[-inset])
    #K <- wilcox.test(pcr[inset], pcr[-inset])
    #K <- var.test(pcr[inset], pcr[-inset])
    data.table(stat=K$statistic,p=K$p.value)
  }) %>% rbindlist
  data.table(pc=pcq,gs=names(hm),p=pall$p,stat=pall$stat)
},mc.cores=6) %>% rbindlist
min.p <- gs.res.t[p!=0,]$p %>% min
gs.res.t[p==0,p:=min.p]
gs.res.t[,p.adj:=p.adjust(gs.res.t$p)]
## or
gs.res.t <- gs.res.t[p.adj<0.05,]
#gs.res.t[,as:=-log10(p.adj) * sign(log(stat))] # for F test
gs.res.t[,as:=-log10(p.adj) * sign(stat)]
mat <- melt(gs.res.t,id.vars=c('pc','gs'),measure.vars='as') %>% dcast(pc~gs+variable,fill=0)
mat <- as.matrix(mat[,-1],rownames=mat$pc)
library(pheatmap)
cnames <- gsub("\\_as|HALLMARK_","",colnames(mat)) %>% gsub("\\_"," ",.)
colnames(mat) <- cnames
pheatmap(mat,colnames=cnames)





M[variable=='PC3' & gene_biotype=='protein_coding',][order(abs(Z),decreasing=TRUE)]
