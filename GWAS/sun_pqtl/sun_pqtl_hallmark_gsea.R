library(biomaRt)
library(pheatmap)
library(cowplot)
library(parallel)
library(grid)
library(gridExtra)

#OUT_DIR <- '/home/ob219/share/as_basis/GWAS/sun_pqtl_filter'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/sun_pqtl/fdr_0.05_by_chr_filtered'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
fs <- list.files(path=OUT_DIR,pattern="*.RDS",full.names=TRUE)
res <- lapply(fs,readRDS) %>% do.call('rbind',.)
res.DT <- data.table(trait=rownames(res),res)
## define a Z score for loading to see if any are significant across pc's
M <- melt(res.DT,id.vars='trait')
M[,Z:=(value-mean(value))/sqrt(var(value)),by='variable']
M[,hgnc:=gsub("([^\\.]+)\\..*","\\1",trait)]

## compile a list of proteins with significant loadings above fdr threshold

M[,p:=pnorm(abs(Z),lower.tail=FALSE) * 2]
M[,p.adj:=p.adjust(p,method="fdr"),by='variable']

M[p.adj<0.05,]



ensembl = useEnsembl(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','entrezgene','chromosome_name'),
filters = 'external_gene_name', values = M$hgnc, mart =ensembl) %>% data.table
genedesc <- genedesc[chromosome_name %in% c(1:22,'X'),]
genedesc <- genedesc[genedesc[ , .I[which.min(entrezgene)], by = 'ensembl_gene_id']$V1]

M <- merge(M,genedesc,by.x='hgnc',by.y='external_gene_name')

Mf <- M[gene_biotype=='protein_coding' & !chromosome_name %in% c('X','Y'),]


universe <- Mf$entrezgene %>% as.integer %>% unique

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

gsea <- function(gl,geneid='entrezgene',test.type=c('t.test','f.test','wilcox.test','bartlett.test'),fdr.method=c('fdr','bonferroni'),fdr.thresh=0.05){
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
      }else if(test.type=='bartlett.test'){
        ## test of variance
        K <- bartlett.test(list(pcr[inset],pcr[-inset]))
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
  gs.res.t[,p.adj:=p.adjust(p,method=fdr.method),by='pc']
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

pdf("~/tmp/hallmark_sun.pdf",paper="a4r")
gsea(geneid='entrezgene',hm,test.type='f.test',fdr.method="fdr",fdr.thresh=0.05) %>% phwrap(.,main="F Test Hallmark")
dev.off()


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
## reactome
f <- scan("~/tmp/c2.cp.reactome.v6.2.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
react <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(react) <- sapply(f,'[[',1)
## remove genesets with no members !
react <- react[sapply(react,length) > 5]

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
gsea.hm.DT[,p.adj:=p.adjust(p,method='fdr')]
gsea.DT <- rbind(gsea.react.DT,gsea.hm.DT)
## test to see what happens if we just look at one PC
gsea.DT[,p.adj:=p.adjust(p,method='fdr')]
gsea.DT <- gsea.DT[pc=='PC3' & p.adj<0.05 & sign>0,]
gsea.DT[,source:='REACTOME']
gsea.DT[grep('^HALLMARK_',gs),source:='HALLMARK']
gsea.DT[,gs:=gsub("(HALLMARK|REACTOME)\\_","",gs) %>% gsub("\\_"," ",.) %>% str_trunc(., 30, "right")]
gsea.DT[,gs:=factor(gs,levels=gsea.DT[order(stat,decreasing=TRUE),]$gs)]

pp <- ggplot(gsea.DT,aes(y=stat,x=gs,color=source)) + geom_point() + geom_errorbar(aes(ymin=lo,ymax=hi)) + coord_flip() +
xlab("Gene Set") + ylab("log(F statistic)") + ggtitle(sprintf("Sun et al. pQTL PC3 FDR<0.05")) + theme(axis.text.y = element_text(angle = 0, vjust = 1, hjust=1,size=8)) +
theme(legend.position="bottom") + guides(color=FALSE)
save_plot(file="~/tmp/gsea_pqtl_pc3.pdf",pp,base_width=10,useDingbats=FALSE)


gsea.react.DT <- newgsea(react)
gsea.hm.DT <- newgsea(hm)
gsea.DT <- rbind(gsea.react.DT,gsea.hm.DT)
#gsea.DT <- gsea.react.DT
## test to see what happens if we just look at one PC
gsea.DT[,p.adj:=p.adjust(p,method='fdr')]


filt.DT <- gsea.DT[pc %in% c('PC1','PC3')  & p.adj<0.01 & sign>0,]
filt.DT[,source:='REACTOME']
filt.DT[grep('^HALLMARK_',gs),source:='HALLMARK']
filt.DT[,gs:=gsub("(HALLMARK|REACTOME)\\_","",gs) %>% gsub("\\_"," ",.) %>% str_trunc(., 70, "right")]
filt.DT[,gs:=factor(gs,levels=filt.DT[order(stat,decreasing=FALSE),]$gs)]

pp <- ggplot(filt.DT,aes(y=pc,x=gs,fill=stat)) + geom_tile() + coord_flip() +
xlab("Component") + ylab("log(F statistic)") + ggtitle(sprintf("Sun et al. pQTL PC3 FDR<0.01")) +
theme(axis.text.y = element_text(angle = 0, vjust = 1, hjust=1,size=10)) + scale_fill_gradient("log(F-statistic)",low='yellow',high='red')
save_plot(file="~/tmp/gsea_sun_pc1_3.pdf",pp,base_width=9,useDingbats=FALSE)





## attempt forest plot of genes

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
basis.DT[,category:='zzz_basis']

forest_plot <- function(proj.dat,basis.dat=basis.DT,pc,fdr_thresh=0.05,theme=NA){
  dat <- proj.dat[variable==pc & p.adj<fdr_thresh,]
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
  dat[,trait:=factor(trait,levels=dat[order(category,delta,decreasing=TRUE),]$trait)]
  #ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point(aes(size=log10(n1))) + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  #coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc)
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }
  ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(sprintf("%s Sun pQTL",pc)) + theme +
  xlab("Trait") + ylab("Change in basis loading from control")
}

control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
Mf <- merge(Mf,control.DT,by.x='variable',by.y='PC')
Mf[,delta:=value-control.loading]
Mf[,variance:=var(delta),by='variable']
Mf[,category:='sun_pqtl']

pdf("~/tmp/forest_plot.pdf",paper="a4r")
lapply(paste('PC',1:10,sep=""),function(pc){
    forest_plot(Mf,pc=pc)
})
dev.off()


pl <- list()
pl[[1]] <- gsea(geneid='entrezgene',hm,test.type='f.test',fdr.method="fdr",fdr.thresh=0.05) %>% phwrap(.,main="F Test Hallmark")
#pl[[2]] <- gsea(lom,test.type='f.test',fdr.method="fdr",fdr.thresh=0.01) %>% phwrap(.,main="F Test Location")
pl[[2]] <- gsea(hm,test.type='bartlett.test',fdr.method="fdr") %>% phwrap(.,main="Bartlett Test Hallmark")
#pl[[4]] <- gsea(lom,test.type='bartlett.test',fdr.method="fdr",fdr.thresh=0.01) %>% phwrap(.,main="F Test Location")
#pl[[2]] <- gsea(hm,test.type='t.test',fdr.method="fdr") %>% phwrap(.,main="T Test Hallmark")
g<-grid.arrange(grobs=pl,ncol=2)
dev.print(pdf,"~/tmp/hallmark_f_eqtlgen.pdf")

## add in reactome just in case

f <- scan("~/tmp/c2.cp.reactome.v6.2.entrez.gmt","character()",sep="\n") %>% strsplit(.,"\t")
react <- lapply(f,function(x) {
  genes <- x[c(-1,-2)] %>% as.integer
  genes <- genes[genes %in% universe]
})
names(react) <- sapply(f,'[[',1)
## remove genesets with no members !
react <- react[sapply(react,length) > 5]

pl <- list()
pl[[1]] <- gsea(react,test.type='f.test',fdr.method="fdr",fdr.thresh=0.05) %>% phwrap(.,main="F Test Reactome")
#pl[[2]] <- gsea(lom,test.type='f.test',fdr.method="fdr",fdr.thresh=0.01) %>% phwrap(.,main="F Test Location")
pl[[2]] <- gsea(react,test.type='bartlett.test',fdr.method="fdr") %>% phwrap(.,main="Bartlett Test Reactome")
#pl[[4]] <- gsea(lom,test.type='bartlett.test',fdr.method="fdr",fdr.thresh=0.01) %>% phwrap(.,main="F Test Location")
#pl[[2]] <- gsea(hm,test.type='t.test',fdr.method="fdr") %>% phwrap(.,main="T Test Hallmark")
g<-grid.arrange(grobs=pl,ncol=2)


## do mean
#gsea(hm,test.type='t.test',fdr.method="fdr") %>% phwrap(.,main="t Test Hallmark")
#dev.print(pdf,"~/tmp/hallmark_t_eqtlgen.pdf")
#library(xtable)
#gsea(lom,test.type='t.test',fdr.method="fdr") %>% xtable


## next perform the same analysis using the t-cell modules from burren et al.

mod.DT <- fread("/home/ob219/share/as_basis/tcell_modules_burren_et_al/13059_2017_1285_MOESM1_ESM")
mod.DT.f <- mod.DT[geneSymbol %in% Mf$hgnc,]

Mf <- Mf[hgnc %in% mod.DT.f$geneSymbol,]
wgcna <- split(mod.DT.f$i,mod.DT.f$all_12)
pl <- list()
pl[[1]] <- gsea(wgcna,geneid='',test.type='f.test',fdr.method="bonferroni") %>% phwrap(.,main="F Test Burren et al")
#pl[[2]] <- gsea(wgcna,geneid='trait',test.type='t.test',fdr.method="bonferroni") %>% phwrap(.,main="t Test Burren et al")
pl[[2]] <- gsea(wgcna,geneid='trait',test.type='bartlett.test',fdr.method="bonferroni") %>% phwrap(.,main="Bartlett Test Burren et al")
g<-grid.arrange(grobs=pl,ncol=2)
dev.print(pdf,"~/tmp/burren_eqtlgen.pdf")


## rather than looking at a module basis use Effron two group model in order to select genes of interest for each PC and then plot genes across all PC's

## this is implemented in locfdr package
library(locfdr)
pc <- Mf[variable=='PC3']
Mf[,fdr:=locfdr(Z,plot=1)$fdr,by='variable']
saveRDS(Mf,"~/tmp/pQTL.RDS")

## below here is not used as yet as really didn't find anything.

## plot a geom tile of Z scores for genes with FDR<0.05
pt <- Mf[fdr<0.01 & variable=='PC3',]
pta <- Mf[trait %in% pt$trait,]
setnames(pta,'variable','pc')
pta <- melt(pta,id.vars=c('external_gene_name','pc'),measure.vars='Z') %>% dcast(.,external_gene_name~pc+variable)
mat <- as.matrix(pta[,-1]) %>% apply(.,2,as.numeric)
rownames(mat) <- pta$external_gene_name
pheatmap(mat)

plot <- Mf[trait %in% pt$trait,]
hc.pc <- dist(mat %>% t) %>% hclust
plot[,external_gene_name:=factor(external_gene_name,levels=hc.gene$labels[hc.gene$order])]
hc.gene <- dist(mat) %>% hclust
plot[,variable:=factor(variable,levels=hc.pc$labels[hc.pc$order] %>% sub("\\_Z","",.))]
plot[]

ggplot(plot,aes(x=external_gene_name,y=variable,fill=Z,label=signif(Z,digits=2))) + geom_tile() + scale_fill_gradientn(colours=c('green','white','red')) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + geom_text() + xlab("Gene") + ylab("PC")

ggplot(plot,aes(x=external_gene_name,y=variable,fill=Z,label=signif(fdr,digits=2))) + geom_tile() + scale_fill_gradientn(colours=c('green','white','red')) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + geom_text() + xlab("Gene") + ylab("PC")

ptpn22 <- fread('/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/sum_stats/eqtlgen/ENSG00000134242.tab')
ptpn22<-ptpn22[p.value<5e-8,]

library(cupcake)

GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'

#man.DT <- fread(SNP_MANIFEST_FILE)

## load data
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
basis.DT <- basis.DT[,.(trait,beta=log(or),p.value,pid,a1,a2)]
basis.DT <- basis.DT[pid %in% ptpn22$pid,]
all.DT <- rbind(ptpn22[,.(trait='PTPN22',beta=or,p.value,pid,a1,a2)],basis.DT)

all.DT <- all.DT[order(pid),]

#M <- melt(all.DT,id.vars=c('pid','trait'),measure.vars='beta')
M <- all.DT
M[,c('chr','pos'):=tstrsplit(pid,':')]

M[,pid:=factor(pid,levels=M[order(chr,pos),]$pid %>% unique)]

ggplot(M,aes(x=pid,y=trait,fill=pmax(p.value,5e-10) %>% -log10(.) * sign(beta),label=signif(exp(beta),digits=2))) + geom_tile() + scale_fill_gradientn("mlp",colours=c('green','white','red')) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + geom_text() + xlab("PID") + ylab("Trait")



## create a heatmap of jia and basis traits to help us interpret

jia.DT <- readRDS("/home/ob219/share/as_basis/GWAS/tmp/jia_bb_summary.RDS")[grep("^jia_",short.trait),.(trait,variable,value,Z)]
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
basis.proj <- readRDS(BASIS_FILE)
bproj <- data.table(trait=rownames(basis.proj$x),basis.proj$x) %>% melt(.,id.vars="trait")
bproj[,Z:=1*sign(value)]
bproj[,basis.trait:=1]
#all <- rbind(bproj,jia.DT)


gi <- rbindlist(list(bproj,jia.DT,M[external_gene_name=='PTPN22',.(trait=external_gene_name,variable,value,Z,basis.trait=0)]),fill=TRUE)
gi[is.na(basis.trait),basis.trait:=0]
gi[,c('trait','variable'):=list(factor(trait,levels=gi[order(basis.trait),]$trait %>% unique),factor(variable,levels=paste('PC',1:11,sep='')))]
gi<-merge(gi[trait=='control',.(variable,tv=value)],gi,by='variable')
gi[,delta:=value-tv]
gi[basis.trait==1,Z:=1*sign(delta)]


ggplot(gi,aes(x=trait,y=variable,fill=Z,label=signif(delta,digits=2))) + geom_tile() + geom_text() + scale_fill_gradientn("Z",colours=c('green','white','red'))




## create a heatmap of all genes

Mf[,raw.p:=pnorm(abs(Z),lower.tail=FALSE) * 2]
Mf[,fdr:=p.adjust(raw.p,method="fdr"),by='variable']
keep.trait <- Mf[fdr<0.05,]$trait
plot.t <- Mf[trait %in% keep.trait,]
ggplot(plot.t,aes(x=external_gene_name,y=variable,fill=Z,label=signif(fdr,digits=1))) + geom_tile() + geom_text() +
scale_fill_gradientn("mlp",colours=c('green','white','blue')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
