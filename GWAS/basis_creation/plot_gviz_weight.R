## SIMULATION WORKBOOK
#Misc code and notes for t1d_bootstrap simulation

## check the current imputed objects and get a list of available SNPs
## MUST BE DONE ON CSD3

# chromosome 1 failed ;(

# In meantime take a look weighting

library(devtools)
install_github('ollyburren/cupcake')
library(cupcake)
library(Gviz)
library(biomaRt)
library(cowplot)


SHRINKAGE_METHOD<-'ws_emp_shrinkage'
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
GWAS_DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_0619.RDS'



## load shrinkage data

shrink.DT <- readRDS(SHRINKAGE_FILE)

## SNPID is 12:6440009 (rs1800693)

#shrink.DT[pid=='12:6440009',]
#shrink.DT[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
#shrink.DT[chr==12 & between(pos,6440009-1e4,6440009+1e4),]$ld.block

#shrink.DT[pid=='1:114377568',]$ld.block



## what about if we heatmap the whole BiomartGeneRegionTrack

basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
#block.reg <- basis.DT[ld.block==2211,]
#block.reg[,c('beta','se.beta'):=list(-log(or),abs(log(or))/qnorm(p.value/2,lower.tail=FALSE))]
basis.DT[,c('chr','position'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]

## project on other data


e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

hgnc_name<-'TNFRSF1A'

target<-BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        external_name = hgnc_name,
        collapseTranscripts="gene",shape = "arrow",
        stackHeight=0.1, filters=list(with_ox_refseq_mrna=T,hgnc_symbol=hgnc_name))

rchr <- gsub("^chr","",unique(seqnames(target)))
rstart<-min(start(range(target)))
rend<-max(end(range(target)))

ideo<-fread("/rds/user/ob219/hpc-work/DATA/ucsc/hg19_cytoBandIdeo.txt")
setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))

ldb<-unique(subset(basis.DT,chr==rchr & between(position,rstart-1e4,rend+1e4))$ld)
raw.dat<-subset(basis.DT,ld.block %in% c(ldb-1,ldb,ldb+1))
raw.dat[,slp:=sign(log(or)) * (-log10(p.value) %>% pmin(.,-log10(5e-15)))]
raw.dat[,beta:=log(or)]
raw.dat[,mlp:=-log10(p.value)]
## for clarity remove some of the traits
raw.dat <- raw.dat[trait %in% c('MS','PBC','T1D'),]
raw.dat <- merge(raw.dat,shrink.DT[,.(pid,ws_emp_shrinkage)],by='pid')
raw.dat[,shrunk_beta:=beta*ws_emp_shrinkage]

melt.DT<-melt(raw.dat,id.vars=c('pid','chr','position','trait'),measure.vars=c('beta','mlp','shrunk_beta'))
beta<-dcast(melt.DT[variable=="beta",],pid+chr+position~trait+variable)
beta.gr <- with(beta,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L)))
mcols(beta.gr)<- as.data.frame(beta[,4:ncol(beta)]) %>% DataFrame

mlp<-dcast(melt.DT[variable=="mlp",],pid+chr+position~trait+variable)
mlp.gr <- with(mlp,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L)))
mcols(mlp.gr)<- as.data.frame(mlp[,4:ncol(mlp)]) %>% DataFrame

shrunk_beta<-dcast(melt.DT[variable=="shrunk_beta",],pid+chr+position~trait+variable)
shrunk_beta.gr <- with(shrunk_beta,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L)))
mcols(shrunk_beta.gr)<- as.data.frame(shrunk_beta[,4:ncol(shrunk_beta)]) %>% DataFrame


#bshrink.gr <- with(shrink.DT[pid %in% beta$pid,],GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(position),width=1L),shrink=ws_emp_shrinkage))
wshrink.gr <- with(shrink.DT[pid %in% beta$pid,],GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(pos),width=1L),wshrink=ws_emp_shrinkage))

## also create an object that is the shrunk OR for the dataset

library('latex2exp')

tracks<-list()
tracks$itr <- IdeogramTrack(genome="hg19", chromosome=seqlevels(wshrink.gr),bands=ideo)
tracks$axis<-GenomeAxisTrack()
tracks$gene <- BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        chromosome = seqlevels(beta.gr), start = min(mlp$position), end = max(mlp$position),
        collapseTranscripts="meta",shape = "arrow", transcriptAnnotation = "symbol",fill='black',col='black',
        stackHeight=0.5, filters=list(with_ox_refseq_mrna=T),cex=5,background.title='black',cex.group=1,protein_coding='black',
        fontcolor.group='black',just.group = "above")
tracks$mlp<- DataTrack(name="-log10(p)",cex=1,mlp.gr,groups=gsub("_mlp","",names(mcols(mlp.gr))),cex.legend=1)
displayPars(tracks$mlp)$background.title <- 'black'
tracks$beta<- DataTrack(name="beta",cex=1,beta.gr,groups=gsub("_beta","",names(mcols(beta.gr))),legend=FALSE)
#displayPars(tracks$combmlp)$background.title <- 'black'
#displayPars(tracks$beta)$background.title <- 'black'
#tracks$ppi <- DataTrack(name="sCVPP",cex=1,ppi.gr,groups=gsub("_ppi","",names(mcols(ppi.gr))),cex.legend=1,legend=FALSE)
#displayPars(tracks$ppi)$cex <- 1.5
#displayPars(tracks$ppi)$jitter.x <- FALSE
#displayPars(tracks$ppi)$background.title <- 'black'
#tracks$shrink <- DataTrack(name="PPA Shrink",bshrink.gr,col='dodgerblue')
#displayPars(tracks$shrink)$background.title <- 'dodgerblue'
tracks$shrink2 <- DataTrack(name="Weight",wshrink.gr,col='black',cex=1)
displayPars(tracks$shrink2)$background.title <- 'black'
tracks$sbeta <- DataTrack(name="Weighted beta",cex=1,shrunk_beta.gr,groups=gsub("_shrunk_beta","",names(mcols(beta.gr))),legend=FALSE)
displayPars(tracks$sbeta)$background.title <- 'firebrick'
#plotTracks(tracks,from=min(raw.dat[ld.block %in% ldb,]$position),to=max(raw.dat[ld.block %in% ldb,]$position))

## zoom into ld block with highest posterior
tmp.DT<-shrink.DT[ld.block %in% ldb,list(tot=sum(ws_ppi)),by=ld.block]
max.ldb <- tmp.DT[which.max(tot),]$ld.block
#plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta')],from=204732511-400000,to=204738683+400000)
#plotTracks(tracks,from=min(raw.dat[ld.block %in% max.ldb,]$position),to=max(raw.dat[ld.block %in% max.ldb,]$position))


pdf("~/tmp/tnfrsf1a_shrinkage_1.pdf",useDingbats=FALSE)
#plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta','ppi','shrink2')],from=min(raw.dat[ld.block %in% max.ldb,]$position),to=max(raw.dat[ld.block %in% max.ldb,]$position))
plotTracks(tracks,from=(shrink.DT[ld.block==2211,]$pos %>% min)-1e4,to=(shrink.DT[ld.block==2211,]$pos %>% max)+1e4)
dev.off()
