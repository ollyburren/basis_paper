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

args <- list(
    snp_support_file='/rds/user/ob219/hpc-work/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab',
    outfile="~/tmp/testbasis.RDS",
    manifest_file='/home/ob219/git/as_basis/manifest/as_manifest_feb_2018_w_ms.csv',
    gwas_data_dir='/rds/user/ob219/hpc-work/as_basis/gwas_stats/filter_feb_2018_w_ms/aligned'
)

## loaf in basis data
basis.DT<-get_gwas_data(args$manifest_file,args$snp_support_file,args$gwas_data_dir)
#SLE data is a bit ropey and contains p-vals that are zero - not sure why but can
#compute using beta and standard error. I checked and makes little different to PCA's
basis.DT[p.value==0,]$p.value<-1e-13
shrink.DT<-compute_shrinkage_metrics(basis.DT)
basis.DT[,c('chr','position'):=tstrsplit(pid,':')]
shrink.DT[,c('chr','position'):=tstrsplit(pid,':')]
basis.DT[,position:=as.numeric(position)]

## for slide one prepare a forest plot for PTPN22

 ptpn22.dat <- basis.DT[pid=='1:114377568',]

 ## compute beta and standard error of beta

 ptpn22.dat[,c('beta','se.beta'):=list(-log(or),abs(log(or))/qnorm(p.value/2,lower.tail=FALSE))]
 ptpn22.dat[,trait.order:=factor(trait,levels=ptpn22.dat[order(beta),]$trait)]


 ## create a forest plot of these


 library(cowplot)

 ppa <- ggplot(ptpn22.dat,aes(x=trait.order,y=beta)) + geom_point() +
 geom_errorbar(aes(ymin=beta-(1.96 * se.beta),ymax=beta+(1.96 * se.beta))) +
 coord_flip() + geom_hline(yintercept=0,color='firebrick1',lty=2) + xlab("Immune-mediated disease") +
 ylab("log(Odds Ratio)") + background_grid(major = "xy", minor = "none") + ggtitle("PTPN22 R620W (rs2476601)")
 save_plot(ppa,file="~/tmp/ptpn22_forest_plot.pdf")

 ## replot higlighting GW sig results

 ppb <- ggplot(ptpn22.dat,aes(x=trait.order,y=beta,color=p.value<5e-8)) + geom_point() +
 geom_errorbar(aes(ymin=beta-(1.96 * se.beta),ymax=beta+(1.96 * se.beta))) +
 coord_flip() + geom_hline(yintercept=0,color='firebrick1',lty=2) + xlab("Immune-mediated disease") +
 ylab("log(Odds Ratio)") + background_grid(major = "xy", minor = "none") + ggtitle("PTPN22 R620W (rs2476601)") +
 scale_colour_manual(values=c('TRUE'='red','FALSE'='lightgrey'),guide=FALSE)
 save_plot(ppb,file="~/tmp/ptpn22_forest_plot_highlight.pdf")

 ## do a tile plot for this data

 ## get a region surrounding our SNP of interest
offset <- 100000
ptpn22.reg <- basis.DT[chr==1 & between(position,114377568-offset,114377568+offset)]
ptpn22.reg[,c('beta','se.beta'):=list(-log(or),abs(log(or))/qnorm(p.value/2,lower.tail=FALSE))]

pos <- ptpn22.reg$position %>% unique
pos <- pos[order(pos)]
idx <- which(pos==114377568)
## take three either side
pos <- pos[(idx-3):(idx+3)]

ptpn22.reg <- ptpn22.reg[position %in% pos,]
## create a naming scheme for this
nas <- ptpn22.reg[,list(pos=position),by=position]
nas[,snp.name:=paste('SNP',-3:+3,sep='')]
nas[snp.name=='SNP0',snp.name:='rs2476601']
nas[snp.name=='SNP1',snp.name:='SNP+1']
nas[snp.name=='SNP2',snp.name:='SNP+2']
nas[snp.name=='SNP3',snp.name:='SNP+3']
ptpn22.reg <- merge(ptpn22.reg,nas,by='position')

ptpn22.reg[,snp.name:=factor(snp.name,levels=nas$snp.name)]
ptpn22.reg[,trait.order:=factor(trait,levels=levels(ptpn22.dat$trait.order))]

offset <- 100000
ptpn22.reg <- basis.DT[chr==1 & between(position,114377568-offset,114377568+offset)]
ptpn22.reg[,c('beta','se.beta'):=list(-log(or),abs(log(or))/qnorm(p.value/2,lower.tail=FALSE))]

pos <- ptpn22.reg$position %>% unique
pos <- pos[order(pos)]
idx <- which(pos==114377568)
## take three either side
pos <- pos[(idx-7):(idx+7)]

ptpn22.reg <- ptpn22.reg[position %in% pos,]

M <- melt(ptpn22.reg,id.vars=c('trait','position'),measure.vars='beta')
M <- dcast(M,trait~position+variable)
traits<-M$trait
M<-M[,-1] %>% as.matrix
rownames(M) <- traits



ppc <- ggplot(ptpn22.reg,aes(y=trait.order,x=snp.name,fill=beta,label=signif(beta,digits=1))) + geom_tile() +
geom_text() +
scale_fill_gradient2(guide=FALSE,low='red',high='blue') + xlab("SNPs") + ylab("Immune-mediated traits") +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
save_plot(ppc,file="~/tmp/region_matrix.pdf")

## lets do prcomp on this tiny region just to illustrate the matrices that we get back



## what about if we heatmap the whole BiomartGeneRegionTrack


ptpn22.reg <- basis.DT[ld.block==140,]
ptpn22.reg[,c('beta','se.beta'):=list(-log(or),abs(log(or))/qnorm(p.value/2,lower.tail=FALSE))]
nas <- ptpn22.reg[,list(pos=position),by=position]
nas[,snp.name:=paste('SNP',1:.N,sep='')]
nas[pos==114377568,snp.name:='rs2476601']
ptpn22.reg <- merge(ptpn22.reg,nas,by='position')

ptpn22.reg[,snp.name:=factor(snp.name,levels=nas$snp.name)]
ptpn22.reg[,trait.order:=factor(trait,levels=levels(ptpn22.dat$trait.order))]

ppd <- ggplot(ptpn22.reg,aes(y=trait.order,x=snp.name,fill=beta,label=signif(beta,digits=2))) + geom_tile() +
scale_fill_gradient2(guide=FALSE,low='red',high='blue') + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
xlab("SNPs") + ylab("Immune-mediated traits")
save_plot(ppd,file="~/tmp/region_matrix_bigger.pdf")

M <- melt(ptpn22.reg,id.vars=c('trait','position'),measure.vars='beta')
M <- dcast(M,trait~position+variable)
traits<-M$trait
M<-M[,-1] %>% as.matrix
## add in the control vector
rbind(M,control=rep(0,ncol(M)))
rownames(M) <- traits
M <- rbind(M,control=rep(0,ncol(M)))

pc <- prcomp(M)

## draw a picture of pc$x -- the principal component scores

pc.score.dt <- data.table(trait=rownames(pc$x),pc$x)

p1.var <- signif(summary(pc)[['importance']][2,]['PC1']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC1',.)
p2.var <- signif(summary(pc)[['importance']][2,]['PC2']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC2',.)

ptitle <- sprintf("chr1:%s-%s",ptpn22.reg$position %>% min %>% format(.,big.mark=",",scientific=FALSE),ptpn22.reg$position %>% max %>% format(.,big.mark=",",scientific=FALSE))

library(ggrepel)
pcaa <- ggplot(pc.score.dt,aes(x=PC1,y=PC2,label=trait)) + geom_point(size=2,color='black') + geom_text_repel() +
geom_hline(yintercept=pc.score.dt[trait=='control']$PC2,lty=2) + geom_vline(xintercept=pc.score.dt[trait=='control']$PC1,lty=2) +
xlab(p1.var) + ylab(p2.var) + guides(color=FALSE) + ggtitle(ptitle)
save_plot(pcaa,file="~/tmp/PCA_noproj.pdf")



## project on other data

bb.DT <- get_gwas_data(args$manifest_file,args$snp_support_file,args$gwas_data_dir,trait_list=c('bb_UC','bb_CD','bb_SLE','bb_MS','bb_RA','bb_CEL','bb_T1D','bb_asthma'))
ptpn22.reg <- bb.DT[ld.block==140,]
ptpn22.reg[,c('beta','se.beta'):=list(-log(or),abs(log(or))/qnorm(p.value/2,lower.tail=FALSE))]
ptpn22.reg[,c('chr','position'):=tstrsplit(pid,':')]
nas <- ptpn22.reg[,list(pos=position),by=position]
nas[,snp.name:=paste('SNP',1:.N,sep='')]
nas[pos==114377568,snp.name:='rs2476601']
ptpn22.reg <- merge(ptpn22.reg,nas,by='position')
M_BB <- melt(ptpn22.reg,id.vars=c('trait','position'),measure.vars='beta')
M_BB <- dcast(M_BB,trait~position+variable)
traits<-M_BB$trait
M_BB<-M_BB[,-1] %>% as.matrix
rownames(M_BB) <- traits

bb_score.dt <- predict(pc,M_BB)
bb_score.dt <- data.table(trait=rownames(bb_score.dt),bb_score.dt)

pc.score.dt[,projected:=FALSE]
bb_score.dt[,projected:=TRUE]


pcab <- ggplot(rbind(pc.score.dt,bb_score.dt),aes(x=PC1,y=PC2,label=trait,color=projected)) + geom_point(size=2) + geom_text_repel() +
geom_hline(yintercept=pc.score.dt[trait=='control']$PC2,lty=2) + geom_vline(xintercept=pc.score.dt[trait=='control']$PC1,lty=2) +
xlab(p1.var) + ylab(p2.var)  + scale_colour_manual(values=c('TRUE'='firebrick1','FALSE'='black'),guide=FALSE) + ggtitle(ptitle)
save_plot(pcab,file="~/tmp/PCA_proj.pdf")

## plot rotations

ev <- data.table(snp=rownames(pc$rotation),pc$rotation)
rot.dat <- melt(ev,id.vars='snp')

rot.dat[,variable:=factor(variable,levels=paste0('PC',11:1))]

basis <- ggplot(rot.dat,aes(x=snp,y=variable,fill=value)) + geom_tile() + scale_fill_gradient2(guide=FALSE,low='red',high='blue') +
theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
xlab("SNPs") + ylab("Principal components")
save_plot(basis,file="~/tmp/basis.pdf")


dat <- rbind(pc.score.dt,bb_score.dt)
dat.mat <- as.matrix(dat[,2:12,with=FALSE])
rownames(dat.mat) <- dat$trait
pdf("~/tmp/hclust_beta.pdf")
dist(dat.mat[,1:2]) %>% hclust %>% plot
dev.off()
## what about with shrinkages ?ws_emp_shrinkage

basis.mat.emp <- create_ds_matrix(basis.DT[ld.block==140,],shrink.DT[ld.block==140,],'ws_emp_shrinkage')
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)

pc.emp.score.dt <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
p1.var <- signif(summary(pc.emp)[['importance']][2,]['PC1']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC1',.)
p2.var <- signif(summary(pc.emp)[['importance']][2,]['PC2']*100,digits=3) %>% sprintf("%s (%.1f%%)",'PC2',.)
ptitle <- sprintf("chr1:%s-%s",ptpn22.reg$position %>% min %>% format(.,big.mark=",",scientific=FALSE),ptpn22.reg$position %>% max %>% format(.,big.mark=",",scientific=FALSE))

library(ggrepel)
pcaa <- ggplot(pc.emp.score.dt,aes(x=PC1,y=PC2,label=trait)) + geom_point(size=2,color='black') + geom_text_repel() +
geom_hline(yintercept=pc.emp.score.dt[trait=='control']$PC2,lty=2) + geom_vline(xintercept=pc.emp.score.dt[trait=='control']$PC1,lty=2) +
xlab(p1.var) + ylab(p2.var) + guides(color=FALSE) + ggtitle(ptitle)

ev <- data.table(snp=rownames(pc.emp$rotation),pc.emp$rotation)
rot.dat <- melt(ev,id.vars='snp')

rot.dat[,variable:=factor(variable,levels=paste0('PC',11:1))]

basis <- ggplot(rot.dat,aes(x=snp,y=variable,fill=value)) + geom_tile() + scale_fill_gradient2(guide=FALSE,low='red',high='blue') +
theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
xlab("SNPs") + ylab("Principal components")


save_plot(pcaa,file="~/tmp/PCA_noproj.pdf")





e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

hgnc_name<-'CTLA4'

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
raw.dat <- raw.dat[!trait %in% c('asthma','CD','UC'),]
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


bshrink.gr <- with(shrink.DT[pid %in% beta$pid,],GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(position),width=1L),shrink=bshrink))
wshrink.gr <- with(shrink.DT[pid %in% beta$pid,],GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(position),width=1L),wshrink=ws_ppi))

## also create an object that is the shrunk OR for the dataset

#compute Wakefields
ppi<-raw.dat[,list(pid=pid,ppi=wakefield_pp(p.value,maf,unique(n),unique(n1/n))),by=c('trait','ld.block')][,.(trait,ppi,pid,ld.block)]
setkeyv(ppi,c('pid','trait'))
setkeyv(raw.dat,c('pid','trait'))
raw.dat <- ppi[raw.dat]
melt.DT<-melt(raw.dat,id.vars=c('pid','chr','position','trait'),measure.vars=c('beta','mlp','ppi'))
ppid<-dcast(melt.DT[variable=="ppi",],pid+chr+position~trait+variable)
ppi.gr <- with(ppid,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L)))
mcols(ppi.gr)<- as.data.frame(ppid[,4:ncol(ppid)]) %>% DataFrame


tracks<-list()
tracks$itr <- IdeogramTrack(genome="hg19", chromosome=seqlevels(comb.gr),bands=ideo)
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
tracks$beta<- DataTrack(name="log(OR)",cex=1,beta.gr,groups=gsub("_beta","",names(mcols(beta.gr))),legend=FALSE)
#displayPars(tracks$combmlp)$background.title <- 'black'
displayPars(tracks$beta)$background.title <- 'black'
tracks$ppi <- DataTrack(name="sCVPP",cex=1,ppi.gr,groups=gsub("_ppi","",names(mcols(ppi.gr))),legend=FALSE)
#displayPars(tracks$ppi)$cex <- 1.5
displayPars(tracks$ppi)$jitter.x <- FALSE
displayPars(tracks$ppi)$background.title <- 'black'
#tracks$shrink <- DataTrack(name="PPA Shrink",bshrink.gr,col='dodgerblue')
#displayPars(tracks$shrink)$background.title <- 'dodgerblue'
tracks$shrink2 <- DataTrack(name="X-disease sCVPP",wshrink.gr,col='black',cex=1)
displayPars(tracks$shrink2)$background.title <- 'black'
tracks$sbeta <- DataTrack(name="Shrunk log(OR)",cex=1,shrunk_beta.gr,groups=gsub("_shrunk_beta","",names(mcols(beta.gr))),legend=FALSE)
displayPars(tracks$sbeta)$background.title <- 'firebrick'
#plotTracks(tracks,from=min(raw.dat[ld.block %in% ldb,]$position),to=max(raw.dat[ld.block %in% ldb,]$position))

## zoom into ld block with highest posterior
tmp.DT<-shrink.DT[ld.block %in% ldb,list(tot=sum(ws_ppi)),by=ld.block]
max.ldb <- tmp.DT[which.max(tot),]$ld.block
#plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta')],from=204732511-400000,to=204738683+400000)
#plotTracks(tracks,from=min(raw.dat[ld.block %in% max.ldb,]$position),to=max(raw.dat[ld.block %in% max.ldb,]$position))

pdf("~/tmp/ctla4_shrinkage_1.pdf",useDingbats=FALSE)
#plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta','ppi','shrink2')],from=min(raw.dat[ld.block %in% max.ldb,]$position),to=max(raw.dat[ld.block %in% max.ldb,]$position))
plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta')],from=204571198-20000,to=204826298+20000)
dev.off()

pdf("~/tmp/ctla4_shrinkage_2.pdf",useDingbats=FALSE)
#plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta','ppi','shrink2')],from=min(raw.dat[ld.block %in% max.ldb,]$position),to=max(raw.dat[ld.block %in% max.ldb,]$position))
plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta','ppi')],from=204571198-20000,to=204826298+20000)
dev.off()

pdf("~/tmp/ctla4_shrinkage_3.pdf",useDingbats=FALSE)
#plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta','ppi','shrink2')],from=min(raw.dat[ld.block %in% max.ldb,]$position),to=max(raw.dat[ld.block %in% max.ldb,]$position))
plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta','shrink2')],from=204571198-20000,to=204826298+20000)
dev.off()

pdf("~/tmp/ctla4_shrinkage_4.pdf",useDingbats=FALSE)
#plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta','ppi','shrink2')],from=min(raw.dat[ld.block %in% max.ldb,]$position),to=max(raw.dat[ld.block %in% max.ldb,]$position))
plotTracks(tracks[names(tracks) %in% c('gene','mlp','beta','sbeta')],from=204571198-20000,to=204826298+20000)
dev.off()
