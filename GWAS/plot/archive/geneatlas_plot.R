library(cowplot)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/17_07_19_0619_summary_results.RDS'
#RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)

#all.traits <- traits<-split(res.DT$trait,res.DT$category) %>% lapply(.,unique)

## add in basis traits for comparison

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
basis.DT[,category:='basis']

BB_LU <- list(
  CD = 'crohns.disease',
  CEL = 'malabsorption.coeliac.disease',
  MS = 'multiple.sclerosis',
  RA = 'rheumatoid.arthritis',
  SLE = 'systemic.lupus.erythematosis.sle',
  T1D = 'type.1.diabetes',
  UC = 'ulcerative.colitis',
  PBC = 'PBC',
  PSC = 'PSC',
  asthma = 'asthma',
  VIT = 'vitiligo'
)

## remove all srd gene atlas and recalculate p.adj as it is currently wrong
res.DT <- res.DT[category!='geneatlas_srd',]
res.DT[,p.adj:=p.adjust(p.value,method='BH')]

category.foc <- 'geneatlas_icd'
#category.foc <- 'geneatlas_cancer'

talk.DT <- res.DT[category %in% c('bb_disease',category.foc),]
talk.DT<-talk.DT[(category %in% talk.DT[p.adj<0.01,]$category) | category==category.foc,]
at <- talk.DT$category %>% unique
at <- at[!at %in% c(category.foc)]
library(RColorBrewer)
cols<-brewer.pal(length(talk.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c(category.foc,at,'basis')
cols[category.foc] <- 'deeppink2'
cols['basis'] <- '#7c0799'
cols['bb_disease'] <- '#05af6e'
talk.DT[,trait:=gsub("^bb_SRD:","",trait)]


forest_plot_focal_merge <- function(proj.dat,basis.dat=basis.DT,pc,focal,title,cat_levels,fdr_thresh=0.05,theme=NA){
  dat <- proj.dat[variable==pc & (p.adj<fdr_thresh | trait %in% focal | trait %in% unlist(BB_LU)),]
  #dat[trait %in% focal,category:='aa_focal_diseases']
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  ## define a thing called trait label
  dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
  #for(i in seq_along(BB_LU)){
  #  tra <- BB_LU[[i]]
  #  dat[trait==tra,trait:=names(BB_LU)[i]]
  #}
  for(i in seq_along(BB_LU)){
    tra <- names(BB_LU)[i]
    dat[trait==tra,trait:=BB_LU[[i]]]
  }
  dat[,category:=factor(category,levels=names(cat_levels))]
  ## sort focal trait first
  foc.dt <- dat[trait %in% focal,.(trait,delta,category,n1)]
  #foc.dt[,tl:=sprintf("%s(%s)",trait,format(n1,big.mark=",",scientific=FALSE))]
  foc <- foc.dt[order(delta,decreasing=TRUE),.(trait,n1)]
  ## next do the rest
  nf.dt <- dat[!trait %in% focal,.(trait,delta,category,n1)][order(category,decreasing=FALSE),][!duplicated(trait),]
  #nf.dt[,tl:=sprintf("%s(%s)",trait,format(n1,big.mark=",",scientific=FALSE))]
  nf.dt[is.na(n1),tl:=trait]
  nfoc <- nf.dt[order(delta,decreasing=TRUE),.(trait,n1)]
  trait.order <- c(nfoc$trait,foc$trait)
  #dat[!duplicated(trait),][order(trait %in% focal,delta,decreasing=TRUE),]$trait
  #dat[,trait:=factor(trait,levels=dat[order(delta,!trait %in% focal,decreasing=TRUE),]$trait)]
  dat[,trait:=factor(trait,levels=trait.order)]
  #rbind(nfoc,foc)[,tl:=sprintf("%s(%s)",trait,format(n1,big.mark=",",scientific=FALSE))]
  ldat <- rbind(nfoc,foc)[,tl:=sprintf("%s",trait)]
  ldat[is.na(n1),tl:=trait]
  levels(dat$trait) <- ldat$tl
  #ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point(aes(size=log10(n1))) + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  #coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc)
  ## next alter the labels to include sample size
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }
  ggplot(dat,aes(x=trait,y=delta,colour=category,alpha=p.adj<fdr_thresh,lty=p.adj<fdr_thresh)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(title) + theme +
  xlab("Trait") + ylab("Change in basis score from control") + guides(alpha=FALSE,lty=FALSE) +
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.5)) + scale_colour_manual("Category",values=cat_levels,labels=c('Focus',"UKBB\nSRD",'Basis')) +
  scale_linetype_manual(values=c('TRUE'=1,'FALSE'=2)) +
  theme(axis.text.y=element_text(size=12),legend.position="bottom")
}

#only for blood traits where lots of things are significant !
talk.DT<-talk.DT[(category==category.foc & p.adj<0.01) | category!=category.foc,]
talk.DT[,trait:=strtrim(trait, 50)]
all.traits <- traits<-split(talk.DT$trait,talk.DT$category) %>% lapply(.,unique)
#pc<-'PC1'
#pp1 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits[category.foc] %>% unlist,title=pc,cat_levels=cols)

pdf(file="~/tmp/geneatlas_160719_fdr01.pdf",paper="a4r",onefile=TRUE)
lapply(paste('PC',1:11,sep=''),function(pc){
  forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits[category.foc] %>% unlist,title=pc,cat_levels=cols,fdr_thresh=0.01)
})
dev.off()

## create a summary heatmap that is clustered

RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/17_07_19_0619_summary_results.RDS'
#RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)[category=='geneatlas_icd',]
## only use single category icd codes at this stage
library(stringr)
pat <- "[A-Z][0-9]{2}"
res.DT[,icd.count:=str_count(trait,pat)]
leaf.DT <- res.DT[str_count(trait,pat)==1,]
node.DT <- res.DT[str_count(trait,pat)>1,]
fail.DT <- res.DT[str_count(trait,pat)<1,]


#keep <- leaf.DT[p.adj<0.01,]$trait %>% unique %>% as.character()
keep <- leaf.DT[p.adj<0.05,]$trait %>% unique %>% as.character()
leaf.DT <- leaf.DT[trait %in% keep,]
mat.dt <- melt(leaf.DT[,.(pc=variable,trait,delta)],measure.vars='delta') %>% dcast(.,trait~pc)
mat <- as.matrix(mat.dt[,-1])
rownames(mat) <- mat.dt$trait
hc.mat <- dist(mat) %>% hclust
leaf.DT[,trait:=factor(trait,levels=hc.mat$labels[hc.mat$order])]
leaf.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]

leaf.DT[,is.sig:='']
leaf.DT[p.adj<0.05,is.sig:='*']
pp1 <- ggplot(leaf.DT[category==category.foc & variable!='PC12',],aes(x=variable,y=trait,fill=delta,label=is.sig)) + geom_tile() +
geom_text(size=6) +
scale_fill_gradient2("Difference\nfrom control") + xlab("Principal Component") +
ylab("Disease/Subtype") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


keep <- node.DT[p.adj<0.05,]$trait %>% unique %>% as.character()
node.DT <- node.DT[trait %in% keep,]
mat.dt <- melt(node.DT[,.(pc=variable,trait,delta)],measure.vars='delta') %>% dcast(.,trait~pc)
mat <- as.matrix(mat.dt[,-1])
rownames(mat) <- mat.dt$trait
hc.mat <- dist(mat) %>% hclust
node.DT[,trait:=factor(trait,levels=hc.mat$labels[hc.mat$order])]
node.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]

node.DT[,is.sig:='']
node.DT[p.adj<0.05,is.sig:='*']
pp1 <- ggplot(node.DT[category==category.foc & variable!='PC12',],aes(x=variable,y=trait,fill=delta,label=is.sig)) + geom_tile() +
geom_text(size=6) +
scale_fill_gradient2("Difference\nfrom control") + xlab("Principal Component") +
ylab("Disease/Subtype") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

### comparing neale et al self reported disease with

library(cowplot)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/17_07_19_0619_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)

ga.srd <- res.DT[category=='geneatlas_srd',.(trait=gsub("GA:","",trait),pc=variable,ga.Z=Z,ga.cases=n1,ga.controls=n0,ga.p.adj=p.adj)]
neale.srd <- res.DT[category=='bb_disease',.(trait=gsub("bb_SRD:","",trait),pc=variable,neale.Z=Z,neale.cases=n1,neale.controls=n0,neale.p.adj=p.adj)]

M <- merge(ga.srd,neale.srd,by=c('trait','pc'),all.y=TRUE)

M[,pc:=factor(pc,levels=paste0('PC',1:12))]
M[,label:='NONE']
#M[abs(neale.Z)>2.8 & abs(ga.Z)<2.8,label:='NEALE']
#M[abs(neale.Z)<2.8 & abs(ga.Z)>2.8,label:='GA']
#M[abs(neale.Z)>2.8 & abs(ga.Z)>2.8,label:='BOTH']
thresh <- 0.01
M[neale.p.adj<thresh & ga.p.adj>thresh, label:='NEALE']
M[neale.p.adj>thresh & ga.p.adj<thresh, label:='GA']
M[neale.p.adj<thresh & ga.p.adj<thresh, label:='BOTH']

## compare cases

pp1 <- ggplot(M,aes(x=neale.cases,y=ga.cases)) + geom_point() +
geom_abline(col='red',lty=2) + xlab("Neale Cases") + ylab("GeneAtlas Cases")
pp2 <- ggplot(M,aes(x=neale.controls,y=ga.controls)) + geom_point() +
geom_abline(col='red',lty=2) + xlab("Neale Controls") + ylab("GeneAtlas Controls")

plot_grid(pp1,pp2)

## what self reported phenotypes differ between the two studies

Mall <- merge(ga.srd,neale.srd,by=c('trait','pc'),all=TRUE)[trait!='unclassifiable']
Mall[is.na(ga.cases),.(trait,neale.cases)] %>% unique
Mall[is.na(neale.cases),.(trait,ga.cases)] %>% unique

library(ggrepel)
M[,dl:='']
M[label!='NONE',dl:=trait]

ggplot(M,aes(x=neale.Z,y=ga.Z,col=label,label=dl)) +
geom_point() + geom_abline(a=0,b=1,col='black',lty=2,alpha=0.3) +
theme_bw() + facet_wrap(~pc,scales='free') + geom_text_repel()
