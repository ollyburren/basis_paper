BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
all.DT <- readRDS("/home/ob219/share/as_basis/GWAS/RESULTS/08_04_19_summary_results.RDS")

## FOREST PLOT CODE

library(cowplot)

res.DT <- all.DT
all.traits <- traits<-split(res.DT$trait,res.DT$category) %>% lapply(.,unique)
## add in basis traits for comparison read in again just in case changed
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
basis.DT[,category:='basis']

## another idea is to collapse biobank and basis traits where applicable

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
  PSA = 'PSA'
)

talk.DT <- res.DT
talk.DT<-talk.DT[category %in% talk.DT[p.adj<0.05,]$category,]
at <- talk.DT$category %>% unique
at <- at[!at %in% c('bowes_jia_2019')]
library(RColorBrewer)
cols<-brewer.pal(length(talk.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c('bowes_jia_2019',at,'basis')
cols['bowes_jia_2019'] <- 'deeppink2'
cols['basis'] <- '#7c0799'
cols['bb_disease'] <- '#05af6e'
talk.DT[,trait:=gsub("^bb_","",trait)]
talk.DT[,trait:=gsub("^jia_","",trait)]

## perhaps the other way around

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
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.5)) + scale_colour_manual("Category",values=cat_levels,labels=c('JIA','UKBB SR Disease','Basis')) +
  scale_linetype_manual(values=c('TRUE'=1,'FALSE'=2)) +
  theme(axis.text.y=element_text(size=12))

}

pc <- 'PC3'
title <- sprintf("JIA subtypes %s (Variance Explained %0.1f%%) OLD Basis",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
pp1 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits['bowes_jia_2019'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols)
