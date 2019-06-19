## FOREST PLOTS FOR particular diseases
library(cowplot)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/19_06_19_0619_summary_results.RDS'
#RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)

all.traits <- traits<-split(res.DT$trait,res.DT$category) %>% lapply(.,unique)

## add in basis traits for comparison

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_vit_t2d.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
basis.DT[,category:='zzz_basis']




## work out which PC's we want to draw forest plots for
#PCs <- res.DT[trait %in% FOCUS.DISEASES & variable!='PC11' & p.adj<SIG.THRESH,]$variable %>% unique

#pc <- 'PC3'

forest_plot_focal <- function(proj.dat,basis.dat=basis.DT,pc,focal,title,fdr_thresh=0.05,theme=NA){
  dat <- proj.dat[variable==pc & (p.adj<fdr_thresh | trait %in% focal),]
  dat[trait %in% focal,category:='aa_focal_diseases']
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
  dat[,trait:=factor(trait,levels=dat[order(category,delta,decreasing=TRUE),]$trait)]
  #ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point(aes(size=log10(n1))) + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  #coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc)
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }
  ggplot(dat,aes(x=trait,y=delta,colour=category,alpha=p.adj<fdr_thresh)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(title) + theme +
  xlab("Trait") + ylab("Change in basis loading from control") + guides(alpha=FALSE) + scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.3))
}

## for talk leave out some of the stuff to make figure clearer

talk.DT <- res.DT[category %in% c('bb_medications','astle_blood','bowes_jia','astle_blood','bb_disease')]



forest_plot_focal(talk.DT,pc='PC1',focal=all.traits[['bowes_jia']],title="JIA subtypes PC1")
dev.print(pdf,file="~/tmp/pc1_plot.pdf",useDingbats=FALSE)
forest_plot_focal(talk.DT,pc='PC3',focal=all.traits[['bowes_jia']],title="JIA subtypes PC3")
dev.print(pdf,file="~/tmp/pc3_plot.pdf",useDingbats=FALSE)


forest_plot_focal(proj.dat=res.DT,pc='PC1',focal=all.traits[['myogen']],title="Myogen Myositis PC1")
forest_plot_focal(res.DT,pc='PC10',focal=all.traits[['myogen']],title="Myogen Myositis PC10",theme=theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3),axis.text.y=element_text(size=9)))
pdf(file="~/tmp/for_wendy_19_12_2018.pdf",paper="a4r",width=14,height=8)
forest_plot_focal(res.DT,pc='PC1',focal=all.traits[['bowes_jia']],title="JIA subtypes PC3")
dev.off()
forest_plot_focal(res.DT,pc='PC4',focal=all.traits[['lyons_vasculitis']],title="Lyons et al. Vasculitis PC4")
forest_plot_focal(res.DT,pc='PC6',focal=all.traits[['lyons_vasculitis']],title="Lyons et al. Vasculitis PC6")


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
  ggplot(dat,aes(x=trait,y=delta,colour=category,)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc) + theme +
  xlab("Trait") + ylab("Change in basis loading from control")
}

## extra slide for talk showing correlation with medications
talk.DT <- res.DT[category %in% c('bb_medications','bb_disease')]
forest_plot_med <- function(proj.dat,basis.dat=basis.DT,pc,fdr_thresh=0.05,theme=NA){
  dat <- proj.dat[variable==pc & (p.adj<fdr_thresh | category=='basis'),]
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  #dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
  dat[,trait:=factor(trait,levels=dat[order(category,delta,decreasing=TRUE),]$trait)]
  #ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point(aes(size=log10(n1))) + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  #coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc)
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }
  ggplot(dat,aes(x=trait,y=delta,colour=category,)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc) + theme +
  xlab("Trait") + ylab("Change in basis loading from control")
}

pdf("~/tmp/medication_example.pdf",onefile=TRUE,paper="a4r",width=14,height=8)
forest_plot_med(talk.DT,pc='PC1')
dev.off()

talk.DT <- res.DT[category %in% c('bb_medications','bb_disease','bb_cancer')]
rbind(talk.DT,basis.DT,fill=TRUE)
basis.DT[,category:='basis']
forest_plot_med(rbind(talk.DT,basis.DT,fill=TRUE),pc='PC5',fdr=0.05)

everything.DT <- res.DT

pdf(file="~/tmp/everything_forest_0619.pdf",width=14,height=20,onefile=TRUE)
lapply(paste('PC',1:13,sep=''),function(pc){
  if(pc != 'PC122'){
    forest_plot(everything.DT,pc=pc)
  }else{
    forest_plot(everything.DT,pc=pc,theme=theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3),axis.text.y=element_text(size=6)))
  }
})
dev.off()

lucy.DT <- res.DT[!category %in% c('lee_CD_prognosis','lyons_vasculitis','bb_medications')]
pdf(file="~/tmp/lucy_forest_PC%d.pdf",paper="a4r",width=14,height=8,onefile=FALSE)
lapply(paste('PC',1:10,sep=''),function(pc){
  if(pc != 'PC10'){
    forest_plot(lucy.DT,pc=pc)
  }else{
    forest_plot(lucy.DT,pc=pc,theme=theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3),axis.text.y=element_text(size=6)))
  }

})
dev.off()

## without medications for clarity

pdf(file="~/tmp/lucy_forest__no_med_PC%d.pdf",paper="a4r",width=14,height=8,onefile=FALSE)
lapply(paste('PC',1:10,sep=''),function(pc){
    forest_plot(lucy.DT[category!='bb_medications',],pc=pc)
})
dev.off()


### for talk clean up and only present a few traits

talk.DT <- res.DT[!category %in% c('lee_CD_prognosis','lyons_vasculitis','bb_medications','estrada_NMO','myogen')]
forest_plot(lucy.DT,pc='PC1')


pdf(file="~/tmp/basis_results.pdf",paper="a4r",width=14,height=8)

  if(pc != 'PC10'){
    forest_plot(res.DT,pc=pc)
  }else{
    forest_plot(res.DT,pc=pc,theme=theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3),axis.text.y=element_text(size=6)))
  }
})
dev.off()





forest_plot_focal <- function(proj.dat,basis.dat=basis.DT,pc,focal,title,cat_levels,fdr_thresh=0.05,theme=NA){
  dat <- proj.dat[variable==pc & (p.adj<fdr_thresh | trait %in% focal),]
  #dat[trait %in% focal,category:='aa_focal_diseases']
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
  dat[,category:=factor(category,levels=names(cat_levels))]
  dat[,trait:=factor(trait,levels=dat[order(category,delta,decreasing=TRUE),]$trait)]
  #ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point(aes(size=log10(n1))) + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  #coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc)
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }
  ggplot(dat,aes(x=trait,y=delta,colour=category,alpha=p.adj<fdr_thresh)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(title) + theme +
  xlab("Trait") + ylab("Change in basis loading from control") + guides(alpha=FALSE) +
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.3)) + scale_colour_manual("Category",values=cat_levels) +
  theme(axis.text.y=element_text(size=8))

}
basis.DT[,category:='basis']
## jia and jdm for talk with Lucy and Claire in January

lucy.DT <- res.DT[!category %in% c('lee_CD_prognosis','lyons_vasculitis')]
lucy.DT<-lucy.DT[category %in% lucy.DT[p.adj<0.05,]$category,]
at <- lucy.DT$category %>% unique
at <- at[!at %in% c('myogen','bowes_jia')]
library(RColorBrewer)
cols<-brewer.pal(length(lucy.DT$category %>% unique)+1, 'Paired')
names(cols) <- c('myogen','bowes_jia',at,'basis')
#lucy.DT[,category:=factor(category,levels=c('myogen','bowes_jia',at,'basis') %>% rev)]

basis.DT[,category:='basis']


pdf(file="~/tmp/lucy_forest_PC%d.pdf",paper="a4r",width=14,height=8,onefile=FALSE)
lapply(paste('PC',1:10,sep=''),function(pc){
    title <- sprintf("JIA/JDM subtypes %s",pc)
    forest_plot_focal(lucy.DT,pc=pc,focal=all.traits[c('myogen','bowes_jia')] %>% unlist,title=title)
})
dev.off()
lucy.no.DT <- lucy.DT[category!='bb_medications']
pdf(file="~/tmp/lucy_forest_no_med_PC%d.pdf",paper="a4r",width=14,height=8,onefile=FALSE)
forest_plot_focal(lucy.no.DT,pc='PC10',focal=all.traits[c('myogen','bowes_jia')] %>% unlist,title="JIA/JDM subtypes PC10 Med removed")
dev.off()

## for immunogenomics (just JIA and selected studies)

talk.DT <- res.DT[!category %in% c('lee_CD_prognosis','lyons_vasculitis','bb_medications','estrada_NMO','myogen','ferreira_asthma','psyc_consortium','tian_infectious_disease','astle_blood','bb_cancer')]
talk.DT<-talk.DT[category %in% talk.DT[p.adj<0.05,]$category,]
at <- talk.DT$category %>% unique
at <- at[!at %in% c('bowes_jia')]
library(RColorBrewer)
cols<-brewer.pal(length(talk.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c('bowes_jia',at,'basis')
cols['bowes_jia'] <- 'deeppink2'
cols['basis'] <- '#7c0799'
cols['bb_disease'] <- '#05af6e'
talk.DT[,trait:=gsub("^bb_","",trait)]




pc <- 'PC1'
title <- sprintf("JIA subtypes %s",pc)
pdf(file="~/tmp/jia_immunogenomics_pc1.pdf",paper="a4r",width=14,height=8,useDingbats=FALSE)
pp1 <- forest_plot_focal(talk.DT,pc=pc,focal=all.traits['bowes_jia'] %>% unlist,title=title,cat_levels=cols)
pp1
dev.off()

pc <- 'PC3'
title <- sprintf("JIA subtypes %s",pc)
pdf(file="~/tmp/jia_immunogenomics_pc3.pdf",paper="a4r",width=14,height=8,useDingbats=FALSE)
pp2 <- forest_plot_focal(talk.DT,pc=pc,focal=all.traits['bowes_jia'] %>% unlist,title=title,cat_levels=cols)
pp2
dev.off()
plot_grid(pp1,pp2)

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
  asthma = 'asthma'
)

## rename the bb traits keeping the category



talk.DT <- res.DT[!category %in% c('lee_CD_prognosis','lyons_vasculitis','bb_medications','estrada_NMO','myogen','ferreira_asthma','psyc_consortium','tian_infectious_disease','astle_blood','bb_cancer')]
talk.DT<-talk.DT[category %in% talk.DT[p.adj<0.05,]$category,]
at <- talk.DT$category %>% unique
at <- at[!at %in% c('bowes_jia')]
library(RColorBrewer)
cols<-brewer.pal(length(talk.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c('bowes_jia',at,'basis')
cols['bowes_jia'] <- 'deeppink2'
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



pc <- 'PC1'
title <- sprintf("JIA subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
pp1 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols)

pc <- 'PC1'
title <- sprintf("%s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
ppz <- forest_plot_focal_merge(talk.DT[category!='bowes_jia'],pc='PC1',focal=all.traits['bb_disease'] %>% unique,title=title,cat_levels=cols)

pc <- 'PC3'
title <- sprintf("JIA subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
pp3 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols)

pc <- 'PC9'
title <- sprintf("JIA subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
pp9 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols)


pdf(file="~/tmp/jia_immunogenomics_pc1.pdf",paper="a4r",width=14,height=8,useDingbats=FALSE)
pp1
dev.off()

pdf(file="~/tmp/jia_immunogenomics_pc3.pdf",paper="a4r",width=14,height=8,useDingbats=FALSE)
pp3
dev.off()

pdf(file="~/tmp/jia_immunogenomics_pc9.pdf",paper="a4r",width=14,height=8,useDingbats=FALSE)
pp9
dev.off()

pdf(file="~/tmp/immunogenomics_pc1.pdf",paper="a4r",width=14,height=8,useDingbats=FALSE)
ppz
dev.off()


## we are interested in antibody status in T1D

t1d.DT <- res.DT[!category %in% c('lee_CD_prognosis','lyons_vasculitis','bb_medications','estrada_NMO','myogen','ferreira_asthma','psyc_consortium','tian_infectious_disease','astle_blood','bb_cancer','bowes_jia')]
t1d.DT<-t1d.DT[category %in% t1d.DT[p.adj<0.05,]$category,]
at <- t1d.DT$category %>% unique
at <- at[!at %in% c('liley_t1d')]
library(RColorBrewer)
cols<-brewer.pal(length(t1d.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c('liley_t1d',at,'basis')
cols['liley_t1d'] <- 'deeppink2'
cols['basis'] <- '#7c0799'
cols['bb_disease'] <- '#05af6e'
t1d.DT[,trait:=gsub("^bb_","",trait)]
t1d.DT[,trait:=gsub("^z_","",trait)]

pc <- 'PC1'
title <- sprintf("T1D Ab subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
pp1 <- forest_plot_focal_merge(t1d.DT,pc=pc,focal=all.traits['liley_t1d'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols)

pc <- 'PC1'
title <- sprintf("T1D Ab subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
pp1 <- forest_plot_focal_merge(t1d.DT,pc=pc,focal=all.traits['liley_t1d'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols)


## KEN STUFF


forest_plot_focal_merge_ken <- function(proj.dat,basis.dat=basis.DT,pc,focal,title,cat_levels,fdr_thresh=0.05,theme=NA){
  dat <- proj.dat[variable==pc & (p.adj<fdr_thresh | trait %in% focal | trait %in% unlist(BB_LU)),]
  cat_levels <- cat_levels[names(cat_levels) %in% c(dat$category,'basis')]
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
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.5)) + scale_colour_manual("Category",values=cat_levels) +
  scale_linetype_manual(values=c('TRUE'=1,'FALSE'=2)) +
  theme(axis.text.y=element_text(size=12))

}

egpa.DT <- res.DT[!category %in% c('liley_t1d','lee_CD_prognosis','bb_medications','estrada_NMO','myogen','ferreira_asthma','psyc_consortium','bb_cancer','bowes_jia')]
egpa.DT<-egpa.DT[category %in% c(egpa.DT[p.adj<0.05,]$category,'lyons_vasculitis'),]
egpa.DT[trait=='mpo_Pos',trait:='EGPA:ANCA_MPO+']
egpa.DT[trait=='anca_Neg',trait:='EGPA:ANCA-']
egpa.DT[trait=='egpa',trait:='EGPA']
egpa.DT[trait=='mpo',trait:='VASC:ANCA_MPO+']
at <- egpa.DT$category %>% unique
at <- at[!at %in% c('lyons_vasculitis')]
library(RColorBrewer)
cols<-brewer.pal(length(egpa.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c('lyons_vasculitis',at,'basis')
cols['lyons_vasculitis'] <- 'deeppink2'
cols['basis'] <- '#7c0799'
cols['bb_disease'] <- '#05af6e'
cols['astle_blood'] <- 'orange2'
cols['tian_infectious_disease'] <- 'dodgerblue'
egpa.DT[,trait:=gsub("^bb_","",trait)]
egpa.DT[,trait:=gsub("^z_","",trait)]

## relable traits

all.traits <- traits<-split(egpa.DT$trait,egpa.DT$category) %>% lapply(.,unique)


# pdf("~/tmp/egap.pdf",onefile=TRUE,paper="a4r",width=14,height=8)
#
# for(pc in paste0('PC',1:10)){
#   title <- sprintf("EGPA subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
#   forest_plot_focal_merge_ken(egpa.DT,pc=pc,focal=all.traits['lyons_vasculitis'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols) %>% print
# }
# dev.off()

pdf("~/tmp/egpa_pc1.pdf",onefile=TRUE,paper="a4r",width=14,height=8)
pc <- 'PC1'
title <- sprintf("EGPA subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
forest_plot_focal_merge_ken(egpa.DT,pc=pc,focal=all.traits['lyons_vasculitis'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols)
dev.off()

pdf("~/tmp/egpa_pc4.pdf",onefile=TRUE,paper="a4r",width=14,height=8)
pc <- 'PC4'
title <- sprintf("EGPA subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
forest_plot_focal_merge_ken(egpa.DT,pc=pc,focal=all.traits['lyons_vasculitis'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols)
dev.off()

pdf("~/tmp/egpa_pc6.pdf",onefile=TRUE,paper="a4r",width=14,height=8)
pc <- 'PC6'
title <- sprintf("EGPA subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
forest_plot_focal_merge_ken(egpa.DT,pc=pc,focal=all.traits['lyons_vasculitis'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols)
dev.off()


pdf("~/tmp/egpa_pc10.pdf",onefile=TRUE,paper="a4r",width=14,height=8)
pc <- 'PC10'
title <- sprintf("EGPA subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
forest_plot_focal_merge_ken(egpa.DT,pc=pc,focal=all.traits['lyons_vasculitis'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols)
dev.off()


## what about PID ?

basis.DT[,category:='basis']

forest_plot_focal_merge_pid <- function(proj.dat,basis.dat=basis.DT,pc,focal,title,cat_levels,fdr_thresh=0.05,theme=NA){
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
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.5)) + scale_colour_manual("Category",values=cat_levels,labels=c('PID','UKBB SR Disease','Astle','Basis')) +
  scale_linetype_manual(values=c('TRUE'=1,'FALSE'=2)) +
  theme(axis.text.y=element_text(size=12))

}


pid.DT <- res.DT[!category %in% c('tian_infectious_disease','lyons_vasculitis','liley_t1d','lee_CD_prognosis','bb_medications','estrada_NMO','myogen','ferreira_asthma','psyc_consortium','bb_cancer','bowes_jia')]
pid.DT<-pid.DT[category %in% c(pid.DT[p.adj<0.05,]$category,'ad-pid'),]
at <- pid.DT$category %>% unique
at <- at[!at %in% c('ad-pid')]
library(RColorBrewer)
cols<-brewer.pal(length(pid.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c('ad-pid',at,'basis')
cols['ad-pid'] <- 'deeppink2'
cols['basis'] <- '#7c0799'
cols['bb_disease'] <- '#05af6e'
cols['astle_blood'] <- 'orange2'
cols['tian_infectious_disease'] <- 'cyan'
pid.DT[,trait:=gsub("^bb_","",trait)]
pid.DT[,trait:=gsub("^z_","",trait)]
pid.DT[trait=="ABDEF",trait:='Ab Deficiency PID']

all.traits <- traits<-split(pid.DT$trait,pid.DT$category) %>% lapply(.,unique)

#pc <- 'PC1'
#title <- sprintf("PID %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
#pp1 <- forest_plot_focal_merge_pid(pid.DT,pc=pc,focal=all.traits['ad-pid'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols)

#pdf("~/tmp/pid.pdf",onefile=TRUE,paper="a4r",width=14,height=8)

#for(pc in paste0('PC',1:10)){
#  title <- sprintf("PID subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
#  forest_plot_focal_merge_pid(pid.DT,pc=pc,focal=all.traits['ad-pid'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols) %>% print
#}
#dev.off()

## do a  nice plot of pc3_plot

pdf("~/tmp/pid_pc3.pdf",onefile=TRUE,paper="a4r",width=14,height=8)
pc <- 'PC3'
title <- sprintf("Antibody Deficiency Associated PID %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
forest_plot_focal_merge_pid(pid.DT,pc=pc,focal=all.traits['ad-pid'] %>% unlist %>% gsub("^z_","",.),title=title,cat_levels=cols)
dev.off()


## clustering by weighted Z score

imp.dt<-data.table(pc=colnames(pc.emp$x),weight=summary(pc.emp)[['importance']][2,])
M <- merge(res.DT[category=='bowes_jia',.(pc=variable,trait,Z)],imp.dt,by='pc')
M[,hval:=Z * weight]
M<-melt(M,id.vars=c('trait','pc'),measure.vars='hval') %>% dcast(.,"trait~pc")
m<-as.matrix(M[,-1])
rownames(m) <- M$trait
dist(m) %>% hclust %>% plot

hinks <- fread("/home/ob219/tmp/hinks_hla_correlation.csv")
hi <- as.matrix(hinks[,-1])
rownames(hi) <- hinks$trait
#pdf("~/tmp/hinks_hla_genetic_correlation.pdf")
dist(hi) %>% hclust %>% plot(.,,main="",sub="")
#dev.off()
