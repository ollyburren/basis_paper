## thesis stuff

## for the thesis we will not mention external datasets

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

RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}



filt.DT <- res.DT[category %in% c('bb_disease','bb_cancer'),]
filt.DT[category=='bb_medications',category:='Medications']
filt.DT[category=='bb_cancer',category:='Cancer']
filt.DT[category=='bb_disease',category:='Disease']
dat <- filt.DT[trait %in% filt.DT[(p.adj<0.05  | trait %in% unlist(BB_LU)),]$trait,]
dat[,trait:=gsub("^bb_","",trait) %>% gsub("\\."," ",.) %>% sapply(.,simpleCap)]
dat[,ci:=1.96 * sqrt(variance)]
dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
## try a heatmap but manual clustering
mat <- melt(dat[,.(pc=variable,trait,delta)],id.vars=c('pc','trait')) %>% dcast(pc~trait)
mdt <- as.matrix(mat[,-1])
rownames(mdt) <- mat$pc
hc <- dist(t(mdt)) %>% hclust
dat[,trait:=factor(trait,levels=hc$labels[hc$order])]
dat[,variable:=factor(variable,levels=paste('PC',1:10,sep=''))]
dat[,is.sig:='']
dat[p.adj<0.05,is.sig:='*']
pp1 <- ggplot(dat[variable!='PC11'],aes(x=variable,y=trait,fill=delta,label=is.sig)) +
geom_tile() + geom_text() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
scale_fill_gradient2("Difference\nfrom control") +
ylab("UKBB Self reported disease/cancer") +
xlab("Principal Component")

save_plot(pp1,file="~/tmp/ukbbsrd.pdf",base_height=8,base_aspect_ratio=0.9)

## separate plot for medications

filt.DT <- res.DT[category %in% c('bb_medications'),]
filt.DT[category=='bb_medications',category:='Medications']
filt.DT[category=='bb_cancer',category:='Cancer']
filt.DT[category=='bb_disease',category:='Disease']
## ignore PC10 as so many astham treatments !

dat <- filt.DT[!variable %in% c('PC10','PC11'),]
dat <- filt.DT[variable %in% c('PC10'),]
dat <- dat[trait %in% dat[(p.adj<0.05  | trait %in% unlist(BB_LU)),]$trait,]
dat[,trait:=gsub("^bb_","",trait) %>% gsub("\\."," ",.) %>% sapply(.,simpleCap)]
dat <- dat[trait!='Vitamin C Product',]
dat[,ci:=1.96 * sqrt(variance)]
dat[,c('lower','upper'):=list(delta-ci,delta+ci)]

## try a heatmap but manual clustering
mat <- melt(dat[,.(pc=variable,trait,delta)],id.vars=c('pc','trait')) %>% dcast(pc~trait)
mdt <- as.matrix(mat[,-1])
rownames(mdt) <- mat$pc
hc <- dist(t(mdt)) %>% hclust
dat[,trait:=factor(trait,levels=hc$labels[hc$order])]
dat[,variable:=factor(variable,levels=paste('PC',1:11,sep=''))]
dat[,is.sig:='']
dat[p.adj<0.05,is.sig:='*']
pp1 <- ggplot(dat,aes(x=variable,y=trait,fill=delta,label=is.sig)) +
geom_tile() + geom_text() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
scale_fill_gradient2("Difference\nfrom control") +
ylab("UKBB Self reported medication") +
xlab("Principal Component")
save_plot(pp1,file="~/tmp/ukbbmed.pdf",base_height=8,base_aspect_ratio=0.9)




pheatmap(mdt)
library(pheatmap)





library(cowplot)

all.traits <- traits<-split(res.DT$trait,res.DT$category) %>% lapply(.,unique)

## add in basis traits for comparison

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
basis.DT[,category:='Basis']








forest_plot_focal_merge <- function(proj.dat,basis.dat=basis.DT,pc,title,cat_levels,fdr_thresh=0.05,theme=NA){
  dat <- proj.dat[variable==pc & (p.adj<fdr_thresh  | trait %in% unlist(BB_LU)),]
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  ## define a thing called trait label
  dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
  for(i in seq_along(BB_LU)){
    tra <- names(BB_LU)[i]
    dat[trait==tra,trait:=BB_LU[[i]]]
  }
  #dat[,category:=factor(category,levels=names(cat_levels))]
  ## sort focal trait first
  ## next do the rest
  nf.dt <- dat[,.(trait,delta,category,n1)][order(category,decreasing=FALSE),][!duplicated(trait),]
  nf.dt[is.na(n1),tl:=trait]
  nfoc <- nf.dt[order(category,delta,decreasing=TRUE),.(trait,n1)]
  trait.order <- nfoc$trait
  dat[,trait:=factor(trait,levels=trait.order)]
  ldat <- nfoc[,tl:=sprintf("%s",trait)]
  ldat[is.na(n1),tl:=trait]
  levels(dat$trait) <- ldat$tl
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }

  ggplot(dat,aes(x=trait,y=delta,colour=category,alpha=p.adj<fdr_thresh,lty=p.adj<fdr_thresh)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(title) + theme +
  xlab("Trait") + ylab("Change in basis score from control") + guides(alpha=FALSE,lty=FALSE) +
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.5)) + #scale_colour_manual("Category",values=cat_levels,labels=c('JIA','UKBB SR Disease','Basis')) +
  scale_linetype_manual(values=c('TRUE'=1,'FALSE'=2)) +
  theme(axis.text.y=element_text(size=12))
}
library(RColorBrewer)
at <- filt.DT$category %>% unique

cols <- c(Disease="#05af6e",Basis="#D95F02",Medications="#7570B3",Cancer="#7c0799")
#cols<-brewer.pal(length(at)+1, 'Dark2')
#names(cols) <- c(at,'Basis')
#cols['Basis'] <- '#7c0799'
#cols['Disease'] <- '#05af6e'
#filt.DT[grepl("^bb",trait),trait:=sub("^bb\\_","",trait)]

forest_plot_focal_merge(filt.DT,pc='PC10',title="BLAH",cat_levels=cols)

forest_plot_focal_merge(filt.DT,pc='PC8',title="BLAH",cat_levels=cols)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}



forest_plot_focal_merge <- function(proj.dat,basis.dat=basis.DT,pc,focal,title,cat_levels,fdr_thresh=0.05,labels,theme=NA,noshowbasis=FALSE){
  dat <- proj.dat[variable==pc & (p.adj<fdr_thresh | trait %in% focal | trait %in% unlist(BB_LU)),]
  #dat[trait %in% focal,category:='aa_focal_diseases']
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  ## define a thing called trait label
  dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
  for(i in seq_along(BB_LU)){
    tra <- names(BB_LU)[i]
    dat[trait==tra,trait:=BB_LU[[i]]]
  }
  dat[,trait:=gsub("\\."," ",trait) %>% sapply(.,simpleCap)]
  dat[trait=='Sys',trait:='sys']
  dat[,category:=factor(category,levels=names(cat_levels))]
  ## sort focal trait first
  foc.dt <- dat[trait %in% focal,.(trait,delta,category,n1)]
  #foc.dt[,tl:=sprintf("%s(%s)",trait,format(n1,big.mark=",",scientific=FALSE))]
  foc <- foc.dt[order(delta,decreasing=TRUE),.(trait,n1)]
  ## next do the rest
  nf.dt <- dat[!trait %in% focal,.(trait,delta,category,n1)][order(category,decreasing=FALSE),][!duplicated(trait),]
  #nf.dt[,tl:=sprintf("%s(%s)",trait,format(n1,big.mark=",",scientific=FALSE))]
  nf.dt[is.na(n1),tl:=trait]
  ## sort diseases and basis together so appear merged then add other traits
  bd <- nf.dt[category %in% c('Basis','Disease')]
  bd <- bd[order(delta,decreasing=TRUE),.(trait,n1)]
  nbd <- nf.dt[!category %in% c('Basis','Disease','aa_focal_diseases')]
  nbd <- nbd[order(category,delta,decreasing=TRUE),.(trait,n1)]
  #nfoc <- nf.dt[order(category,delta,decreasing=TRUE),.(trait,n1)]
  trait.order <- c(nbd$trait,bd$trait,foc$trait)
  #dat[!duplicated(trait),][order(trait %in% focal,delta,decreasing=TRUE),]$trait
  #dat[,trait:=factor(trait,levels=dat[order(delta,!trait %in% focal,decreasing=TRUE),]$trait)]
  dat[,trait:=factor(trait,levels=trait.order)]
  #rbind(nfoc,foc)[,tl:=sprintf("%s(%s)",trait,format(n1,big.mark=",",scientific=FALSE))]
  #ldat <- rbind(nfoc,foc)[,tl:=sprintf("%s",trait)]
  #ldat[is.na(n1),tl:=trait]
  #levels(dat$trait) <- ldat$tl
  #ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point(aes(size=log10(n1))) + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  #coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc)
  ## next alter the labels to include sample size
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }
  if(noshowbasis==TRUE)
    dat <- dat[category!='Basis',]
  ggplot(dat,aes(x=trait,y=delta,colour=category,alpha=p.adj<fdr_thresh,lty=p.adj<fdr_thresh)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(title) + theme +
  xlab("Trait") + ylab("Change in basis score from control") + guides(alpha=FALSE,lty=FALSE) +
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.8)) + scale_colour_manual("Category",values=cat_levels,labels=labels) +
  scale_linetype_manual(values=c('TRUE'=1,'FALSE'=2)) +
  theme(axis.text.y=element_text(size=12))
}


filt.DT <- res.DT[category %in% c('bb_disease','bb_cancer','bb_medications','bowes_jia'),]
filt.DT[category=='bb_medications',category:='Medications']
filt.DT[category=='bb_cancer',category:='Cancer']
filt.DT[category=='bb_disease',category:='Disease']
filt.DT[category=='bowes_jia',category:='JIA']
library(cowplot)
filt.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
at <- filt.DT$category %>% unique
at <- at[!at %in% c('bowes_jia')]
library(RColorBrewer)
cols<-brewer.pal(length(filt.DT$category %>% unique), 'Dark2')
names(cols) <- at
cols['JIA'] <- 'deeppink2'
cols['Basis'] <- '#7c0799'
cols['Cancer'] <- 'dodgerblue1'
#cols['Disease'] <- '#05af6e'
cols['Medications'] <- '#05af6e'
filt.DT[,trait:=gsub("^bb_","",trait)]
filt.DT[,trait:=gsub("^jia_","",trait)]

## appendix plots
for(pc in paste('PC',c(1,5,6),sep='')){
  labs1 <- c('JIA',"UKBB\nSRM","UKBB\nSRD","UKBB\nSRC",'Basis')
  title <- pc
  fname <- file.path('~/tmp/pc_app/',sprintf("appendix_%s.pdf",pc))
  pp1 <- forest_plot_focal_merge(filt.DT,pc=pc,focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols,labels=labs1)
  save_plot(pp1 + theme(legend.position="bottom",axis.text.y=element_text(size=13)),file=fname,base_width=8,base_height=9,useDingbats=FALSE)
}

for(pc in paste('PC',c(2,3,4,7,8,9),sep='')){
  labs1 <- c('JIA',"UKBB\nSRM","UKBB\nSRD",'Basis')
  title <- pc
  fname <- file.path('~/tmp/pc_app/',sprintf("appendix_%s.pdf",pc))
  pp1 <- forest_plot_focal_merge(filt.DT,pc=pc,focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols,labels=labs1)
  save_plot(pp1 + theme(legend.position="bottom",axis.text.y=element_text(size=13)),file=fname,base_width=8,base_height=9,useDingbats=FALSE)
}

labs1 <- c('JIA',"UKBB\nSRM","UKBB\nSRD",'Basis')
title <- 'PC10'
fname <- file.path('~/tmp/pc_app/',sprintf("appendix_%s.pdf",pc))
pp1 <- forest_plot_focal_merge(filt.DT,pc=pc,focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols,labels=labs1)
save_plot(pp1 + theme(legend.position="bottom",axis.text.y=element_text(size=9)),file=fname,base_width=8,base_height=9,useDingbats=FALSE)

pp2 <- forest_plot_focal_merge(filt.DT,pc='PC3',focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols,labels=labs1)
save_plot(pp2 + theme(legend.position="bottom",axis.text.y=element_text(size=13)),file="~/tmp/pc3_thesis.pdf",base_height=8)

pp3 <- forest_plot_focal_merge(filt.DT,pc='PC6',focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols,labels=labs1)
save_plot(pp3 + theme(legend.position="bottom",axis.text.y=element_text(size=13)),file="~/tmp/pc6_thesis.pdf",base_height=8)


pp4 <- forest_plot_focal_merge(filt.DT,pc='PC9',focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols,labels=labs1)
save_plot(pp4 + theme(legend.position="bottom",axis.text.y=element_text(size=13)),file="~/tmp/pc9_thesis.pdf",base_height=8)

pp5 <- forest_plot_focal_merge(filt.DT,pc='PC10',focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols,labels=labs1)
save_plot(pp5 + theme(legend.position="bottom",axis.text.y = element_text(size=11,angle = 0, vjust = 1, hjust=1)),file="~/tmp/pc10_thesis.pdf",base_height=10,base_aspect=0.9)


#pp1 + theme(legend.position="bottom")

filt_med.DT <- res.DT[category %in% c('bb_medications','bowes_jia'),]
filt_med.DT[category=='bb_medications',category:='Medications']
filt_med.DT[category=='bb_cancer',category:='Cancer']
filt_med.DT[category=='bb_disease',category:='Disease']
filt_med.DT[category=='bowes_jia',category:='JIA']
library(cowplot)
filt_med.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
at <- filt_med.DT$category %>% unique
at <- at[!at %in% c('bowes_jia')]
library(RColorBrewer)
cols<-brewer.pal(length(filt_med.DT$category %>% unique), 'Dark2')
names(cols) <- at
cols['JIA'] <- 'deeppink2'
cols['Basis'] <- '#7c0799'
cols['Medications'] <- '#05af6e'
filt_med.DT[,trait:=gsub("^bb_","",trait)]
filt_med.DT[,trait:=gsub("^jia_","",trait)]

labs2 <- c('JIA','UKBB SR Medications')

#title <- sprintf("JIA subtypes %s (Variance Explained %0.1f%%)",pc,(summary(pc.emp)[['importance']][2,][pc] * 100) %>% signif(.,digits=2))
title <- ''
pp2 <- forest_plot_focal_merge(filt_med.DT,pc=pc,focal=all.traits['bowes_jia'] %>% unlist %>% gsub("^jia_","",.),title=title,cat_levels=cols,labels=labs2,noshowbasis=TRUE)

plot_grid(pp1 + theme(legend.position="bottom"),pp2 + theme(legend.position="bottom"),ncol=2,labels="auto")
