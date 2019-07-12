library(cowplot)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/12_07_19_0619_summary_results.RDS'
#RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)



## add in basis traits for comparison

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
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

N_LU <- list(
  anca_Neg = 'EGPA_Lyons:ANCA-',
  mpo_Pos = 'EGPA_Lyons:MPO+',
  #mpo = 'AAV_Lyons:MPO+',
  mpo_meta = 'AAV_Wong:MPO+',
  pr3_meta = 'AAV_Wong:PR3+'
  #egpa = 'EGPA_Lyons:Overall'
)

category.foc <- c('Focal')
talk.DT <- res.DT[trait %in% names(N_LU),category:='Focal']
talk.DT <- res.DT[category %in% c('bb_disease','Focal'),]
talk.DT<-talk.DT[(category %in% talk.DT[p.adj<0.05,]$category) | category==category.foc,]
## rename with nicer labels
for(i in seq_along(N_LU)){
  talk.DT[trait==names(N_LU)[i],trait:=N_LU[[i]]]
}
at <- talk.DT$category %>% unique
at <- at[!at %in% c(category.foc)]
library(RColorBrewer)
cols<-brewer.pal(length(talk.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c(category.foc,at,'basis')
cols[category.foc] <- 'deeppink2'
cols['basis'] <- '#7c0799'
cols['bb_disease'] <- '#05af6e'
talk.DT[,trait:=gsub("^bb_SRD:","",trait)]

all.traits <- traits<-split(talk.DT$trait,talk.DT$category) %>% lapply(.,unique)


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

pc <- 'PC1'

pp1 <- forest_plot_focal_merge(talk.DT[trait %in% c('AAV_Wong:MPO+','AAV_Wong:PR3+') | category!='Focal',],pc='PC1',focal=all.traits[category.foc] %>% unlist,title=pc,cat_levels=cols,fdr_thresh=0.01)
save_plot(file="~/tmp/smith_090719_forest.pdf",pp1,base_height=9)

## we particularly care about PC1 for AAV MPO vs PR3


pdf(file="~/tmp/smith_090719.pdf",paper="a4r",onefile=TRUE)
lapply(paste('PC',1:12,sep=''),function(pc){
  forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits[category.foc] %>% unlist,title=pc,cat_levels=cols,fdr_thresh=0.01)
})
dev.off()



## want a heatmap of all the changes

FDR <- 0.01

talk.DT[,significant:=p.adj<FDR]
talk.DT[,pc:=factor(variable,levels=paste('PC',1:12,sep=''))]
talk.DT[,is.sig:='']
talk.DT[p.adj<FDR,is.sig:='*']
talk.DT[,trait.label:=sprintf("%s (%d)",trait,n1)]

pp1 <- ggplot(talk.DT[category=='Focal' & pc!='PC12',],aes(x=pc,y=trait.label,fill=delta,label=is.sig)) + geom_tile() +
geom_text(size=10) +
scale_fill_gradient2("Difference\nfrom control") + xlab("Principal Component") +
ylab("Disease/Subtype") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

save_plot(file="~/tmp/smith_090719_vasc_heatmap.pdf",pp1,base_width=11)


## next want to plot forest plots for astle data (top 13 cell types) with EGPA and MPO plotted on top


N_LU <- list(
  anca_Neg = 'EGPA_Lyons:ANCA-',
  mpo_Pos = 'EGPA_Lyons:MPO+',
  mpo_meta = 'AAV_Wong:MPO+'
)

category.foc <- c('Focal')
talk.DT <- res.DT[trait %in% names(N_LU),category:='Focal']
main.bc <- list(platelet=c('pdw','mpv','plt'),rbc=c('irf','ret','rdw','hct','mch'),lymphoid=c('mono','baso','eo','neut','lymph'))
talk.DT <- res.DT[category %in% c('Focal') | trait %in% unlist(main.bc) ,]
talk.DT<-talk.DT[(category %in% talk.DT[p.adj<0.05,]$category) | category==category.foc,]
## rename with nicer labels
for(i in seq_along(N_LU)){
  talk.DT[trait==names(N_LU)[i],trait:=N_LU[[i]]]
}
at <- talk.DT$category %>% unique
at <- at[!at %in% c(category.foc)]
library(RColorBrewer)
cols<-brewer.pal(length(talk.DT$category %>% unique)+1, 'Dark2')
names(cols) <- c(category.foc,at,'basis')
cols[category.foc] <- 'deeppink2'
cols['basis'] <- '#7c0799'
cols['bb_disease'] <- '#05af6e'
talk.DT[,trait:=gsub("^bb_SRD:","",trait)]

all.traits <- traits<-split(talk.DT$trait,talk.DT$category) %>% lapply(.,unique)

pc <- 'PC4'
pp4 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits[category.foc] %>% unlist,title=pc,cat_levels=cols,fdr_thresh=0.01) +
theme(legend.position = "none")
#save_plot(file="~/tmp/smith_090719_forest.pdf",pp1,base_height=9)
pc <- 'PC7'
pp7 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits[category.foc] %>% unlist,title=pc,cat_levels=cols,fdr_thresh=0.01) +
theme(legend.position = "none")
pc <- 'PC10'
pp10 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits[category.foc] %>% unlist,title=pc,cat_levels=cols,fdr_thresh=0.01) +
theme(legend.position = "none")
pc <- 'PC11'
pp11 <- forest_plot_focal_merge(talk.DT,pc=pc,focal=all.traits[category.foc] %>% unlist,title=pc,cat_levels=cols,fdr_thresh=0.01) +
theme(legend.position = "none")

plot_grid(pp4,pp7,pp10,pp11,nrow=2) %>% save_plot(file="~/tmp/astle_vasc.pdf",.,base_height=7,base_aspect_ratio=2)
