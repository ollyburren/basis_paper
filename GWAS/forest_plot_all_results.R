## FOREST PLOTS FOR particular diseases
library(cowplot)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/19_12_18_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)

all.traits <- traits<-split(res.DT$trait,res.DT$category) %>% lapply(.,unique)

## add in basis traits for comparison

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
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
  ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(title) + theme +
  xlab("Trait") + ylab("Change in basis loading from control")
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
  ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc) + theme +
  xlab("Trait") + ylab("Change in basis loading from control")
}

pdf(file="~/tmp/basis_results.pdf",paper="a4r",width=14,height=8)
lapply(paste('PC',1:10,sep=''),function(pc){
  if(pc != 'PC10'){
    forest_plot(res.DT,pc=pc)
  }else{
    forest_plot(res.DT,pc=pc,theme=theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3),axis.text.y=element_text(size=6)))
  }
})
dev.off()
