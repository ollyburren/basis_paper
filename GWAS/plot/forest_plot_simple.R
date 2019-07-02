## FOREST PLOTS FOR particular diseases
library(cowplot)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/24_06_19_0619_summary_results.RDS'
#RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)

all.traits <- traits<-split(res.DT$trait,res.DT$category) %>% lapply(.,unique)

## add in basis traits for comparison

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
basis.DT[,category:='basis']


forest_plot <- function(proj.dat,basis.dat=basis.DT,pc,fdr_thresh=0.05,theme=NA){
  message(fdr_thresh)
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

everything.DT <- res.DT

pdf(file="~/tmp/everything_forest_0619_pah.pdf",width=14,height=20,onefile=TRUE)
lapply(paste('PC',1:12,sep=''),function(pc){
  if(pc != 'PC122'){
    forest_plot(everything.DT,pc=pc)
  }else{
    forest_plot(everything.DT,pc=pc,theme=theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3),axis.text.y=element_text(size=6)))
  }
})
dev.off()

everything.DT[trait=='jia_case_19',trait:='jia_combined_19']

pdf(file="~/tmp/everything_forest_0619_1FDR.pdf",width=14,height=20,onefile=TRUE)
lapply(paste('PC',1:12,sep=''),function(pc){
  if(pc != 'PC122'){
    forest_plot(everything.DT,pc=pc,fdr_thresh=0.01)
  }else{
    forest_plot(everything.DT,pc=pc,fdr_thresh=0.01,theme=theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3),axis.text.y=element_text(size=6)))
  }
})
dev.off()


## tSNE stuff - not used but here for future refeence.
## perplexity here is difficult to tune due to the small number of dimensions
if(FALSE){
library(Rtsne)

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
) %>% unlist
bb.keep<-paste('bb_SRD',BB_LU,sep=':')
bb.keep<-c(bb.keep,'bb_SRD:ankylosing.spondylitis')

res.DT <- res.DT[trait!='bb_SRM:vitamin.c.product']

dat <- melt(res.DT[trait %in% (res.DT[p.adj<1,]$trait %>% unique),],id.vars=c('trait','variable'),measure.vars='delta') %>% dcast(trait~variable)

mat <- dat[,-1] %>% as.matrix
rownames(mat) <- dat$trait
set.seed(9)
tsne_model_1 = Rtsne(dist(mat), check_duplicates=FALSE, pca=TRUE, perplexity=10, theta=0.5, dims=2)
d_tsne_1 = data.table(tsne_model_1$Y)

lookup <- merge(data.table(trait=dat$trait,order=1:nrow(dat)),unique(res.DT[,.(trait,category)]),by='trait')
pdat <- cbind(d_tsne_1,lookup)

library(cowplot)
library(ggrepel)

ggplot(pdat[category %in% c('lyons_egpa','bowes_jia_2019','estrada_NMO') | trait %in% bb.keep], aes(x=V1, y=V2,color=category,label=trait)) +
  geom_point() + geom_text_repel()
}
