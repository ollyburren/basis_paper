## analyse jia projections
library(cowplot)
DATA.DIR <- '/home/ob219/share/as_basis/ichip/individual_data/individual_proj'
all.files <- list.files(path=DATA.DIR,pattern="*.RDS",full.names=TRUE)
BASIS_FILE <- '/home/ob219/share/as_basis/ichip/support/basis_ic.RDS'
VARIANCE_FILE <- '/home/ob219/share/as_basis/ichip/support/shrink_av_ichip.RDS'

res <- lapply(all.files,function(f){
  trait <- basename(f) %>% gsub("\\_projection\\.RDS","",.)
  dat <- readRDS(f)
  t.DT <- data.table(ind=rownames(dat),dat)
  melt.DT <- melt(t.DT,id.vars='ind')
  melt.DT[,trait:=trait]
}) %>% rbindlist

res <- res[!trait %in%  c('JIA_undefined','JIA_unknown'),]

## get mean and variance across trait and pc

summ.DT <- res[,list(mean.load=mean(value),var.load=var(value)),by=c('trait','variable')]

pc.emp <- readRDS(BASIS_FILE)
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
tmp <- basis.DT[trait=='control',] %>% t
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))

var.DT <- readRDS(VARIANCE_FILE)
summ.DT[,variable:=factor(variable,levels=paste0('PC',1:13))]
summ.DT <- merge(summ.DT,ctrl.DT,by='variable')
summ.DT <- merge(summ.DT,var.DT,by.x='variable',by.y='pc')
summ.DT[,Z:=(mean.load-control.loading)/sqrt(mfactor)]
summ.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
#bb.DT.m[,p.adj:=p.adjust(p.value),by='variable']
summ.DT[,p.adj:=p.adjust(p.value),by='variable']
summ.DT[,variable:=factor(variable,levels=paste0('PC',1:13))]
summ.DT[,short.trait:=substr(trait,1,15),]


pd <- position_dodge(0.1)
pa <- ggplot(summ.DT[!trait %in% c('jiaUnA','jiamissing'),],aes(x=variable,y=mean.load-control.loading,group=trait,col=short.trait)) + geom_point(position=pd) +
geom_line(position=pd) + guides(size=FALSE)

jia.sum <- readRDS("~/share/as_basis/GWAS/tmp/jia_plot.RDS")
pb <- ggplot(jia.sum[!trait %in% c('jiaUnA','jiamissing'),],aes(x=variable,y=value-control.loading,group=trait,col=short.trait)) + geom_point(position=pd) +
geom_line(position=pd) + guides(size=FALSE)
## comparison between summary and individual data
plot_grid(pa + ggtitle("mean.individual"),pb + ggtitle("summary stats"), nrow=2)


## next do t.test of sys and era and the rest to see which PC's are important to discern between

g1 <- c('ERA','systemic')
diff <- lapply(paste('PC',1:13,sep=''),function(y){
  t.tmp <- t.test(res[variable==y & trait %in% g1,]$value,res[variable==y & !trait %in% g1,]$value)
  data.table(pc=y,p.value=t.tmp$p.value,tstat=t.tmp$statistic)
}) %>% rbindlist()

diff[,p.adj:=p.adjust(p.value)]

## other class is
g1 <- c('jPsA','ext_oligo')
g2 <- c('RFneg_poly','RFpos_poly','pers_oligo')
diff2 <- lapply(paste('PC',1:13,sep=''),function(y){
  t.tmp <- t.test(res[variable==y & trait %in% g1,]$value,res[variable==y & trait %in% g2,]$value)
  data.table(pc=y,p.value=t.tmp$p.value,tstat=t.tmp$statistic)
}) %>% rbindlist()

diff2[,p.adj:=p.adjust(p.value)]
