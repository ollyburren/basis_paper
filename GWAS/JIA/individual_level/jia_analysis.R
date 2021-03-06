library(cowplot)

## analyse jia projections

DATA.DIR <- '/home/ob219/share/as_basis/GWAS/individual_data/individual_proj'
all.files <- list.files(path=DATA.DIR,pattern="*.RDS",full.names=TRUE)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'

res <- lapply(all.files,function(f){
  trait <- basename(f) %>% gsub("\\_projection\\.RDS","",.)
  dat <- readRDS(f)
  t.DT <- data.table(ind=rownames(dat),dat)
  melt.DT <- melt(t.DT,id.vars='ind')
  melt.DT[,trait:=trait]
}) %>% rbindlist

res <- res[!trait %in%  c('jiaUnA','jiamissing','raj_cd14','raj_cd4'),]

t.test(res[variable=='PC3' & trait %in% c('jiaERA','jiasys')]$value,res[variable=='PC3' & !trait %in% c('jiaERA','jiasys')]$value)
t.test(res[variable=='PC3' & trait=='jiaERA']$value,res[variable=='PC3' & !trait %in% c('jiaERA','jiasys')]$value)
t.test(res[variable=='PC3' & trait=='jiasys']$value,res[variable=='PC3' & !trait %in% c('jiaERA','jiasys')]$value)

## look at all of them

traits <- res$trait %>% unique

all.compare <- lapply(paste('PC',1:10,sep=''),function(PC){
  message(PC)
  lapply(traits,function(tra){
    lapply(traits,function(tra2){
      tt <- t.test(res[variable==PC & trait==tra]$value,res[variable==PC & trait==tra2]$value)
      data.table(pc=PC,trait1=tra,trait2=tra2,p=tt$p.value,t.stat=tt$statistic)
    }) %>% rbindlist
  }) %>% rbindlist
}) %>% rbindlist

all.compare <- all.compare[trait1 != trait2,]
## get rid of reciprocal comparisons
all.compare <- all.compare[which(!duplicated(abs(t.stat))),]
all.compare[abs(t)]
all.compare[,fdr:=p.adjust(p,method="fdr")]
all.compare[fdr<0.05,]

all.compare.rest <- lapply(paste('PC',1:10,sep=''),function(PC){
  message(PC)
  lapply(traits,function(tra){
      tt <- t.test(res[variable==PC & trait==tra]$value,res[variable==PC & trait!=tra]$value)
      data.table(pc=PC,trait1=tra,p=tt$p.value,t.stat=tt$statistic)
  }) %>% rbindlist
}) %>% rbindlist

all.compare.rest[,fdr:=p.adjust(p,method="bonferroni"),by=pc]
all.compare.rest[fdr<0.05,]

## get mean and variance across trait and pc

summ.DT <- res[,list(mean.load=mean(value),var.load=var(value)),by=c('trait','variable')]

pc.emp <- readRDS(BASIS_FILE)
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
tmp <- basis.DT[trait=='control',] %>% t
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))

#var.DT <- readRDS(VARIANCE_FILE)
summ.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
summ.DT <- merge(summ.DT,ctrl.DT,by='variable')
#summ.DT <- merge(summ.DT,var.DT,by.x='variable',by.y='pc')
summ.DT[,Z:=(mean.load-control.loading)/sqrt()]
summ.DT[,Z:=(mean.load-control.loading)/sqrt(mfactor)]
summ.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
#bb.DT.m[,p.adj:=p.adjust(p.value),by='variable']
summ.DT[,p.adj:=p.adjust(p.value),by='variable']
summ.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
summ.DT[,short.trait:=substr(trait,1,15),]

summ.DT[,Subtype:=gsub("^jia","",trait)]
pd <- position_dodge(0.1)
#pa <- ggplot(summ.DT[!trait %in% c('jiaUnA','jiamissing','raj_cd14','raj_cd4'),],aes(x=variable,y=mean.load-control.loading,group=Subtype,col=Subtype)) + geom_point(position=pd) +
#geom_line(position=pd) + guides(size=FALSE) + xlab("Principal Component") + ylab(expression(Delta~"Control Loading")) +

pa <- ggplot(summ.DT[!trait %in% c('jiaUnA','jiamissing','raj_cd14','raj_cd4'),],aes(x=variable,y=mean.load-control.loading,group=Subtype,col=Subtype)) + geom_point(size=2,position=pd) +
geom_line(position=pd) + ylab(expression(Delta*"Control Loading")) + xlab("Principal Component") + geom_hline(yintercept=0,color="black") +
background_grid(major = "xy", minor = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
save_plot(pa,file="~/tmp/ind_jia_line.pdf",base_aspect=1.3)


jia.sum <- readRDS("~/share/as_basis/GWAS/tmp/jia_plot.RDS")
jia.sum[,Subtype:=gsub("jia\\_","",trait)]
pb <- ggplot(jia.sum[!trait %in% c('jiaUnA','jiamissing'),],aes(x=variable,y=value-control.loading,group=Subtype,col=Subtype)) + geom_point(position=pd) +
geom_line(position=pd) + guides(size=FALSE) +  xlab("Principal Component") + ylab(expression(Delta~"Control Loading"))
## comparison between summary and individual data
plot_grid(pa + ggtitle("Genotypes") + geom_hline(yintercept=0,color='black'),pb + ggtitle("Summary statistics") + geom_hline(yintercept=0,color='black'), nrow=2)
dev.print(pdf,"~/tmp/gt_vs_summ.pdf")


## next do t.test of sys and era and the rest to see which PC's are important to discern between

g1 <- c('jiaERA','jiasys')
diff <- lapply(paste('PC',1:11,sep=''),function(y){
  t.tmp <- t.test(res[variable==y & trait %in% g1,]$value,res[variable==y & !trait %in% g1,]$value)
  data.table(pc=y,p.value=t.tmp$p.value,tstat=t.tmp$statistic)
}) %>% rbindlist()

diff[,p.adj:=p.adjust(p.value)]

## other class is
g1 <- c('jiaPsA','jiaEO')
g2 <- c('jiaRFneg','jiaRFpos','jiaPO')
diff2 <- lapply(paste('PC',1:11,sep=''),function(y){
  t.tmp <- t.test(res[variable==y & trait %in% g1,]$value,res[variable==y & trait %in% g2,]$value)
  data.table(pc=y,p.value=t.tmp$p.value,tstat=t.tmp$statistic)
}) %>% rbindlist()

diff2[,p.adj:=p.adjust(p.value)]
