library(cowplot)
data.dir <- '/home/ob219/share/as_basis/GWAS/individual_data/individual_proj/'
files <- list.files(path=data.dir,full.names=TRUE)
files <- files[grep("^jia",basename(files))]
dttypes <- gsub("jia([^\\_]+)\\_.*","\\1",basename(files))
names(files) <- dttypes
all.res <- lapply(seq_along(files),function(i){
  tmp <- readRDS(files[i]) %>% data.table
  tmp[,trait:=names(files)[i]]
}) %>% rbindlist
all.res <- all.res[!trait %in% c('missing','UnA'),]
all.res[,id:=1:.N]
m <- melt(all.res,id.vars=c('id','trait'))


BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
#basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
#basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
#basis.DT[,delta:=(value-control.loading)]
m <- merge(m,control.DT,by.x='variable',by.y='PC')
m[,mc:=value-control.loading]
ms <- m[,list(mean.score=mean(mc),sd.score=sd(mc)),by=c('trait','variable')]
ms[,ci:=1.96 * sd.score]
ms[,c('lower','upper'):=list(mean.score-ci,mean.score+ci)]

## plot pc3
plot.dat <- ms[variable=='PC3',]
plot.dat[,trait:=factor(trait,levels=plot.dat[order(mean.score,decreasing=FALSE),]$trait)]

pp<-ggplot(plot.dat,aes(x=trait,y=mean.score)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) +
xlab("JIA subtype") + ylab("Mean change in basis\nloading from control")
save_plot(pp,file="~/tmp/jia_individual_gt_pc3.pdf")

##


ggplot(m[variable %in% c('PC1','PC2','PC3','PC4'),],aes(y=value,x=trait)) + geom_violin() + facet_wrap(~variable,scales="free_x")
