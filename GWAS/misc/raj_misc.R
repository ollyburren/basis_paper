## how many SNPs do we have for Raj


all.lms.cd4 <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd4/regression_models.RDS")
names(all.lms.cd4)<-paste('PC',1:11,sep="")
all.lms.cd14 <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd14/regression_models.RDS")
names(all.lms.cd14)<-paste('PC',1:11,sep="")

for.reg.cd4 <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd4/cd4.RDS")
for.reg.cd14 <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd14/cd14.RDS")

## for each model we can obtain t-statistics for the beta coefficients as p-values.
## we can the check for inflation by plotting as a qqplot trellis plot.

pc.p.cd4 <- lapply(seq_along(names(all.lms.cd4)),function(i){
  pc<-names(all.lms.cd4)[i]
  all.p <- sapply(all.lms.cd4[[i]],function(x){
    summary(x)$coefficient["value",3]
  })
  DT <- data.table(PC=pc,probe=names(for.reg.cd4)[grep('^P\\_',names(for.reg.cd4))],p.coeff=all.p,cell.type='CD4+')
})


pc.p.cd14 <- lapply(seq_along(names(all.lms.cd14)),function(i){
  pc<-names(all.lms.cd14)[i]
  all.p <- sapply(all.lms.cd14[[i]],function(x){
    summary(x)$coefficient["value",3]
  })
  DT <- data.table(PC=pc,probe=names(for.reg.cd14)[grep('^P\\_',names(for.reg.cd14))],p.coeff=all.p,cell.type='CD14+')
})

## do qq stuff for both

saveRDS(list(cd4=pc.p.cd4,cd14=pc.p.cd14),"~/share/as_basis/GWAS/raj/all_lm.p.RDS")

for.qq.cd4<-rbindlist(lapply(pc.p.cd4,function(x){
  points <- qqnorm(x$p.coeff,plot.it=FALSE)
  data.table(pc=unique(x$PC),x1=points[[1]],y1=points[[2]],cell.type='CD4+')
}))



for.qq.cd14<-rbindlist(lapply(pc.p.cd14,function(x){
  points <- qqnorm(x$p.coeff,plot.it=FALSE)
  data.table(pc=unique(x$PC),x1=points[[1]],y1=points[[2]],cell.type='CD14+')
}))

for.qq <- rbind(for.qq.cd4,for.qq.cd14)

for.qq$pc <- factor(for.qq$pc,levels=paste0('PC',1:11))


library(cowplot)
gpl <- ggplot(for.qq,aes(x=x1,y=y1,color=cell.type)) +
geom_abline(intercept=0,slope=1,color='dodgerblue',lty=2) +
geom_point(size=0.1) + facet_wrap(~pc,ncol=3,nrow=4) +
ylab("Observed -log10(P)") + xlab("Expected -log10(P)") +
scale_color_discrete("Cell Type") + theme(legend.position="bottom") +
guides(colour = guide_legend(override.aes = list(size=5)))

save_plot(gpl,file="~/tmp/raj_qq.pdf",base_height=7)
save_plot(gpl,file="~/tmp/raj_qq.png",base_height=7)

## what about top genes

ffdr <- rbind(pc.p.cd14 %>% rbindlist,pc.p.cd4 %>% rbindlist)

ffdr[cell.type=='CD4+',df:=211]
ffdr[cell.type=='CD14+',df:=209]
ffdr[,p:=pt(abs(p.coeff),lower.tail=FALSE,df=df) * 2]

ffdr[,p.adj:=p.adjust(p,method="BH"),by=c('PC','cell.type')]
ffdr[p.adj<0.05,]

## pull out the relevant models and coeff

pids <- names(for.reg.cd14)[grep('^P\\_',names(for.reg.cd14))]
idx.pc6 <- which(pids == ffdr[p.adj<0.05 & PC=='PC6',]$probe)
idx.pc10 <- which(pids == ffdr[p.adj<0.05 & PC=='PC10',]$probe)

library(xtable)
rbind(summary(all.lms.cd14[[6]][[idx.pc6]])$coefficient[2,],
summary(all.lms.cd14[[10]][[idx.pc10]])$coefficient[2,]) %>% xtable


## here look at the expected loading and see if it is different from controls
files <- list.files(path="/home/ob219/share/as_basis/GWAS/individual_data/individual_proj/",pattern="*.RDS",full.names=TRUE)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])

## look at bootstrapping jia

jia.f <- files[grep("jia",files)]
library(parallel)
hg.bs <- mclapply(jia.f,function(f){
  trait <- basename(f)
  dat <- readRDS(f)
  abs <- lapply(1:1000,function(i){
    bs <- dat[sample.int(nrow(dat),size=211,replace=TRUE),]
    lapply(paste0('PC',1:11),function(pc){
      data.table(pc=pc,av=mean(bs[,pc]),sd=sd(bs[,pc],trait)
      #data.table(pc=pc,t=t.test(bs[,pc],mu=control.DT[PC==pc,]$control.loading)$statistic,trait=trait)
    }) %>% rbindlist
  }) %>% rbindlist
  abs[,list(mean.t=mean(t),sd.t=sd(t)),by=c('pc','trait')]
},mc.cores=8) %>% rbindlist

sd.factor <- qt(0.05/2,lower.tail=FALSE,df=211)

hg.bs[,c('ci.low','ci.high'):=list(mean.t-(sd.factor * sd.t),mean.t+(sd.factor * sd.t))]
hg.bs <- hg.bs[grep("missing|UnA",invert=TRUE,trait),]


#hg <- lapply(files[!files %in% jia.f],function(f){
hg <- lapply(files,function(f){
  trait <- basename(f)
  dat <- readRDS(f)
  lapply(paste0('PC',1:11),function(pc){
    csc <- dat[,pc] - control.DT[PC==pc,]$control.loading
    data.table(pc=pc,av=mean(csc),sdc=sd(csc),trait)
    #data.table(pc=pc,p=t.test(dat[,pc],mu=control.DT[PC==pc,]$control.loading)$p.value,trait=trait)
  }) %>% rbindlist
}) %>% rbindlist
hg[,study:='Raj']
hg[grep("^jia",trait),study:='JIA']
hg[,trait:=gsub("jia","",trait) %>% sub("_projection.RDS","",.)]
hg <- hg[!trait %in% c('missing','UnA')]
#hg[,p.adj:=p.adjust(p,method="bonferroni"),by='pc']
hg[,pc:=factor(pc,levels=paste0('PC',1:11))]
hg[,c('ci.low','ci.high'):=list(av-(1.96 * sdc),av+(1.96 * sdc))]

library(latex2exp)

pp1 <- ggplot(hg,aes(x=pc,y=av,fill=trait,ymin=ci.low,ymax=ci.high)) +
geom_errorbar(position="dodge") +
geom_bar(stat="identity",position = "dodge") + #coord_cartesian(ylim=c(0,7)) +
#geom_hline(yintercept=-log10(0.05),lty=2,col='red') + #ylab(TeX('$\\log_{10}(p_{adjust})$')) +
ylab("t-statistic") +
xlab("Principal Component") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pp2 <- ggplot(hg[study!='JIA',],aes(x=pc,y=mlp.adj,fill=trait)) +
geom_bar(stat="identity",position = "dodge") + coord_cartesian(ylim=c(0,7)) +
geom_hline(yintercept=-log10(0.05),lty=2,col='red') + ylab(TeX('$\\log_{10}(p.adj)$')) +
xlab("Principal Component") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pp3 <- plot_grid(pp1,pp2,nrow=2,labels="auto")
save_plot(pp3,file="~/tmp/compare_grs.pdf",base_height=7)



## why are cd4 and cd14 different - genotypes should be the same ?

cd14.files <- list.files(path='/home/ob219/share/as_basis/GWAS/individual_data/filtered_gt/raj/cd14',pattern="^chr*",full.names=TRUE)
library(annotSnpStats)
cd14.cs <- lapply(cd14.files,function(f){
  cs <- readRDS(f) %>% col.summary
  data.table(rs=rownames(cs),cs)
}) %>% rbindlist

cd4.files <- list.files(path='/home/ob219/share/as_basis/GWAS/individual_data/filtered_gt/raj/cd4',pattern="^chr*",full.names=TRUE)
cd4.cs <- lapply(cd4.files,function(f){
  cs <- readRDS(f) %>% col.summary
  data.table(rs=rownames(cs),cs)
}) %>% rbindlist

M<-merge(cd4.cs[,.(rs,cd4.maf=MAF)],cd14.cs[,.(rs,cd14.maf=MAF)],by='rs')

samp.cd14 <- readRDS(cd14.files[[1]]) %>% samples
samp.cd4 <- readRDS(cd4.files[[1]]) %>% samples
