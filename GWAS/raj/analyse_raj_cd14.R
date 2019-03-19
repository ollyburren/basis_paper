library(cowplot)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'
dat <- readRDS("/home/ob219/share/as_basis/GWAS/individual_data/individual_proj/raj_cd14_projection.RDS")
t.DT <- data.table(ind=rownames(dat),dat)
res <- melt(t.DT,id.vars='ind')
pc.emp <- readRDS(BASIS_FILE)
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x)
tmp <- basis.DT[trait=='control',] %>% t
ctrl.DT <- data.table(variable=rownames(tmp)[-1],control.loading=as.numeric(tmp[-1,1]))


res <- merge(res,ctrl.DT,by='variable')
res[,delta:=value-control.loading]
res[,variable:=factor(variable,levels=paste('PC',1:11,sep=''))]
ggplot(res,aes(x=variable,y=delta,group=ind)) + geom_line(alpha=0.3)

(load("/home/ob219/share/Projects/twas/raj-cd14-expression.RData"))
expr.DT<-data.table(t(expr))
setnames(expr.DT,paste('P',rownames(expr),sep='_'))
expr.DT[,sample:=colnames(expr)]
## load the transform matrix
translate <- readRDS('/home/ob219/share/Projects/twas/model_output/trans_raj-cd14.rds')
expr.DT[,sample_id:=translate[sample]]
for.reg <- merge(res,expr.DT,by.x='ind',by.y='sample_id')
setkey(expr.DT,sample_id)
setkey(res,ind)

ggplot(for.reg[variable=='PC1',],aes(x=value,y=P_7896908)) + geom_point() + geom_smooth(method='lm')
saveRDS(for.reg,file="/home/ob219/share/as_basis/GWAS/raj/cd14/cd14.RDS")
for.reg <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd14/cd14.RDS")

probes <- names(for.reg)[grep("P\\_[0-9]",names(for.reg))]

## first we attempt to fit a set of linear models where PC is the explantory model
library(parallel)
all.lms <- mclapply(paste('PC',1:11,sep=""),function(pc){
  message(pc)
  dat <- for.reg[variable==pc,]
  lapply(probes,function(p){
    #dat <- for.reg[variable==p,(pc,p,with=FALSE]
    #setnames(dat,c('PC','expression'))
    lm(sprintf("%s~value",p),dat)
  })
},mc.cores=8)
## large amount of reg
saveRDS(all.lms,"/home/ob219/share/as_basis/GWAS/raj/cd14/regression_models.RDS")
all.lms <- readRDS("/home/ob219/share/as_basis/GWAS/raj/cd14/regression_models.RDS")

names(all.lms)<-paste('PC',1:11,sep="")
## for each model we can obtain t-statistics for the beta coefficients as p-values.
## we can the check for inflation by plotting as a qqplot trellis plot.

## should create a library for this.
get_qq_dt <- function(x,l=0.99,minx=TRUE){
  n <- length(x)
  ## expected
  q <- -log10((n:1)/(n+1))
  ## observed
  x <- sort(x)
  ## l is a parameter that allows us to only plot a subset of the data
  n1 <- round(l*n)
  if(minx)
    return(data.table(expected=c(0,q[n1:n]),observed=c(0,x[n1:n])))
  return(data.table(expected=q,observed=x))
}

pc.p <- lapply(seq_along(names(all.lms)),function(i){
  pc<-names(all.lms)[i]
  all.p <- sapply(all.lms[[i]],function(x){
    summary(x)$coefficient["value",3]
  })
  DT <- data.table(PC=pc,probe=names(for.reg)[grep('^P\\_',names(for.reg))],p.coeff=all.p)
})


pc3.tstat<-pc.p[[3]]
library(locfdr)
pc3.tstat[,fdr:=locfdr(p.coeff)]

## for qq.plot

for.qq<-rbindlist(lapply(pc.p,function(x){
  points <- qqnorm(x$p.coeff,plot.it=FALSE)
  data.table(pc=unique(x$PC),x1=points[[1]],y1=points[[2]])
}))

for.qq$pc <- factor(for.qq$pc,levels=paste0('PC',1:11))

#for.qq <- rbindlist(lapply(pc.p,function(x) cbind(unique(x$PC),get_qq_dt(x$p.coeff))))
#for.qq$V1<-factor(for.qq$V1,levels=paste0('PC',1:11))

library(cowplot)
gpl <- ggplot(for.qq,aes(x=x1,y=y1)) +
geom_abline(intercept=0,slope=1,color='dodgerblue',lty=2) +
geom_point(size=0.5) + facet_wrap(~pc,ncol=3,nrow=4) +
ylab("Observed -log10(P)") + xlab("Expected -log10(P)") +
scale_color_manual(guide=FALSE,values=c('TRUE'='firebrick','FALSE'='black'))

## do fdr ?

pc.ps <- rbindlist(pc.p)
pc.ps[,p.val:=pt(abs(p.coeff),df=209,lower.tail=FALSE) * 2]
pc.ps[,adj.p:=p.adjust(p.val,method="BH"),by='PC']
pc.ps[adj.p<0.05,]


############STOP BELOW HERE IS OLD !!!!!




## first we attempt to fit a set of linear models where PC is the explantory model

## just do PC3 at first
dat <-  for.reg[variable=='PC3',]
pc3.reg <- mclapply(probes,function(p){
  #dat <- for.reg[variable==p,(pc,p,with=FALSE]
  #setnames(dat,c('PC','expression'))
  lm(sprintf("%s~value",p),dat)
},mc.cores=8)
pc3.p <- sapply(pc3.reg,function(x){
  summary(x)$coefficient["value",3]
})
DT.pc3 <- data.table(PC='PC3',probe=names(for.reg)[grep('^P\\_',names(for.reg))],p.coeff=pc3.p)


library(parallel)
all.lms <- mclapply(paste('PC',1:11,sep=""),function(pc){
  dat <- for.reg[variable==pc,]
  lapply(probes,function(p){
    #dat <- for.reg[variable==p,(pc,p,with=FALSE]
    #setnames(dat,c('PC','expression'))
    lm(sprintf("%s~value",p),dat)
  })
},mc.cores=8)
## large amount of reg
saveRDS(all.lms,"/home/ob219/share/as_basis/GWAS/raj/cd4/regression_models.RDS")


names(all.lms)<-paste('PC',1:11,sep="")
## for each model we can obtain t-statistics for the beta coefficients as p-values.
## we can the check for inflation by plotting as a qqplot trellis plot.

## should create a library for this.
get_qq_dt <- function(x,l=0.99,minx=TRUE){
  n <- length(x)
  ## expected
  q <- -log10((n:1)/(n+1))
  ## observed
  x <- sort(x)
  ## l is a parameter that allows us to only plot a subset of the data
  n1 <- round(l*n)
  if(minx)
    return(data.table(expected=c(0,q[n1:n]),observed=c(0,x[n1:n])))
  return(data.table(expected=q,observed=x))
}

pc.p <- lapply(seq_along(names(all.lms)),function(i){
  pc<-names(all.lms)[i]
  all.p <- sapply(all.lms[[i]],function(x){
    summary(x)$coefficient["value",3]
  })
  DT <- data.table(PC=pc,probe=names(for.reg)[grep('^P\\_',names(for.reg))],p.coeff=all.p)
})


pc3.tstat<-pc.p[[3]]
library(locfdr)
pc3.tstat[,fdr:=locfdr(p.coeff)]

## for qq.plot

for.qq <- rbindlist(lapply(pc.p,function(x) cbind(unique(x$PC),get_qq_dt(-log10(x$p.coeff)))))
for.qq$V1<-factor(for.qq$V1,levels=paste0('PC',1:10))

library(cowplot)
gpl <- ggplot(for.qq,aes(x=expected,y=observed)) +
geom_abline(intercept=0,slope=1,color='dodgerblue',lty=2) +
geom_point(size=0.5) + facet_wrap(~V1,ncol=3,nrow=4) +
ylab("Observed -log10(P)") + xlab("Expected -log10(P)") +
scale_color_manual(guide=FALSE,values=c('TRUE'='firebrick','FALSE'='black'))



(load("/home/ob219/share/Projects/twas/raj-cd14-expression.RData"))
expr.DT<-data.table(expr)
orig <- fread("zcat /home/ob219/share/Data/expr/Raj/GSE56034-CD14/GSE56034_GSM.ImmVarCD14.EU.PC20.txt.gz")
oval<-melt(orig,id.var='ID_REF')
oval[,ID_REF:=as.character(ID_REF)]
expr.DT[,ID_REF:=rownames(expr)]
exval <- melt(expr.DT,id.var='ID_REF')
setnames(exval,c('variable','value'),c('my.variable','my.value'))
M <- merge(oval,exval,by.x=c('ID_REF','variable'),by.y=c('ID_REF','my.variable'))
