library(devtools)
library(parallel)
OUT.DIR <- '/home/ob219/share/as_basis/GWAS/por_calibrations/'
load_all("~/git/cupcake")
nsteps = 1e+4

## here we wish to see what happens when we alter the various parameters

## first what happens with constant target.or and prob but samples size ?


sample.sizes <- c(2500,5000,10000)
target.or <- 2
target.prob <- 0.01
results <- list()
results.ss <- lapply(sample.sizes, function(ss){
  message(ss)
  fname <- sprintf("por_%d_%.1f_%.2f.RDS",ss,target.or,target.prob) %>% file.path(OUT.DIR,.)
  af <- seq(0.01,0.99,length.out=999)
  sprintf("ss_%d",ss)
  res <- mclapply(af,lor_f,n=ss,target.or=target.or,target.prob=target.prob,n.steps=nsteps,mc.cores=8)
  results <- data.table(cbind(sample.size=ss,target.or=target.or,target.prob=target.prob,f=af,do.call("rbind",res)))
  results
}) %>% rbindlist


library(cowplot)
m.ss <- melt(results.ss,id.vars=c('f','sample.size','target.or','target.prob'))
m.ss[,group:=paste(variable,sample.size,sep=' ')]
m.ss[,sample.size:=factor(sample.size,levels=sample.sizes)]



## now what happens if we fix sample size and alter target.or

ss <- 2500
target.ors <- c(2,3,4)
target.prob <- 0.01
results.tor <- lapply(target.ors, function(tor){
  message(tor)
  af <- seq(0.01,0.99,length.out=999)
  res <- mclapply(af,lor_f,n=ss,target.or=tor,target.prob=target.prob,n.steps=nsteps,mc.cores=8)
  results <- data.table(cbind(sample.size=ss,target.or=tor,target.prob=target.prob,f=af,do.call("rbind",res)))
  results
}) %>% rbindlist


m.tor <- melt(results.tor,id.vars=c('f','sample.size','target.or','target.prob'))
m.tor[,group:=paste(variable,target.or,sep=' ')]
m.tor[,target.or:=factor(target.or,levels=c(2,3,4))]







## now what happens if we fix sample size,target.or and change target.prob

ss <- 2500
tor <- 2
target.probs <- c(0.01,0.05,0.1)
results.prob <- lapply(target.probs, function(prob){
  message(tor)
  af <- seq(0.01,0.99,length.out=999)
  res <- mclapply(af,lor_f,n=ss,target.or=tor,target.prob=prob,n.steps=nsteps,mc.cores=8)
  results <- data.table(cbind(sample.size=ss,target.or=tor,target.prob=prob,f=af,do.call("rbind",res)))
  results
}) %>% rbindlist


m.prob <- melt(results.prob,id.vars=c('f','sample.size','target.or','target.prob'))
m.prob[,group:=paste(variable,target.prob,sep=' ')]
m.prob[,target.prob:=factor(target.prob,levels=c(0.01,0.05,0.1))]

ppss <- ggplot(m.ss,aes(x=f,y=value,col=sample.size,lty=variable)) + geom_line(size=1) +
labs(x="Allele Frequency",y="Posterior log(OR)") + background_grid(major = "xy", minor = "none") +
scale_linetype_manual(name="GT",labels=c("0/0","0/1 or 1/0","1/1"),
values=c("11"=6,"01"=1,"00"=3)) + scale_color_discrete(name="N") #+ theme(legend.position="top")

pptor <- ggplot(m.tor,aes(x=f,y=value,lty=variable,col=target.or)) + geom_line(size=1) +
labs(x="Allele Frequency",y="Posterior log(OR)") + background_grid(major = "xy", minor = "none") +
scale_linetype_manual(name="GT",labels=c("0/0","0/1 or 1/0","1/1"),
values=c("11"=6,"01"=1,"00"=3)) + scale_color_discrete(name="Targ. OR") #+ theme(legend.position="top")

ppprob <- ggplot(m.prob,aes(x=f,y=value,lty=variable,col=target.prob)) + geom_line(size=1) +
labs(x="Allele Frequency",y="Posterior log(OR)") + background_grid(major = "xy", minor = "none") +
scale_linetype_manual(name="GT",labels=c("0/0","0/1 or 1/0","1/1"),
values=c("11"=6,"01"=1,"00"=3)) + scale_color_discrete(name="Targ. Prob.") #+ theme(legend.position="top")
save_plot(plot_grid(ppss,ppprob,labels=c('a','b')),file="~/tmp/por_calli_curves.pdf",base_width=11)




if(FALSE){
  a0 <- 499.900000
  b0 <- 4499.100000
  a1 <- 16.659065
  b1 <- 146.083910
  p <- seq(0.005,0.2,length.out=9999)
  pdf("~/tmp/explain_por.pdf")
  plot(p,dbeta(p,shape1=a0,shape2=b0),type='l',xlab="MAF",ylab="Density", main="MAF 0.1 P(OR>2)=0.01")
  lines(p,dbeta(p,shape1=a1,shape2=b1),col='red')
  ## add expected values
  abline(v=a0/(a0+b0),lty=2)
  abline(v=a1/(a1+b1),col='red',lty=2)
  dev.off()
}
