library(devtools)
library(parallel)
install_github('ollyburren/cupcake')
library(cupcake)
OUT.DIR <- '/home/ob219/rds/hpc-work/as_basis/support/'
#load_all("~/git/cupcake")
target.or <- 2
target.prob <- 0.01
n.sample <- 2500 # number of controls
nsteps = 1e+4 # larger values make more accurate approximations to the integral
fname <- sprintf("por_%d_%.1f_%.2f.RDS",n.sample,target.or,target.prob) %>% file.path(OUT.DIR,.)
af <- seq(0.01,0.99,length.out=999)
results <- mclapply(af,lor_f,n=n.sample,target.or=target.or,target.prob=target.prob,n.steps=nsteps,mc.cores=6)
results <- data.table(cbind(f=af,do.call("rbind",results)))
saveRDS(results,file=fname)


if(FALSE){
  library(ggplot2)
  library(cowplot)
  colnames(results) <- make.unique(colnames(results))
  m <- melt(results,"f")
  ppf<-ggplot(m,aes(x=f,y=value,col=variable,group=variable)) +
    geom_point() +
    labs(x="Allele Frequency",y="Posterior log(OR)") + background_grid(major = "xy", minor = "none") +
    scale_color_manual(name="Genotype",labels=c("0/0","0/1 or 1/0","1/1"),
    values=c("11"='firebrick',"01"='black',"00"='dodgerblue')) +
    theme(legend.position=c(0.1,0.8))
  ppf
}

if(FALSE){
  ## can we plot the two beta distributions ?
  lor_f(0.1,n=2500,nsim=nsims,target.or=target.or,target.prob=target.prob,n.steps=nsteps)
  a0 <- 499.900000
  b0 <- 4499.100000
  a1 <- 16.659065
  b1 <- 146.083910
  p <- seq(0.005,0.2,length.out=9999)
  pdf("~/tmp/explain_por.pdf")
  plot(p,dbeta(p,shape1=a0,shape2=b0),type='l',xlab="MAF",ylab="Density", main="MAF 0.1 P(OR>2)=0.01")
  lines(p,dbeta(p,shape1=a1,shape2=b1),col='red')
  ## add expected values
  abline(v=a0/(a0+b0))
  abline(v=a1/(a1+b1),col='red')
  dev.off()
}
