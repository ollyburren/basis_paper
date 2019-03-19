library("cowplot")
library("ggrepel")

## plot both scree and biplot

sbi <- function(file){
  pc.emp <- readRDS(file)
  vexp <- summary(pc.emp)[['importance']][2,]
  PC1.var<-signif(vexp["PC1"]*100,digits=3)
  PC2.var<-signif(vexp["PC2"]*100,digits=3)
  M <- cbind(as.data.table(pc.emp$x),trait=rownames(pc.emp$x))
  scp <- cbind(data.table(vexp),pcs=factor(names(vexp),levels=names(vexp)))
  scp[,cs:=cumsum(vexp)]
  scp$group=1
  text.size <- 20
  ## do a scree plot
  ppl <- ggplot(scp,aes(x=pcs,y=vexp,group=group)) + geom_point() + geom_line() + ylab("Variance Explained") + xlab("Principal Components") +
  theme(axis.text=element_text(size=text.size),axis.title=element_text(size=text.size),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ppr<-ggplot(M,aes(x=PC1,y=PC2,label=trait)) + geom_point(size=3) + geom_text_repel(size=7) + # hjust = 0, nudge_x = 0.005)  +
  scale_color_discrete(guide=FALSE) + scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) +
  xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +  background_grid(major = "xy", minor = "none") +
  theme(axis.text=element_text(size=text.size),axis.title=element_text(size=text.size))
  plot_grid(ppl, ppr, labels = "auto",label_size=text.size)
}


## code to plot the beta basis
BASIS_BETA_FILE <- '/home/ob219/share/as_basis/GWAS/support/basis_beta_gwas.RDS'
BASIS_NOSHRINK_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_noshrink_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'

beta <- sbi(BASIS_BETA_FILE)
save_plot("~/tmp/basis_beta_plot.pdf",beta,ncol = 2,base_height=5)

gamma <- sbi(BASIS_NOSHRINK_FILE)
save_plot("~/tmp/basis_gamma_plot.pdf",gamma,ncol = 2,base_height=5)

shrunk <- sbi(BASIS_FILE)
save_plot("~/tmp/basis_shrunk_plot.pdf",shrunk,ncol = 2,base_height=5)
