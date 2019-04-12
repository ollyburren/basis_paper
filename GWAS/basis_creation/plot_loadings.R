library(cowplot)

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)


## let's plot all the loadings

load.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)

M <- melt(load.DT,id.vars='pid')
M[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
M <- M[order(chr),as.numeric]

pos.DT<-M[,.(pid,chr,pos)] %>% unique
pos.DT <- pos.DT[order(chr,pos)]
pos.DT <- pos.DT[,list(pid=pid,minpos=pos-min(pos)+1),by=chr]
pos.DT[,cpos:=cumsum(minpos)]
M <- merge(M[,.(pid,variable,value)],pos.DT[,.(pid,chr,cpos)],by='pid')
M[,evenchr:=chr%%2]
bl <- ggplot(M[abs(value)>0.0001 & variable!='PC11',],aes(x=cpos,y=value,color=(chr %% 2)!=0)) + geom_point(size=0.5) + facet_wrap(~variable,nrow=5) +
guides(color=FALSE) + ylab("Principal Component Loading") + xlab("Genome Position") + scale_colour_manual(values=c('TRUE'='firebrick1','FALSE'='dodgerblue'))
save_plot(bl,file="~/tmp/loadings.pdf",base_width=10,base_height=7)

## there may be some residual contamination with the distal end of MHC region ?
