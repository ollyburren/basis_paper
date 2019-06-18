## plot JIA results for thesis


library(cowplot)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/09_04_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)



jiat <- res.DT[(category=='wong_aav' & grepl("meta",trait)) | category=='lyons_egpa',]
jiat[category=='lyons_egpa',trait:=paste('egpa',trait,sep=':')]
jiat[,c('ci.lo','ci.hi'):=list(delta-(sqrt(variance) * 1.96),delta+(sqrt(variance) * 1.96)),]
jiat[,variable:=factor(variable,levels=paste0('PC',1:11))]
jiat[,trait:=sub("^jia_","",trait)]

pd <- position_dodge(0.2)
pa <- ggplot(jiat,aes(x=variable,y=delta,group=trait,col=trait)) + geom_point(position=pd) +
geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi),alpha=0.3,position=pd) + geom_line(position=pd) +
ylab(expression(Delta*~"control score")) + geom_hline(yintercept=0,lty=2,col='black') + xlab('Principal component') +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + scale_colour_discrete("Disease Subtype")
pa

## plot geom_tile

bb <- jiat[,.(variable,trait,delta,p.value)]
bb[,p.adj:=p.adjust(p.value,method="bonferroni"),by='variable']
## order the traits by the clustering
bb[,trait:=factor(trait,levels=c('RFpos','ERA','sys','RFneg','EO','PO','PsA'))]
bb[,pc:=factor(variable,levels=paste0('PC',1:11))]
bb[,is.signif:='']
bb[p.adj<0.05,is.signif:='*']
#bb[,Z.plot:=0]
#bb[p.adj<0.05,Z.plot:=Z]

bbplot <- ggplot(bb,aes(x=pc,y=trait,fill=delta,label=is.signif))  +
geom_tile(color='black') + geom_text(cex=10) + scale_fill_gradient2("Difference\nfrom control") + xlab("Principal Component") +
ylab("Trait") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
save_plot(bbplot,file="~/tmp/jia_heatmap.pdf",base_width=6)

## what do signatures for PsA vs Psoriasis look like ?

jiat <- res.DT[trait %in% c('bb_psoriasis','jia_PsA','bb_rheumatoid.arthritis','jia_RFpos','bb_ankylosing.spondylitis','jia_ERA'),]
jiat[,variable:=factor(variable,levels=paste0('PC',1:11))]
jiat[,trait:=sub("^jia_","",trait)]
jiat[trait=='bb_psoriasis',trait:='Psoriasis']
jiat[trait=='bb_rheumatoid.arthritis',trait:='RA']
jiat[trait=='bb_ankylosing.spondylitis',trait:='AS']

bb <- jiat[,.(variable,trait,delta,p.adj)]
#bb[,p.adj:=p.adjust(p.value,method="bonferroni"),by='variable']
## order the traits by the clustering
bb[,trait:=factor(trait,levels=c('PsA','Psoriasis','RFpos','RA','ERA','AS'))]
bb[,pc:=factor(variable,levels=paste0('PC',1:11))]
bb[,is.signif:='']
bb[p.adj<0.05,is.signif:='*']
bb[,extreme:=sign(delta)]
bbplot2 <- ggplot(bb,aes(x=pc,y=trait,fill=extreme,label=is.signif))  +
geom_tile(color='black') + geom_text(cex=10) + scale_fill_gradient2("Difference\nfrom control") + xlab("Principal Component") +
ylab("Trait") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

jiat <- res.DT[trait %in% c('bb_psoriasis','bb_rheumatoid.arthritis','bb_ankylosing.spondylitis') | category=='bowes_jia',]
jiat[,extreme:=sign(delta)]

M <- melt(jiat[,.(pc=variable,trait,extreme)],id.vars=c('trait','pc'),measure.vars=c('extreme')) %>% dcast(trait~pc)

pcs <- c('PC4','PC6')

all.m <- lapply(1:nrow(M),function(i){
  lapply(1:3,function(j){
      tp <- paste(M$trait[i],M$trait[j],sep='-')
      message(tp)
      if(M$trait[i] < M$trait[j])
        return(data.table(tp=tp,match.count=-1))
      mcount <- ((rbind(M[i,pcs,with=FALSE],M[j,pcs,with=FALSE]) %>% as.matrix %>% colSums %>% abs)  == 2) %>% sum
      data.table(tp=tp,match.count=mcount)
  }) %>% rbindlist
}) %>% rbindlist

for(i in nrow(M)){

}

bZ <- dcast(forhc[variable=='Z'],Subtype~pc)
bZm <- bZ[,-1] %>% as.matrix
rownames(bZm) <- bZ$Subtype
hc <- dist(bZm) %>% hclust
pred.DT[,Subtype:=factor(Subtype,levels=hc$labels[hc$order])]
pred.DT[,plotZ:=0]
pred.DT[p.adj<0.05,plotZ:=Z]
pc <- ggplot(pred.DT,aes(x=variable,y=Subtype,fill=plotZ,label=signif(value,digits=1))) +
geom_tile(color='black') + geom_text(cex=2.1) + scale_fill_gradient2("Z score") + xlab("Principal Component") +
ylab("Trait") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


## do a line plot comparing PsA with UKBB PSO and PsA


plot.DT <- res.DT[trait %in% c('bb_psoriasis','bb_psoriatic.arthropathy','jia_PsA'),]
plot.DT[,c('ci.lo','ci.hi'):=list(delta-(sqrt(variance) * 1.96),delta+(sqrt(variance) * 1.96)),]
plot.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
plot.DT[,trait:=sub("^jia_","",trait)]

pd <- position_dodge(0.2)
pa <- ggplot(plot.DT,aes(x=variable,y=delta,group=trait,col=trait)) + geom_point(position=pd) +
geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi),alpha=0.3,position=pd) + geom_line(position=pd) +
ylab(expression(Delta*~"control score")) + geom_hline(yintercept=0,lty=2,col='black') + xlab('Principal component') +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + scale_colour_discrete("JIA Subtype")

plot.DT <- res.DT[trait %in% c('bb_ankylosing.spondylitis','jia_ERA'),]
plot.DT[,c('ci.lo','ci.hi'):=list(delta-(sqrt(variance) * 1.96),delta+(sqrt(variance) * 1.96)),]
plot.DT[,variable:=factor(variable,levels=paste0('PC',1:11))]
plot.DT[,trait:=sub("^jia_","",trait)]

pd <- position_dodge(0.2)
pa <- ggplot(plot.DT,aes(x=variable,y=delta,group=trait,col=trait)) + geom_point(position=pd) +
geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi),alpha=0.3,position=pd) + geom_line(position=pd) +
ylab(expression(Delta*~"control score")) + geom_hline(yintercept=0,lty=2,col='black') + xlab('Principal component') +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + scale_colour_discrete("JIA Subtype")
