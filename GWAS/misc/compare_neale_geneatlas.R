library(cowplot)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/03_10_13_traits_0919_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)

ga.srd <- res.DT[category=='geneatlas_srd',.(trait=gsub("GA:","",trait),pc=variable,ga.Z=Z,ga.cases=n1,ga.controls=n0,ga.delta=delta)]
neale.srd <- res.DT[category=='bb_disease',.(trait=gsub("bb_SRD:","",trait),pc=variable,neale.Z=Z,neale.cases=n1,neale.controls=n0,neale.delta=delta)]

M <- merge(ga.srd,neale.srd,by=c('trait','pc'))

M[,pc:=factor(pc,levels=paste0('PC',1:14))]

plot(M$neale.delta,M$ga.delta)
abline(a=0,b=1,col='red',lty=2)


M[,label:='NONE']
M[abs(neale.Z)>2.8 & abs(ga.Z)<2.8,label:='NEALE']
M[abs(neale.Z)<2.8 & abs(ga.Z)>2.8,label:='GA']
M[abs(neale.Z)>2.8 & abs(ga.Z)>2.8,label:='BOTH']

## compare cases

ggplot(M,aes(x=neale.cases,y=ga.cases)) + geom_point() + geom_abline(b=1,col='red',lty=2)


ggplot(M,aes(x=neale.Z,y=ga.Z,col=label)) +
geom_point() + geom_abline(a=0,b=1,col='red',lty=2) +
theme_bw() + facet_wrap(~pc,scales='free')

qqnorm(M$neale.Z)
