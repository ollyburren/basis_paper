## read in LD hub and make a nice heatmap

library(xlsx)

dat.rg <- read.xlsx(file="/home/ob219/tmp/LD-Hub_genetic_correlation_221x221_no_ENIGMA.xlsx",2)
dat.rg<-data.table(dat.rg)
setnames(dat.rg,"NA.",'study')


M <- melt(dat.rg,id.vars='study')

M <- M[value!='/']
M[,rG:=as.numeric(value)]
M <- M[!is.na(rG)]

keep <- c('RA'='RA_GWASmeta_European_v2.txt.gz.uniform.txt.noMHC.sumstats_deGC.gz','PBC'='PASS_Primary_biliary_cirrhosis.noMHC.sumstats.gz',
SLE='PASS_Lupus.noMHC.sumstats.gz',AST='gabriel_results.txt.noMHC.sumstats.gz',UC='EUR.UC.gwas.assoc.gz.noMHC.sumstats.gz',CD='EUR.CD.gwas.assoc.gz.noMHC.sumstats.gz')

MK <- M[(study %in% keep) & (variable %in% keep),]

dat.p <- read.xlsx(file="/home/ob219/tmp/LD-Hub_genetic_correlation_221x221_no_ENIGMA.xlsx",3)
dat.p<-data.table(dat.p)
#setnames(dat.p,"NA.",'study')
M <- melt(dat.p,id.vars='study')
M[,p:=as.numeric(value)]
M <- M[!is.na(p),]
M[,Z:=]

keep <- c('RA_GWASmeta_European_v2.txt.gz.uniform.txt.noMHC.sumstats_deGC.gz','PASS_Primary_biliary_cirrhosis.noMHC.sumstats.gz',
'PASS_Lupus.noMHC.sumstats.gz','gabriel_results.txt.noMHC.sumstats.gz','EUR.UC.gwas.assoc.gz.noMHC.sumstats.gz','EUR.CD.gwas.assoc.gz.noMHC.sumstats.gz')

## try all results as this has the info I need
dat.all <- read.xlsx(file="/home/ob219/tmp/LD-Hub_genetic_correlation_221x221_no_ENIGMA.xlsx",1)
dat.all <- data.table(dat.all)
setnames(dat.all,"NA.",'study')
M <- melt(dat.all,id.vars='study')
M <- M[value!=1,]
M[,value:=gsub("^[ ]*","",value)]
M[,c('study1','study2','rg','se','z','p','h2_obs','h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se'):=tstrsplit(value,"[ ]+")]
MK <- M[(study %in% keep) & (variable %in% keep),]
MK[,rg:=as.numeric(rg)]
MK[,z:=as.numeric(z)]
MK[,p:=as.numeric(p)]

for(k in names(keep)){
  MK[study1==keep[k],study1:=k]
  MK[study2==keep[k],study2:=k]
}

sorder <- c("RA","PBC","SLE","AST","UC","CD") %>% rev

MK[,c('study1','study2'):=list(factor(study1,levels=sorder),factor(study2,levels=sorder))]
library(cowplot)
MK[rg==1,z:=0]
#MK[rg==1,rg:=NA]
MK[,p.adj:=p.adjust(p,method="bonferroni")]
MK[p.adj>0.05,z:=0]
pp1 <- ggplot(MK,aes(x=study1,y=study2,fill=z,label=signif(rg,digits=1))) + geom_tile(color="black") +
scale_fill_gradient2("Z",space="Lab") + geom_text() +
xlab("Trait") + ylab("Trait")
save_plot(pp1,file="~/tmp/ldhub_imd.pdf")
