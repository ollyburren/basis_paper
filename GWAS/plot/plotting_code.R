library(ggplot2)
library(cowplot)
library(ggrepel)
library(data.table)
library(magrittr)

BASIS_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/basis_june10k.RDS'
VARIANCE_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/analytical_variances_june10k.RDS'
SUMMARY_STATS_PROJ_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/summary_stat_proj_june10k.RDS'
IND_PROJ_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/ind_proj_june10k.RDS'
IND_PROJ_DATA_DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/individual_proj'
MANIFEST <- '/home/ob219/git/as_basis/manifest/as_manifest_july.tsv'

## comparison plot of jia individual projections and


if(!file.exists(IND_PROJ_FILE)){
  fs <- list.files(path=IND_PROJ_DATA_DIR,pattern="*.RDS",full.names=TRUE)
  all.ind.res <- lapply(fs,function(f){
    cm <- readRDS(f)
    data.table(trait=basename(f) %>% gsub("_projection.RDS","",.),
          individual=rownames(cm),cm)
  }) %>% rbindlist
  saveRDS(all.ind.res,file=IND_PROJ_FILE)
}

basis <- readRDS(BASIS_FILE)
proj.summary <- readRDS(SUMMARY_STATS_PROJ_FILE)
proj.individual <- readRDS(IND_PROJ_FILE)


## create a custom colour palette using rcolorbrewer

cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99')
names(cols) <- rownames(basis$x)

## plot basis traits for orientation

proj.basis <- data.table(trait=rownames(basis$x),individual='BASIS',basis$x)

all<-rbindlist(list(proj.basis,proj.summary,proj.individual))
## for biobank we want certain traits to have the same color

ml<-list(
  CD = 'bb_CD',
  CEL = 'bb_CEL',
  MS = 'bb_MS',
  RA = 'bb_RA',
  SLE = 'bb_SLE',
  T1D = 'bb_T1D',
  UC = 'bb_UC',
  PBC = 'PBC',
  PSC = 'PSC',
  asthma = 'bb_asthma'
)
all[,cat:='Other']
for(i in seq_along(ml)) {
    all[trait %in% c(names(ml)[i], ml[i]), cat:=names(ml)[i]]
}

## ok next add in confidence intervals for everything


dp <- function(b,DT,c1='PC1',c2='PC2',error_bar=TRUE,coords=FALSE,labels){
  ## get variance explained and generate labels
  p1.var <- signif(summary(b)[['importance']][2,][c1]*100,digits=3) %>% sprintf("%s (%.1f%%)",c1,.)
  p2.var <- signif(summary(b)[['importance']][2,][c2]*100,digits=3) %>% sprintf("%s (%.1f%%)",c2,.)
  keep.cols <- c('trait','cat',c1,c2,paste(c(c1,c2),'ci',sep='_'))
  pc<-dcast(DT,trait+cat~variable)[,keep.cols,with=FALSE]
  setnames(pc,c('t','cat','pcx','pcy','cix','ciy'))
  pc[t %in% rownames(b$x),c('cix','ciy'):=list(NA,NA)]
  y.int <- pc[t=='control',]$pcx
  x.int <- pc[t=='control',]$pcy
  if(!missing(labels)){
      pc[!t %in% labels,t:='']
  }else{
    pc[t=='control',t:='']
  }
  pc[,t:=gsub("_"," ",t)]
  ggp <- ggplot(pc,aes(x=pcx,y=pcy,color=cat,label=t,
    ymin=pcy-ciy,ymax=pcy+ciy,xmin=pcx-cix,xmax=pcx+cix,alpha=t!='')) +
  geom_point(size=1) + geom_text_repel(box.padding=0.5)  + scale_color_discrete(guide=FALSE) +
  xlab(p1.var) + ylab(p2.var) +
  geom_hline(yintercept = x.int,col='black',size=0.5,lty=2,alpha=0.5) +
  geom_vline(xintercept = y.int,col='black',size=0.5,lty=2,alpha=0.5) +
  scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) +
  background_grid(major = "xy", minor = "none")
  if(coords){
    #xlims<-pc[,list(xmin=min(pcx-cix,na.rm=TRUE),xmax=max(pcx + cix,na.rm=TRUE))] %>% t %>% as.vector
    xlims<-pc[!is.na(cix),list(xmin=min(pcx,na.rm=TRUE),xmax=max(pcx,na.rm=TRUE))] %>% t %>% as.vector
    #ylims<-pc[,list(xmin=min(pcy-ciy,na.rm=TRUE),xmax=max(pcy + ciy,na.rm=TRUE))] %>% t %>% as.vector
    ylims<-pc[!is.na(ciy),list(xmin=min(pcy,na.rm=TRUE),xmax=max(pcy,na.rm=TRUE))] %>% t %>% as.vector
    ggp <- ggp + coord_cartesian(xlim=xlims,ylim=ylims,expand=TRUE)
  }
  if(error_bar)
    ggp <- ggp + geom_errorbar(alpha=0.3) + geom_errorbarh(alpha=0.3)
  ggp
}

var <- readRDS(VARIANCE_FILE)
man <- fread(MANIFEST)[,.(trait,cases,controls)]
man <- rbind(man,data.table(trait='control',cases=0,controls=0))
M <- melt(all,id.vars=c('trait','individual','cat'),measure.vars=paste0('PC',1:11))


M <- merge(M,var,by.x='variable',by.y='pc')
M <- merge(M,man,by.x='trait',by.y='trait')


compute_ci <- function(m,case,ctrl){
  ## do using logs to avoid overflow
  f <- exp(log(case+ctrl) - (log(case) + log(ctrl)))
  v <- f * m
  sqrt(v) * 1.96
}
Mt <- copy(M)
Mt[,c('variable','value'):=list(paste(variable,'ci',sep='_'),compute_ci(mfactor,cases,controls))]
M <- rbind(M,Mt)

dat <- M[individual=='BASIS' | grepl("^bb",trait) | trait=='control',]
## biobank plot
kl <- c(rownames(basis$x),unlist(ml,use.names=FALSE)) %>% unique
kl <- c(kl,'bb_hyperthyroidism_thyrotoxicosis','bb_hypothyroidism_myxoedema','bb_psoriatic_arthropathy','bb_colitis','bb_AS')
bb_plot <- dp(basis,dat,labels=kl)
save_plot("~/tmp/bb_plot.pdf",bb_plot,base_height=6)

## JIA summary stats
dat <- M[individual=='BASIS' | grepl("^jia",trait,ignore.case=TRUE) | trait=='control',]
dat[,trait:=gsub("^jia_","",trait)]
jia_sum_plot <- dp(basis,dat,c2='PC3')

## JIA individual stats

M <- melt(all,id.vars=c('trait','individual','cat'),measure.vars=paste0('PC',1:11))[!is.na(individual),]
Mt <- M[,list(individual=head(individual,n=1),value=sd(value) * 1.96),by=c('trait','cat','variable')]
Mt[,variable:=paste(variable,'ci',sep='_')]
## take the mean for plotting
M <- M[,list(individual=head(individual,n=1),value=mean(value)),by=c('trait','cat','variable')]
M <- rbind(M,Mt)

dat <- M[individual=='BASIS' | grepl("^jia",trait,ignore.case=TRUE) | trait=='control',]
dat[,trait:=gsub("^jia_","",trait)]

## work out coords for plotting
dat[,trait:=gsub("^jia","",trait)]
jia_ind_plot <- dp(basis,dat,c2='PC3',coords=TRUE,error_bar=FALSE)
comp_plot <- plot_grid(jia_sum_plot,jia_ind_plot,labels=c('A','B'))
save_plot("~/tmp/ind_vs_summ_jia.pdf",comp_plot,base_width=10)

## TODO add in the hclust stuff 
