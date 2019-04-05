library(cowplot)
library(scales)
library(devtools)
load_all("~/git/cupcake")
SHRINKAGE_METHOD<-'ws_emp_shrinkage'
OUTDIR <- '/home/ob219/share/as_basis/GWAS/variance_simulations_bootstrap_basis'
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'
BASIS_FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/basis_ws.RDS"
SHRINKAGE_FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/basis/shrink_ws.RDS"
REF_GT_DIR <- '/home/ob219/share/as_basis/GWAS/ctrl_gt/by.chr'
SNP_MANIFEST_FILE <- "/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/support/snp_support_bootstrap_USE.tab"
## compile callibrations
pc.emp <- readRDS(BASIS_FILE)
af <- list.files(path=OUTDIR,pattern="*.RDS",full.names=TRUE)
byexp<-split(af,gsub("\\_[a-z]+.RDS$","",basename(af)))

all.res <- lapply(names(byexp),function(x){
  all.res <- lapply(byexp[[x]],readRDS) %>% do.call('rbind',.)
  var <- apply(all.res,2,var)
  tmp<-strsplit(x,"\\_") %>% unlist
  data.table(total=tmp[1],ncases=tmp[2],variance=var,pc=names(var))
}) %>% rbindlist

CQF.lt <- qchisq(0.025, 500-1, lower.tail=FALSE)
CQF.ut <- qchisq(0.025, 500-1, lower.tail=TRUE)
## an approximate way to do this is variance * sqrt(2/(n-1)) ## I should ask Chris how this works
all.res[,c('ci.lower','ci.upper'):=list((variance * (500-1)/CQF.lt),((variance * (500-1))/CQF.ut)) ]
#library(scales)
all.res[,ci:=variance * sqrt(2/(500-1))]
all.res[,total:=factor(total,levels=sort(unique(as.numeric(total))))]
pd <- position_dodge(width=0.3)
## next load in the true variance that we use for basis
#all.vars <- readRDS(VARIANCE_FILE)
### create variances for the bootstrap thingy
w.DT <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
man.DT <- fread(SNP_MANIFEST_FILE)
shrink.DT <- readRDS(SHRINKAGE_FILE)
analytical.vars <- compute_proj_var(man.DT,w.DT,shrink.DT,REF_GT_DIR,SHRINKAGE_METHOD,quiet=FALSE)
all.vars <- data.table(pc=names(analytical.vars),mfactor=analytical.vars)
## make one for each of the configurations that were used for the main simulation
confs <- unique(all.res[,.(total,ncases)])
confs[,total:=as.character(total) %>% as.numeric]
confs[,ncases:=as.character(ncases) %>% as.numeric]
vars <- lapply(1:nrow(confs),function(i){
  message(i)
  data.table(pc=all.res$pc %>% unique,scale.factor=all.vars$mfactor,ncases=confs$ncases[i],total=confs$total[i])
}) %>% rbindlist
vars[,size.factor:=total/(ncases * (total-ncases))]
vars[,variance:=size.factor * scale.factor]
vars[,total:=factor(total,levels=levels(all.res$total))]

pd <- position_dodge(width=0.3)
ppf <- ggplot(all.res[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) + xlab("Cases") + scale_color_discrete("Sample Size") +
ylab("Variance (PC1 Loadings)") + scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
theme(legend.position =c(0.65,0.85)) + background_grid(major = "xy", minor = "none") +
geom_point(data=vars[pc=='PC1' & total!=350000,],aes(x=as.numeric(ncases),y=variance,color=total,fill=total),position=pd,pch=17,size=5,alpha=0.4,inherit.aes=FALSE) +
scale_fill_discrete(guide=FALSE)
### what happens if we replot using the new routine





vars2.sf<-compute_proj_var2(man.DT,w.DT,shrink.DT,REF_GT_DIR,'ws_emp_shrinkage',quiet=TRUE)
all.vars2 <- data.table(pc=names(vars2.sf),mfactor=vars2.sf)
vars2 <- lapply(1:nrow(confs),function(i){
  message(i)
  data.table(pc=all.res$pc %>% unique,scale.factor=all.vars2$mfactor,ncases=confs$ncases[i],total=confs$total[i])
}) %>% rbindlist
vars2[,size.factor:=total/(ncases * (total-ncases))]
vars2[,variance:=size.factor * scale.factor]
vars2[,total:=factor(total,levels=levels(all.res$total))]

samps<-list.files(path='/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/sims_ws',pattern="wtccc*",full.names=TRUE)

all.res.alt <- lapply(samps,function(path){
  message(path)
  vars <- lapply(list.files(path=path,pattern="^proj*",full.names=TRUE),readRDS) %>% do.call('rbind',.) %>% apply(.,2,var)
  cases <- gsub("wtccc\\_([^_]+)\\_(.*)","\\1",basename(path))
  controls <- gsub("wtccc\\_([^_]+)\\_([^_]+)\\_.*","\\2",basename(path))
  data.table(cases=cases,controls=controls,vars,pc=names(vars))
}) %>% rbindlist
all.res.alt[,total:=as.numeric(cases)+as.numeric(controls)]
all.res.alt[,total:=factor(total,levels=c(levels(all.res$total),3800))]
all.res.alt[,ncases:=as.numeric(cases)]
all.res.alt<-all.res.alt[order(total),]
pd <- position_dodge(width=0.1)
ppf <- ggplot(all.res[pc=='PC1' ,],aes(x=as.numeric(ncases),y=variance,group=total,color=total,ymin=variance-ci.lower, ymax=variance+ci.upper)) +
geom_point(position=pd) + geom_line(position=pd) + geom_errorbar(width=.1,position=pd) + xlab("Cases") + scale_color_discrete("Sample Size") +
ylab("Variance (PC1 Loadings)") + scale_x_continuous(trans="log10",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
theme(legend.position =c(0.65,0.85)) + background_grid(major = "xy", minor = "none") +
geom_point(data=vars[pc=='PC1',],aes(x=as.numeric(ncases),y=variance,color=total,fill=total),position=pd,pch=17,size=5,alpha=0.4,inherit.aes=FALSE) +
geom_point(data=vars2[pc=='PC1',],aes(x=as.numeric(ncases),y=variance,color=total,fill=total),position=pd,pch=18,size=5,alpha=0.4,inherit.aes=FALSE) +
geom_point(data=all.res.alt[pc=='PC1',],aes(x=as.numeric(ncases),y=vars,color=total,fill=total),position=pd,pch=19,size=5,alpha=0.4,inherit.aes=FALSE) +
scale_fill_discrete(guide=FALSE)

## load in bootstrap stuff 
