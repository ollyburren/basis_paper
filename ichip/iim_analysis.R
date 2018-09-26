DATA.DIR <- '/home/ob219/share/as_basis/ichip/'
iim.dat <- file.path(DATA.DIR,'individual_projections/iim_projection.RDS') %>% readRDS
iim.DT <- data.table(ind_id=rownames(iim.dat),iim.dat)
iim.sample.DT <- file.path(DATA.DIR,'ind_sample_info/iim_samples_info.RDS') %>% readRDS

## make a long list

iim.DT <- melt(iim.DT,id.vars='ind_id')

M <- merge(iim.DT,iim.sample.DT,by.x='ind_id',by.y='ID')

by.dx <- M[,list(mean.load=mean(value),sample.count=.N),by=c('variable','NewDx')]

library(cupcake)
BASIS_FILE <- '/home/ob219/share/as_basis/ichip/support/basis_ic.RDS'
pc.emp <- readRDS(BASIS_FILE)


# SNP_MANIFEST_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/snp_manifest/august_iim.tab'
# ICHIP_DATA_DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/aligned/'
# MANIFEST_FILE <- '/home/ob219/git/as_basis/manifest/as_manifest_ichip.tsv'
# SHRINKAGE_METHOD<-'ws_emp'
# basis.DT<-get_gwas_data(MANIFEST_FILE,SNP_MANIFEST_FILE,ICHIP_DATA_DIR,filter_snps_by_manifest=TRUE)
# shrink.DT<-compute_shrinkage_metrics(basis.DT)
# ## need to add control where beta is zero
# basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
# ## need to add control where beta is zero
# basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
# pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
pc.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
control.DT <- pc.DT[trait=='control',]

by.dx <- by.dx[,.(trait=NewDx,variable,value=mean.load)]
plot.DT <- rbind(by.dx,control.DT)



library(cowplot)
ggplot(plot.DT,aes(x=variable,y=value,color=trait,group=trait)) + geom_point() + geom_line(alpha=0.5,lty=2)

## quick test for difference between JDM and DM

lapply(paste('PC',1:9,sep=''),function(pc){
  tmp <- M[variable==pc,]
  t.test(tmp[NewDx=='JDM',]$value,tmp[NewDx=='DM',]$value)$p.value
})

source('https://raw.githubusercontent.com/chr1swallace/random-functions/master/R/Cochran-Armitage.test.mod.R')
## need to remove strata for which there is no data
## from sample manifest select ones we are interested in
idx <- which(iim.sample.DT$NewDx %in% c('JDM','DM'))
## get a list of study centres that have no cases for either JDM or DM
cod <- (table(iim.sample.DT[idx,]$origin,iim.sample.DT[idx,]$NewDx)==0) %>% rowSums
# M contains a merge between individual projection loadings and sample file
# we filter to get traits to compare and those studies with at least one case
M.filt <- M[origin %in% names(cod)[cod!=1] & NewDx %in% c('JDM','DM'),]

## perform CA
lapply(paste('PC',1:9,sep=''),function(pc){
  idx <- which(M.filt$variable==pc)
  expo <- M.filt[idx,]$value ## PC projection loading from individual
  ccd <- M.filt[idx,]$NewDx ## Whether JDM or DM
  stratum <- M.filt[idx,]$origin ## Study centre - note centres with no cases for a particular type are prefiltered
  data.table(pc=pc,ca.pval=Cochran.Armitage.test(expo,ccd,stratum)$p.value)
}) %>% rbindlist


## there appears to be no real difference - thus for future analysis lump JDM and DM and compare with PM
M.c <- copy(M)
M.c[NewDx %in% c('JDM','DM'),NewDx:='JDM_DM']
idx <- which(M.c$NewDx %in% c('JDM_DM','PM') & M.c$variable=='PC1')
## get a list of study centres that have no cases for either JDM,DM or PM
cod <- (table(M.c[idx,]$origin,M.c[idx,]$NewDx)==0) %>% rowSums
# M contains a merge between individual projection loadings and sample file
# we filter to get traits to compare and those studies with at least one case
M.filt <- M.c[origin %in% names(cod)[cod!=1] & NewDx %in% c('JDM_DM','PM'),]

## perform CA
lapply(paste('PC',1:9,sep=''),function(pc){
  idx <- which(M.filt$variable==pc)
  expo <- M.filt[idx,]$value ## PC projection loading from individual
  ccd <- M.filt[idx,]$NewDx ## Whether JDM or DM
  stratum <- M.filt[idx,]$origin ## Study centre - note centres with no cases for a particular type are prefiltered
  data.table(pc=pc,ca.pval=Cochran.Armitage.test(expo,ccd,stratum)$p.value)
}) %>% rbindlist


## alternative method is to use control data to estimate distribution and convert to z score

ctrl <- readRDS("/home/ob219/share/as_basis/ichip/individual_projections/ctrl_projection.RDS")
ctrl.DT <- data.table(ind_id=rownames(ctrl),ctrl)
ctrl.DT <- melt(data.table(ctrl.DT),id.vars='ind_id')
ctrl.DT[,c('fid','origin','NewDx'):=list(1,'control','control')]

all <- rbind(ctrl.DT,M)[,category:=ifelse(NewDx=='control','control','iim')]
bla<-dcast(all,ind_id+NewDx~variable)
#bla<-dcast(all,ind_id+origin~variable)

pc.emp.DT <- data.table(trait=pc.emp$x %>% rownames,pc.emp$x)

#ggplot(bla,aes(x=PC1,y=PC2,color=origin)) + geom_point(alpha=0.5) + geom_point(data=pc.emp.DT,col='black') +
#geom_text(data=pc.emp.DT,aes(x=PC1,y=PC2,label=trait),inherit.aes=FALSE)


#ggplot(bla[origin %in% bla[,.(.N),by=origin][N>100,]$origin,],aes(x=PC1,y=PC2,color=origin)) + geom_point(alpha=0.5) + geom_point(data=pc.emp.DT,col='black') +
#geom_text(data=pc.emp.DT,aes(x=PC1,y=PC2,label=trait),inherit.aes=FALSE)

ggplot(bla,aes(x=PC1,y=PC2,color=NewDx)) + geom_point(alpha=0.5) + geom_point(data=pc.emp.DT,col='black') +
geom_text(data=pc.emp.DT,aes(x=PC1,y=PC2,label=trait),inherit.aes=FALSE)




plot.DT <- melt(all,id.vars=c('ind_id','category'),measure.vars='value') %>% dcast(ind_id+category~variab)




paste('PC',1:13,sep=''))

ggplot(all[variable=='PC1',],aes(x=category,y=value)) + geom_jitter()


proj.mean <- as.matrix(ctrl.DT) %>% apply(.,2,mean)
proj.sd <- as.matrix(ctrl.DT) %>% apply(.,2,var)
null.ss <- data.table(variable=names(proj.mean),proj.mean=proj.mean,proj.sd=proj.sd)
merge(by.dx,null.ss,by='variable')[,.(variable,NewDx,Z=(mean.load-proj.mean)/proj.sd)]

## bootstraps mean and variance estimation -not fast but works
getDistro <- function(ctrl.DT,n,sims=1000){
  res <- lapply(1:sims,function(i){
    tmp <- ctrl.DT[sample.int(nrow(ctrl.DT),n),]
    proj.mean <- as.matrix(tmp) %>% apply(.,2,mean)
    proj.sd <- as.matrix(tmp) %>% apply(.,2,var)
    data.table(variable=names(proj.mean),proj.mean=proj.mean,proj.sd=proj.sd)
  }) %>% rbindlist
  res[,list(e.proj.mean=mean(proj.mean),e.proj.sd=mean(proj.sd),n=n),by='variable']
}

foo <- lapply(split(by.dx,by.dx$NewDx),function(x){
  getDistro(ctrl.DT,x$sample.count %>% unique)
}) %>% rbindlist


ggplot(foo[variable=='PC1',],aes(y=e.proj.sd,x=n,group=variable)) + geom_point() + geom_line()
