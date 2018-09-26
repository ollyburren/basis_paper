DATA.DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/'
jia.dat <- file.path(DATA.DIR,'individual_proj/jia_projection.RDS') %>% readRDS
jia.DT <- data.table(ind_id=rownames(jia.dat),jia.dat)
jia.sample.DT <- file.path(DATA.DIR,'ind_sample_info/jia_samples_info.RDS') %>% readRDS

## make a long list

jia.DT <- melt(jia.DT,id.vars='ind_id')[,ind_id:=as.numeric(ind_id)]

M <- merge(jia.DT,jia.sample.DT,by.x='ind_id',by.y='index')[,.(variable,value,sample_id,ilar_final)]

by.dx <- M[,list(mean.load=mean(value)),by=c('variable','ilar_final')]

library(cupcake)
SNP_MANIFEST_FILE <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/snp_manifest/ichip_september.tab'
ICHIP_DATA_DIR <- '/home/ob219/rds/hpc-work/as_basis/gwas_stats/ichip/aligned/'
MANIFEST_FILE <- '/home/ob219/git/as_basis/manifest/as_manifest_ichip.tsv'
SHRINKAGE_METHOD<-'ws_emp'
basis.DT<-get_gwas_data(MANIFEST_FILE,SNP_MANIFEST_FILE,ICHIP_DATA_DIR,filter_snps_by_manifest=TRUE)
shrink.DT<-compute_shrinkage_metrics(basis.DT)
## need to add control where beta is zero
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
pc.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
control.DT <- pc.DT[trait=='control',]

by.dx <- by.dx[,.(trait=ilar_final,variable,value=mean.load)]

plot.DT <- merge(by.dx,control.DT,by.x='variable',by.y='variable')[,.(value=value.x-value.y,variable,trait=trait.x)]
plot.DT <- plot.DT[!trait %in% c('JIA_unknown','JIA_undefined')]

#plot.DT <- rbind(by.dx,control.DT)



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
idx <- which(jia.sample.DT$NewDx %in% c('JDM','DM'))
## get a list of study centres that have no cases for either JDM or DM
cod <- (table(jia.sample.DT[idx,]$origin,jia.sample.DT[idx,]$NewDx)==0) %>% rowSums
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
