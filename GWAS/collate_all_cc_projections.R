## overall file to collate all projection results

## load in UKBB results

BB_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/ss_shrink_2018'
ukbb <- lapply(list.files(path=BB_DIR,pattern="*.RDS",full.names=TRUE),readRDS) %>% rbindlist
ukbb <- melt(ukbb,id.var='trait')
## need to add case and control numbers
get_phenotype_annotation <- function(filter='20002\\_'){
  bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
  pheno <- fread(bb_phenofile)
  setnames(pheno,names(pheno) %>% make.names)

  #med <- pheno[grepl("2000[23]\\_",Phenotype.Code) & Sex=='both_sexes',]
  med <- pheno[grepl(filter,Phenotype.Code) & Sex=='both_sexes',]
  ## load in phenotype file
  P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/phenotypes.both_sexes.tsv')
  P<-P[,.(phenotype,variable_type,non_missing=n_non_missing,cases=n_cases,controls=n_controls)]
  med[,phe:=make.names(Phenotype.Description) %>% gsub("(Non.cancer.illness.code..self.reported..)|(Treatment.medication.code..)|(Cancer.code..self.reported..)","",.)]
  med <- med[,.(Phenotype.Code,phe)]
  med <- merge(P,med,by.x='phenotype',by.y='Phenotype.Code')[,.(code=phenotype,phe,n1=as.numeric(cases),n0=as.numeric(controls))]
  med
}

anno <- get_phenotype_annotation(".*")[!is.na(n1) & !is.na(n0) & phe != 'unclassifiable',]
## compose an annotation bar

bb_cod.sr_disease <- fread("/home/ob219/share/as_basis/GWAS/bb_projections/coding6.tsv")
bb_cod.sr_cancer <- fread("/home/ob219/share/as_basis/GWAS/bb_projections/coding3.tsv")
bb_cod.sr_medication <- fread("/home/ob219/share/as_basis/GWAS/bb_projections/coding4.tsv")

ukbb.traits <- ukbb$trait %>% unique %>% gsub("^bb\\_","",.)
anno.bb <- anno[phe %in% ukbb.traits]
anno.bb[,c('clade','coding'):=tstrsplit(code,'_') %>% lapply(.,as.numeric)]
## note only self reported data has a coding, drugs and cancer must have a separate classification
## and there are collisions so need to exclude
#root.bb <- merge(anno.bb[clade=='20002',],bb_cod,by='coding')
## fill in parent nodes

get_hier <- function(DT,id,l,top.node=TRUE){
  #message(head(DT))
  if(missing(l)){
    l <- DT[coding==id,]$meaning
    parent_id <- DT[coding==id,]$parent_id
  }else{
    l <- c(l,DT[node_id==id,]$meaning)
    parent_id <- DT[node_id==id,]$parent_id
  }

  if(length(parent_id)==0 || parent_id==0){
    if(top.node){
      message(sprintf("Top node %d",parent_id))
      tail(l,n=1) %>% return
    }else{
      paste(l,sep=':',collapse=':') %>% return
    }
  }else{
    get_hier(DT,parent_id,l,top.node)
  }
}

#root.bb[,bread.crumbs:=sapply(node_id,get_hier,FALSE)]
#srd <- anno.bb[clade=='20002',]
#srd[,top.node:=sapply(coding,get_hier,DT=bb_cod.sr_disease,top.node=TRUE)]
#anno.bb[clade=='20002',top.node:=sapply(coding,get_hier,DT=bb_cod.sr_disease,top.node=TRUE)]
#anno.bb[clade=='20001',top.node:=sapply(coding,get_hier,DT=bb_cod.sr_cancer,top.node=TRUE)]
anno.bb[clade=='20002',top.node:='bb_disease']
anno.bb[clade=='20001',top.node:='bb_cancer']
anno.bb[clade=='20003',top.node:='bb_medications']

ukbb <- merge(ukbb,anno.bb[,.(phe,n1,n0,category=top.node)],by.x='trait',by.y='phe')
ukbb <- ukbb[,trait:=paste('bb',trait,sep='_')]

##remove medication for the time being
#ukbb <- ukbb[!is.na(category),]

## load in jia sub types

jia <- readRDS("/home/ob219/share/as_basis/GWAS/jia_projections/summary/jia_2019.RDS")
sample.DT <- fread('/home/ob219/share/Data/GWAS/jia-mar-2019/summary-stats-samplecount-mar2019.csv')
setnames(sample.DT,'n','n1')
sub.DT <- data.table(idx=0:9,subtype=c('case','sys','PO','EO','RFneg','RFpos','ERA','PsA','undiff','missing'))
samp.DT <- merge(sub.DT,sample.DT[,.(ilar_pheno,n1)],by.x='idx',by.y='ilar_pheno')
samp.DT[,n0:=9196]
jia <- melt(jia,id.var='trait')
samp.DT[,subtype:=sprintf("jia_%s_19",subtype)]
jia <- merge(jia,samp.DT[,.(subtype,n1,n0)],by.x='trait',by.y='subtype')
jia[,category:='bowes_jia_2019']

## code below for the old data
#jia <- readRDS("/home/ob219/share/as_basis/GWAS/jia_projections/summary/jia.RDS")
### get sample counts and merge
#library(annotSnpStats)
#(load('~/share/Data/GWAS/JIA-2017-data/annotsnpstats-22.RData'))
#samp.DT <- samples(G) %>% data.table
#samp.DT <- samp.DT[,list(n1=.N),by=alt_ilar_code]
#samp.DT <- samp.DT[,n0:=5181][!is.na(alt_ilar_code),]
#samp.DT[,alt_ilar_code:=paste('jia',alt_ilar_code,sep='_')]
#jia <- melt(jia,id.var='trait')
#jia <- merge(jia,samp.DT,by.x='trait',by.y='alt_ilar_code')
#jia[,category:='bowes_jia']

## load in tian_infectious_disease

tian <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/tian_infectious_disease.RDS")
## get sample counts and merge
tian_samples <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/sample_size.RDS")
tian <- melt(tian,id.var='trait')
tian <- merge(tian,tian_samples,by='trait')
tian[,category:='tian_infectious_disease']

## load in ferreira_asthma

ferriera <- readRDS("/home/ob219/share/as_basis/GWAS/ferreira_projections/ferreira_asthma.RDS")
ferreira_samples <- readRDS("/home/ob219/share/as_basis/GWAS/ferreira_projections/sample_size.RDS")
ferriera <- melt(ferriera,id.var='trait')
ferriera <- merge(ferriera,ferreira_samples,by='trait')
ferriera[,category:='ferreira_asthma']


## load in astle data

## first get sample sizes - note that because of assumptions in converting linear coefficient
## to OR the cases and control size are equal to the median/mean of the sample size i.e round(ss/2)
ASTLE_DIR <- '/home/ob219/share/as_basis/GWAS/astle'
astle <- lapply(list.files(path=ASTLE_DIR,pattern="*.RDS",full.names=TRUE),readRDS) %>% do.call('rbind',.)
astle <- data.table(trait=rownames(astle),astle)
astle <- melt(astle,id.var='trait')


DATA_DIR <- '/home/ob219/share/Data/GWAS-summary/blood-ukbiobank-2016-12-12'
afiles <- list.files(path=DATA_DIR,pattern="*.gz$",full.names=FALSE)
astle_samples <- data.table(trait = gsub("(.*)\\_build37\\_[0-9]+\\_20161212.tsv.gz","\\1",afiles),
ss = as.numeric(gsub("(.*)\\_build37\\_([0-9]+)\\_20161212.tsv.gz","\\2",afiles)))
astle_samples[,c('n0','n1','ss'):=list(round(ss/2),round(ss/2),NULL)]
astle <- merge(astle,astle_samples)
astle[,category:='astle_blood']

## load in ig titres

eff <- readRDS('/home/ob219/share/as_basis/GWAS/effrosyni_ig/effrosyni_ig.RDS')
eff <- melt(eff,id.var='trait')
eff[,c('n0','n1'):=list(3969,3969)]
eff[,category:='Ig.titre']

## load in CD prognosis
cd.prog <- readRDS('/home/ob219/share/as_basis/GWAS/cd_prognosis/cd_prognosis.RDS')
cd.prog <- data.table(trait=rownames(cd.prog),cd.prog)
cd.prog <- melt(cd.prog,id.var='trait')
cd.prog[,c('n0','n1'):=list(389 + 583,669 + 1093)]
cd.prog[,category:='lee_CD_prognosis']

## load in bp,adhd, and scz

files <- list.files(path="/home/ob219/share/as_basis/GWAS/psych/",pattern="*.RDS",full.names=TRUE)

psy <- lapply(files,readRDS) %>% do.call('rbind',.)
psy <- data.table(trait=rownames(psy),psy)
psy <- melt(psy,id.var='trait')
psy[trait=='ADHD',c('n0','n1'):=list(34194,19099)]
psy[trait=='SCZ',c('n0','n1'):=list(32541,33426)]
psy[trait=='BIP',c('n0','n1'):=list(21524,20129)]
psy[,category:='psyc_consortium']

## have to calculate case and control numbers from samples

# jdm <- fread("jdmall.txt",header=FALSE)
# pm <- fread("pmall.txt",header=FALSE)
# dm <- fread("dmall.txt",header=FALSE)
# ## expect controls to be shared between traits
# control.ids <- union(intersect(dm$V1,jdm$V1),intersect(pm$V1,dm$V1))
# jdm[,affected:=0]
# jdm[!V1 %in% control.ids,affected:=1]
# dm[,affected:=0]
# dm[!V1 %in% control.ids,affected:=1]
# pm[,affected:=0]
# pm[!V1 %in% control.ids,affected:=1]

## myogen_myositis
myogen <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis.RDS')
myogen <- data.table(trait=rownames(myogen),myogen)
myogen <- melt(myogen,id.var='trait')
myogen[,c('n0','n1'):=list(4724,1711)]
myogen[,category:='myogen']

## jdm

jdm <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/jdm_myositis.RDS')
jdm <- data.table(trait=rownames(jdm),jdm)
jdm <- melt(jdm,id.var='trait')
jdm[,c('n0','n1'):=list(4724,473)]
jdm[,category:='myogen']

dm <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/dm_myositis.RDS')
dm <- data.table(trait=rownames(dm),dm)
dm <- melt(dm,id.var='trait')
dm[,c('n0','n1'):=list(4724,705)]
dm[,category:='myogen']

pm <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/pm_myositis.RDS')
pm <- data.table(trait=rownames(pm),pm)
pm <- melt(pm,id.var='trait')
pm[,c('n0','n1'):=list(4724,533)]
pm[,category:='myogen']

myogen <- rbindlist(list(myogen,jdm,dm,pm))

## myogen_myositis
# myogen.flip <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis_flip.RDS')
# myogen.flip <- data.table(trait=rownames(myogen.flip),myogen.flip)
# myogen.flip <- melt(myogen.flip,id.var='trait')
# myogen.flip[,c('n0','n1'):=list(4726,1818)]
# myogen.flip[,category:='myogen']

## NMO - Neuromyelitis Optica
files <- list.files(path="/home/ob219/share/as_basis/GWAS/nmo/",pattern="*.RDS",full.names=TRUE)
nmo <- lapply(files,readRDS) %>% do.call('rbind',.)
#nmo <- readRDS('/home/ob219/share/as_basis/GWAS/nmo/nmo.RDS')
nmo <- data.table(trait=rownames(nmo),nmo)
nmo <- melt(nmo,id.var='trait')
nmo[trait=='NMO_combined',c('n0','n1'):=list(1244,215)]
nmo[trait=='NMO_IgGPos',c('n0','n1'):=list(1244,66+66)]
nmo[trait=='NMO_IgGNeg',c('n0','n1'):=list(1244,20+63)]
nmo[,category:='estrada_NMO']

## load in egpa vasc

files <- list.files(path="/home/ob219/share/as_basis/GWAS/lyons_egpa/projections",pattern="*.RDS",full.names=TRUE)
egpa <- lapply(files,readRDS) %>% do.call('rbind',.)
egpa <- data.table(trait=rownames(egpa),egpa)
egpa <- melt(egpa,id.var='trait')
egpa.samp <- readRDS('/home/ob219/share/as_basis/GWAS/lyons_egpa/sample.RDS')
egpa <- merge(egpa,egpa.samp,by='trait')
egpa[,category:='lyons_egpa']


## load in aav (vasculitis) from Limy Wong

aav <- readRDS("/home/ob219/share/as_basis/GWAS/wong_aav/projections/aav_2019.RDS")
aav <- melt(aav,id.var='trait')
aav.samp <- fread("/home/ob219/share/Data/GWAS-summary/aav_limy_wong/aav_sample_size.txt")
aav <- merge(aav,aav.samp[,.(label,n0,n1)],by.x='trait',by.y='label')
aav[,category:='wong_aav']


## load in abdef

abdef <- readRDS('/home/ob219/share/as_basis/GWAS/abdef/abdef.RDS')
abdef <- data.table(trait=rownames(abdef),abdef)
abdef <- melt(abdef,id.var='trait')
abdef[,c('n0','n1'):=list(9225,733)]
abdef[,category:='ad-pid']

## load in t1d ab data
t1d <- readRDS('/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/liley_t1d/proj.RDS')
## remove age as this would need rescaling to be able to interpret
t1d <- t1d[trait!='z_age']
#t1d <- data.table(trait=rownames(t1d),t1d)
t1d <- melt(t1d,id.var='trait')

t1d[trait=='z_t1d',c('n0','n1'):=list(8825,5908)]
## data from Vincent Plagnol paper
t1d[trait=='z_tpo',c('n0','n1'):=list(5908*(1-0.12),5908*0.12)]
t1d[trait=='z_gad',c('n0','n1'):=list(3208*(1-0.5),3208*0.5)]
t1d[trait=='z_ia2',c('n0','n1'):=list(3197*(1-0.59),3197*0.59)]
t1d[trait=='z_pca',c('n0','n1'):=list(2240*(1-0.1),2240*0.1)]

t1d[,category:='liley_t1d']

psa <- readRDS("/home/ob219/share/as_basis/GWAS/psa_projections/summary/bowes_psa.RDS")
psa <- melt(psa,id.var='trait')
psa[,c('n0','n1','category'):=list(4596,1805,'bowes_psa')]

## for some reason we swapped from n1 to n0 halfway through to n0 n1
## we need to fix otherwise everything gets swapped and case
## sizes become control sizes

ncolorder<-c('trait','variable','value','n0','n1')
setcolorder(ferriera,ncolorder)
setcolorder(tian,ncolorder)
setcolorder(jia,ncolorder)
setcolorder(ukbb,ncolorder)


all.proj <- list(
  ferriera=ferriera,
  tian=tian,
  jia=jia,
  ukbb=ukbb,
  astle=astle,
  eff=eff,
  cd.prog=cd.prog,
  psy=psy,
  myogen=myogen,
  nmo=nmo,
  egpa=egpa,
  abdef=abdef,
  liley=t1d,
  aav=aav,
  psa=psa
) %>% rbindlist
all.proj[,n:=n1+n0]

## load in basis and variance
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_june.RDS'
var.DT <- readRDS(VARIANCE_FILE)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,value:=(value-control.loading)]
basis.DT[,Z:=sign(value)]
basis.DT[,p.adj:=1]
basis.DT[,control.loading:=NULL]

## compute the variance of a projection
all.DT <- merge(all.proj,var.DT,by.x='variable',by.y='pc')
all.DT <- merge(all.DT,control.DT,by.x='variable',by.y='PC')
all.DT[,variance:=((log(n)-(log(n1) + log(n-n1)))+ log(mfactor)) %>% exp]
## add in control loading
all.DT[,Z:=(value-control.loading)/sqrt(variance)]
all.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
all.DT[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']
all.DT[,delta:=value-control.loading]

saveRDS(all.DT,'/home/ob219/share/as_basis/GWAS/RESULTS/09_04_19_summary_results.RDS')
