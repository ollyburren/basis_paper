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
ukbb <- merge(ukbb,anno[,.(phe,n1,n0)],by.x='trait',by.y='phe')
ukbb <- ukbb[,trait:=paste('bb',trait,sep='_')]


## load in jia sub types

jia <- readRDS("/home/ob219/share/as_basis/GWAS/jia_projections/summary/jia.RDS")
## get sample counts and merge
library(annotSnpStats)
(load('~/share/Data/GWAS/JIA-2017-data/annotsnpstats-22.RData'))
samp.DT <- samples(G) %>% data.table
samp.DT <- samp.DT[,list(n1=.N),by=alt_ilar_code]
samp.DT <- samp.DT[,n0:=5181][!is.na(alt_ilar_code),]
samp.DT[,alt_ilar_code:=paste('jia',alt_ilar_code,sep='_')]
jia <- melt(jia,id.var='trait')
jia <- merge(jia,samp.DT,by.x='trait',by.y='alt_ilar_code')

## load in tian_infectious_disease

tian <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/tian_infectious_disease.RDS")
## get sample counts and merge
tian_samples <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/sample_size.RDS")
tian <- melt(tian,id.var='trait')
tian <- merge(tian,tian_samples,by='trait')

## load in ferreira_asthma

ferriera <- readRDS("/home/ob219/share/as_basis/GWAS/ferreira_projections/ferreira_asthma.RDS")
ferreira_samples <- readRDS("/home/ob219/share/as_basis/GWAS/ferreira_projections/sample_size.RDS")
ferriera <- melt(ferriera,id.var='trait')
ferriera <- merge(ferriera,ferreira_samples,by='trait')


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

## load in ig titres

eff <- readRDS('/home/ob219/share/as_basis/GWAS/effrosyni_ig/effrosyni_ig.RDS')
eff <- melt(eff,id.var='trait')
eff[,c('n0','n1'):=list(3969,3969)]

## load in CD prognosis
cd.prog <- readRDS('/home/ob219/share/as_basis/GWAS/cd_prognosis/cd_prognosis.RDS')
cd.prog <- data.table(trait=rownames(cd.prog),cd.prog)
cd.prog <- melt(cd.prog,id.var='trait')
cd.prog[,c('n0','n1'):=list(389 + 583,669 + 1093)]

## load in bp,adhd, and scz

files <- list.files(path="/home/ob219/share/as_basis/GWAS/psych/",pattern="*.RDS",full.names=TRUE)

psy <- lapply(files,readRDS) %>% do.call('rbind',.)
psy <- data.table(trait=rownames(psy),psy)
psy <- melt(psy,id.var='trait')
psy[trait=='ADHD',c('n0','n1'):=list(34194,19099)]
psy[trait=='SCZ',c('n0','n1'):=list(32541,33426)]
psy[trait=='BIP',c('n0','n1'):=list(21524,20129)]

## myogen_myositis
myogen <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis.RDS')
myogen <- data.table(trait=rownames(myogen),myogen)
myogen <- melt(myogen,id.var='trait')
myogen[,c('n0','n1'):=list(4726,1818)]

## NMO - Neuromyelitis Optica
files <- list.files(path="/home/ob219/share/as_basis/GWAS/nmo/",pattern="*.RDS",full.names=TRUE)
nmo <- lapply(files,readRDS) %>% do.call('rbind',.)
#nmo <- readRDS('/home/ob219/share/as_basis/GWAS/nmo/nmo.RDS')
nmo <- data.table(trait=rownames(nmo),nmo)
nmo <- melt(nmo,id.var='trait')
nmo[trait=='NMO_combined',c('n0','n1'):=list(1244,215)]
nmo[trait=='NMO_IgGPos',c('n0','n1'):=list(1244,66+66)]
nmo[trait=='NMO_IgGNeg',c('n0','n1'):=list(1244,20+63)]

## load in egpa vasc

files <- list.files(path="/home/ob219/share/as_basis/GWAS/vasc/projections",pattern="*.RDS",full.names=TRUE)
vasc <- lapply(files,readRDS) %>% do.call('rbind',.)
vasc <- data.table(trait=rownames(vasc),vasc)
vasc <- melt(vasc,id.var='trait')
vasc.samp <- readRDS('/home/ob219/share/as_basis/GWAS/vasc/sample.RDS')
vasc <- merge(vasc,vasc.samp,by='trait')

## load in abdef

abdef <- readRDS('/home/ob219/share/as_basis/GWAS/abdef/abdef.RDS')
abdef <- data.table(trait=rownames(abdef),abdef)
abdef <- melt(abdef,id.var='trait')
abdef[,c('n0','n1'):=list(9225,733)]


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
  vasc=vasc,
  abdef=abdef
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

## how many traits have a pc score that is less than FDR 5%
keep.traits <- all.DT[p.adj<0.05,]$trait %>% unique
out.DT <- all.DT[trait %in% keep.traits,]
setnames(out.DT,'variable','PC')
#out.DT[,delta:=value-control.loading]
## vitamin c is in twice remove !
out.DT <- out.DT[trait!='bb_vitamin.c.product',]
## work out hclust order
deltas <- melt(out.DT,id.vars=c('trait','PC'),measure.vars='delta') %>% dcast(.,trait~PC+variable)
delta.mat <- as.matrix(deltas[,-1])
rownames(delta.mat) <- deltas$trait
delta.hc <- dist(delta.mat) %>% hclust
out.DT[,trait:=factor(trait,levels=delta.hc$labels[delta.hc$order])]
out.DT[,PC:=factor(PC,levels=paste('PC',1:11,sep=""))]

library(cowplot)
ggplot(out.DT[p.adj<0.05 & !PC %in% c('PC10','PC11'),],aes(x=trait,y=PC,fill= (pmax(Z,-11) %>% pmin(.,11)),label=signif(delta,digits=1))) + geom_tile(color='black') +
geom_text(size=5,angle=90) + scale_fill_gradientn("Z",colours=c('green','white','red')) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#dev.print(pdf,"~/tmp/all_cc.pdf")

#out.DT[p.adj>0.05, Z:=0]

#ggplot(out.DT[PC != 'PC11',],aes(x=trait,y=PC,fill= (pmax(Z,-11) %>% pmin(.,11)),label=signif(delta,digits=1))) + geom_tile(color='black') +
#geom_text(size=5,angle=90) + scale_fill_gradientn("Z",colours=c('green','white','red')) +
#theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## print a similar for projected traits

## work out hclust order for basis
if(FALSE){
setnames(basis.DT,'variable','PC')
deltas.basis <- melt(basis.DT,id.vars=c('trait','PC'),measure.vars='value') %>% dcast(.,trait~PC+variable)
delta.basis.mat <- as.matrix(deltas.basis[,-1])
rownames(delta.basis.mat) <- deltas.basis$trait
delta.basis.hc <- dist(delta.basis.mat) %>% hclust
basis.DT[,trait:=factor(trait,levels=delta.basis.hc$labels[delta.basis.hc$order])]
basis.DT[,PC:=factor(PC,levels=paste('PC',1:11,sep=""))]



ggplot(basis.DT[trait!='control',],aes(x=trait,y=PC,fill= Z,label=signif(value,digits=1))) + geom_tile(color='black') +
geom_text(size=5,angle=90) + scale_fill_gradientn("Z",colours=c('green','white','red')) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.print(pdf,"~/tmp/basis.pdf")
}
