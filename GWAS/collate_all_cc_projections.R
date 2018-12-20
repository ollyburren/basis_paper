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
jia[,category:='bowes_jia']

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

## myogen_myositis
myogen <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis.RDS')
myogen <- data.table(trait=rownames(myogen),myogen)
myogen <- melt(myogen,id.var='trait')
myogen[,c('n0','n1'):=list(4726,1818)]
myogen[,category:='myogen']

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

files <- list.files(path="/home/ob219/share/as_basis/GWAS/vasc/projections",pattern="*.RDS",full.names=TRUE)
vasc <- lapply(files,readRDS) %>% do.call('rbind',.)
vasc <- data.table(trait=rownames(vasc),vasc)
vasc <- melt(vasc,id.var='trait')
vasc.samp <- readRDS('/home/ob219/share/as_basis/GWAS/vasc/sample.RDS')
vasc <- merge(vasc,vasc.samp,by='trait')
vasc[,category:='lyons_vasculitis']

## load in abdef

abdef <- readRDS('/home/ob219/share/as_basis/GWAS/abdef/abdef.RDS')
abdef <- data.table(trait=rownames(abdef),abdef)
abdef <- melt(abdef,id.var='trait')
abdef[,c('n0','n1'):=list(9225,733)]
abdef[,category:='ad-pid']

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
  #myogen=rbind(myogen,myogen.flip),
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

saveRDS(all.DT,'/home/ob219/share/as_basis/GWAS/RESULTS/19_12_18_summary_results.RDS')


if(FALSE){
## how many traits have a pc score that is less than FDR 5%
keep.traits <- all.DT[p.adj<0.05 & !variable %in% c('PC10','PC11'),]$trait %>% unique
out.DT <- all.DT[trait %in% keep.traits,]
setnames(out.DT,'variable','PC')
#out.DT[,delta:=value-control.loading]
## vitamin c is in twice remove !
out.DT <- out.DT[trait!='bb_vitamin.c.product',]
## work out hclust order

get_hclust_order <- function(DT){
  deltas <- melt(DT,id.vars=c('trait','PC'),measure.vars='delta') %>% dcast(.,trait~PC+variable)
  delta.mat <- as.matrix(deltas[,-1])
  rownames(delta.mat) <- deltas$trait
  delta.hc <- dist(delta.mat) %>% hclust
  delta.hc$labels[delta.hc$order]
}
out.DT[,PC:=factor(PC,levels=paste('PC',1:11,sep=""))]
med.DT <- out.DT[category=='medications' & !PC %in% c('PC10','PC11'),]
med.DT[,trait:=factor(trait,levels=get_hclust_order(med.DT))]
med.DT[p.adj>=0.05,delta:=NA]
disease.DT <- out.DT[category!='medications' & !PC %in% c('PC10','PC11'),]
disease.DT[,trait:=factor(trait,levels=get_hclust_order(disease.DT))]
disease.DT[p.adj>=0.05,delta:=NA]

# deltas <- melt(out.DT,id.vars=c('trait','PC'),measure.vars='delta') %>% dcast(.,trait~PC+variable)
# delta.mat <- as.matrix(deltas[,-1])
# rownames(delta.mat) <- deltas$trait
# delta.hc <- dist(delta.mat) %>% hclust

#out.DT[,trait:=factor(trait,levels=delta.hc$labels[delta.hc$order])]


library(cowplot)
#ggplot(out.DT[p.adj<0.05 & !PC %in% c('PC10','PC11'),],aes(x=trait,y=PC,fill= (pmax(Z,-11) %>% pmin(.,11)),label=signif(delta,digits=1))) + geom_tile(color='black') +
pa <- ggplot(med.DT,aes(x=trait,y=PC,fill= (pmax(Z,-11) %>% pmin(.,11)),label=signif(delta,digits=1))) + geom_tile(color='black') +
geom_text(size=5,angle=90) + scale_fill_gradientn("Z",colours=c('green','white','red')) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill=FALSE) + ggtitle("Medication")

pb <- ggplot(disease.DT,aes(x=trait,y=PC,fill= (pmax(Z,-11) %>% pmin(.,11)),label=signif(delta,digits=1))) + geom_tile(color='black') +
geom_text(size=5,angle=90) + scale_fill_gradientn("Z",colours=c('green','white','red')) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +  ggtitle("Disease") +
theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plot_grid(pa,pb)

stop()

## category stuff
ss.plot <- out.DT[PC=='PC1',.(trait,n0=n0,n1=n1,category)]
ss.plot <- melt(ss.plot,id.vars=c('trait','category'))

#ggplot(ss.plot[variable=='n1',],aes(x=trait,y=value,fill=category)) + geom_bar(stat="identity",col='black') +
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  scale_y_log10(
#    breaks = scales::trans_breaks("log10", function(x) 10^x),
#    labels = scales::trans_format("log10", scales::math_format(10^.x))
#  )

## we are interested in this which appears to separate mpo_Pos and anca_Neg
int<-all.DT[trait %in% c('mpo_Pos','anca_Neg') & variable=='PC6',]
int[,ci:=1.96 * sqrt(variance)]

ggplot(int,aes(x=trait,y=value)) + geom_point() +  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.1,)


T.test <- function(n, mean, sd) {
  s <- sum((n - 1) * sd^2) / (sum(n) - 2) # weighted variance
  t <- sqrt(prod(n) / sum(n)) * (diff(mean) / sqrt(s)) # t statistic
  df <- sum(n) - 2  # degrees of freedom
  p <- (1 - pt(abs(t), df)) * 2 # p value
  c(t = t, p = p)
}

## alternative method just simulate

m1=int$value[1]; s1=int$variance[1]; n1=800
m2=int$value[2]; s2=int$variance[2]; n2=800

# method 1: normal samples, rescaled to match original:
z1 = rnorm(n1)
x1 = scale(z1)*s1+m1
z2 = rnorm(n2)
x2 = scale(z2)*s2+m2
t.test(x1,x2) # Welch-Satterthwaite
plot.it <- rbind(data.table(trait=int$trait[1],loadings=x1[,1]),data.table(trait=int$trait[2],loadings=x2[,1]))

ggplot(plot.it,aes(x=trait,y=loadings)) + geom_boxplot()
}



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
