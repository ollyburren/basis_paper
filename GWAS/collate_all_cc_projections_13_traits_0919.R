## overall file to collate all projection results

## load in UKBB results

ukbb <- readRDS('/home/ob219/share/as_basis/GWAS/RESULTS/uk_bb_for_fdr_13_traits_0919.RDS')
ukbb <- data.table(trait=rownames(ukbb),ukbb)
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
  #med[,phe:=make.names(Phenotype.Description) %>% gsub("(Non.cancer.illness.code..self.reported..)|(Treatment.medication.code..)|(Cancer.code..self.reported..)","",.)]
  med[,phe:=make.names(Phenotype.Description) %>% gsub("Cancer.code..self.reported..","SRC:",.)]
  med[,phe:=gsub("Treatment.medication.code..","SRM:",phe)]
  med[,phe:=gsub("Non.cancer.illness.code..self.reported..","SRD:",phe)]
  med[,phe:=gsub("Diagnoses...main.ICD10..","ICD10:",phe)]
  med <- med[,.(Phenotype.Code,phe)]
  med <- merge(P,med,by.x='phenotype',by.y='Phenotype.Code')[,.(code=phenotype,phe,n1=as.numeric(cases),n0=as.numeric(controls))]
  med
}

anno <- get_phenotype_annotation(".*")[!is.na(n1) & !is.na(n0) & phe != 'unclassifiable',]
## compose an annotation bar
ukbb.traits <- ukbb$trait %>% unique %>% gsub("^bb\\_","",.)
anno.bb <- anno[phe %in% ukbb.traits]
anno.bb[,c('clade','coding'):=tstrsplit(code,'_') %>% lapply(.,as.numeric)]
anno.bb[clade=='20002',top.node:='bb_disease']
anno.bb[clade=='20001',top.node:='bb_cancer']
anno.bb[clade=='20003',top.node:='bb_medications']

ukbb <- merge(ukbb,anno.bb[,.(phe,n1,n0,category=top.node)],by.x='trait',by.y='phe')
ukbb <- ukbb[,trait:=paste('bb',trait,sep='_')]


## read in roslin geneatlas results
GENEATLAS <- '/home/ob219/share/as_basis/GWAS/geneatlas/13_traits_0919_tmp'
ga <- lapply(list.files(path=GENEATLAS,pattern="*.RDS",full.names=TRUE),readRDS) %>% rbindlist
ga <- melt(ga,id.var='trait')
ga[,trait:=paste('GA',trait,sep=':')]
meta.dt <- fread("~/tmp/41588_2018_248_MOESM3_ESM.csv")
meta.dt <- meta.dt[Category=='Binary',.(ID,Description=paste('GA',make.names(Description),sep=':'),Cases,Controls=round(Cases/Sample)-Cases,prop=Sample)]
srd.idx <- grep("^selfReported",meta.dt$ID)
icd.idx <- grep("^clinical_c",meta.dt$ID)
cancer.idx <- grep("^cancer_c",meta.dt$ID)
ga.srd <- merge(ga,meta.dt[srd.idx,.(Description,n1=Cases,n0=Controls,category='geneatlas_srd')],by.x='trait',by.y='Description')
ga.icd <- merge(ga,meta.dt[icd.idx,.(Description,n1=Cases,n0=Controls,category='geneatlas_icd')],by.x='trait',by.y='Description')
ga.cancer <- merge(ga,meta.dt[cancer.idx,.(Description,n1=Cases,n0=Controls,category='geneatlas_cancer')],by.x='trait',by.y='Description')
#lefthandedness etc.
#ga.other <- merge(ga,meta.dt[-(c(srd.idx,icd.idx,cancer.idx)),.(Description,n1=Cases,n0=Controls,category='geneatlas_other')],by.x='trait',by.y='Description')
ga <- rbindlist(list(ga.srd,ga.icd,ga.cancer))
ga[,n:=n0+n1]

## get sample counts etc from previous results

pr <- readRDS('/home/ob219/share/as_basis/GWAS/RESULTS/04b_09_19_0619_summary_results.RDS')
pr <- pr[,.(trait,n0,n1,sdy,n,category)] %>% unique

## load in non ukbb

nonukbb <- readRDS('/home/ob219/share/as_basis/GWAS/RESULTS/non_uk_bb_for_fdr_13_traits_0919.RDS')
nonukbb <- data.table(trait=rownames(nonukbb),nonukbb)
nonukbb <- melt(nonukbb,id.var='trait')

nonukbb[trait=='li_ankspond',trait:='li_as']
nonukbb[trait=='myositis_myogen_ssimp',trait:='gwas_dmjdmpm_new_ssimp']
nonukbb[trait=='renton_mg',trait:='renton_mg_combined']



nonukbb<-merge(nonukbb,pr,by='trait',nonukbb.x=TRUE)
all.proj <- rbindlist(list(ukbb[,n:=n0+n1],nonukbb,ga),fill=TRUE)

## load in lmm egpa results
n0 <- 6688
sc <- list(all.egpa=534,
  anca.negative.egpa=358,
  mpo.anca.positive.egpa=159)

egpa.lmm <- lapply(names(sc),function(x){
  ## note that the labelling is wrong these are actually 13 trait basis projections
  f <- sprintf("%s_0619.RDS",x) %>% file.path('/home/ob219/share/as_basis/GWAS/lyons_egpa_lmm/projections/',.)
  readRDS(f)
}) %>% do.call('rbind',.)
egpa.lmm <- data.table(trait=rownames(egpa.lmm),egpa.lmm)
egpa.lmm[,trait:=sprintf("%s_lmm",trait)]
egpa.lmm <- melt(egpa.lmm,id.var='trait')

tpr <- lapply(names(sc),function(x){
  data.table(trait=sprintf("%s_lmm",x),n0=n0,n1=sc[[x]],sdy=NA,n=sc[[x]]+n0,category='lyons_egpa_lmm')
}) %>% rbindlist
egpa.lmm<-merge(egpa.lmm,tpr,by='trait')

all.proj <- rbindlist(list(all.proj,egpa.lmm),fill=TRUE)

## add in astle

ASTLE_DIR <- '/home/ob219/share/as_basis/GWAS/astle/13_traits_1019/'
astle <- lapply(list.files(path=ASTLE_DIR,pattern="*.RDS",full.names=TRUE),readRDS) %>% do.call('rbind',.)
astle <- data.table(trait=rownames(astle),astle)
astle <- melt(astle,id.var='trait')

## for analysis we will stick to just the 13 main blood types mentioned in Astle et al.
keep <- c('pdw','mpv','plt','irf','ret','rdw','hct','mch','mono','baso','eo','neut','lymph')
astle <- astle[trait %in% keep,]
DATA_DIR <- '/home/ob219/share/Data/GWAS-summary/blood-ukbiobank-2016-12-12'
afiles <- list.files(path=DATA_DIR,pattern="*.gz$",full.names=FALSE)
#astle_samples <- data.table(trait = gsub("(.*)\\_build37\\_[0-9]+\\_20161212.tsv.gz","\\1",afiles),
#ss = as.numeric(gsub("(.*)\\_build37\\_([0-9]+)\\_20161212.tsv.gz","\\2",afiles)))
#astle_samples[,c('n0','n1','ss'):=list(round(ss/2),round(ss/2),NULL)]
## here we assume that the trait outcome has been standardised such that sdY==1.
astle_samples <- data.table(trait = gsub("(.*)\\_build37\\_[0-9]+\\_20161212.tsv.gz","\\1",afiles),
n0 = as.numeric(gsub("(.*)\\_build37\\_([0-9]+)\\_20161212.tsv.gz","\\2",afiles)),
n1 = 0, sdy=1)
astle <- merge(astle,astle_samples)
astle[,category:='astle_blood']
astle[,n:=n0]

all.proj <- rbind(all.proj,astle,fill=TRUE)

## add in tian

tian <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/tian_infectious_disease_13_traits_0919.RDS")
## get sample counts and merge
tian_samples <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/sample_size.RDS")
tian <- melt(tian,id.var='trait')
tian <- merge(tian,tian_samples,by='trait')
tian[,category:='tian_infectious_disease']
tian[,n:=n0 + n1]

all.proj <- rbind(all.proj,tian,fill=TRUE)


## renton mg imputed
n_controls <- 1977
cases <- list(YoungOnset=235,LateOnset=737,Overall=972)
DATA_DIR <- '/home/ob219/share/as_basis/GWAS/renton_mg/'
renton_ssimp <- lapply(names(cases),function(trait){
  file <- sprintf("%s_ssimp_13_traits_0919.RDS",trait) %>% file.path(DATA_DIR,.)
  tmp <- readRDS(file)
  tmp <- data.table(trait=trait,tmp)
  tmp <- tmp[,c('n1','n0'):=list(cases[[trait]],n_controls)]
}) %>% rbindlist()
renton_ssimp <- melt(renton_ssimp,id.vars=c('trait','n0','n1'))
setcolorder(renton_ssimp,c('trait','variable','value','n0','n1'))
renton_ssimp[,n:=n0+n1]
renton_ssimp[,category:='renton_mg_ssimp']
renton_ssimp[,trait:=sprintf("%s_ssimp",trait)]

all.proj <- rbind(all.proj,renton_ssimp,fill=TRUE)

## nmo imputed
## arterido psa imputed
## myositis imputed


VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_gwas_13_traits_0919.RDS'
var.DT <- readRDS(VARIANCE_FILE)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,value:=(value-control.loading)]
basis.DT[,Z:=sign(value)]
basis.DT[,p.adj:=1]
basis.DT[,control.loading:=NULL]


### new method for computing variance using seb
variance.dir <- '/home/ob219/share/as_basis/GWAS/seb_proj_var_13_traits_0919'
files <- list.files(path=variance.dir,pattern='*.RDS',full.names=TRUE)
all.var <- lapply(files,function(f){
  tmp <- readRDS(f)
  ta<- basename(f) %>% gsub("\\.RDS","",.)
  tmp.dt <- data.table(pc=names(tmp),trait=ta,seb.var=tmp)
}) %>% rbindlist

all.var[,trait:=gsub("UKBB_NEALE:","bb_",trait)]
all.var[,trait:=gsub("TIAN:","",trait)]
all.var[trait=='renton_mg',trait:='renton_mg_combined']
all.var[trait=='li_ankspond',trait:='li_as']
all.var[,trait:=gsub("^ASTLE:","",trait)]
setnames(all.var,'pc','variable')

## compute the variance of a projection
all.DT <- merge(all.proj,var.DT,by.x='variable',by.y='pc')
all.DT <- merge(all.DT,control.DT,by.x='variable',by.y='PC')
all.DT[is.na(sdy),variance:=((log(n)-(log(n1) + log(n-n1)))+ log(mfactor)) %>% exp]
all.DT[!is.na(sdy),variance:=(log(sdy^2/n) + log(mfactor)) %>% exp]

##some will be missing need to add these - Limy, mahajan_t2d and all Astle.
all.DT<-merge(all.DT,all.var,by=c('variable','trait'),all.x=TRUE)
## need Chris' basis-sparse-0.999.RData file to do this

## add in the imputation variance for myogen
myogen.ssimp.var <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_empirical_variances_13_traits_0919.RDS')
ssimp.id <- which(grepl("_ssimp$",all.DT$trait) & all.DT$category == 'myogen')
keep <- all.DT[ssimp.id,]
keep <- merge(keep,myogen.ssimp.var[,.(trait,variable,new.variance=value)],by=c('trait','variable'))
keep[,variance:=new.variance]
keep$new.variance <- NULL
all.DT <- rbind(all.DT[-ssimp.id,],keep)

## add in the imputation variance for arterido
psa_arterido.ssimp.var <- readRDS('/home/ob219/share/as_basis/GWAS/psa_aterido/psa_aterido_empirical_variances_13_traits_0919.RDS')
ssimp.id <- which(grepl("_ssimp$",all.DT$trait) & all.DT$category == 'psa_aterido')
keep <- all.DT[ssimp.id,]
keep <- merge(keep,psa_arterido.ssimp.var[,.(trait,variable,new.variance=value)],by=c('trait','variable'))
keep[,variance:=new.variance]
keep$new.variance <- NULL
all.DT <- rbind(all.DT[-ssimp.id,],keep)

## add in the imputation variance for nmo
nmo.ssimp.var <- readRDS('/home/ob219/share/as_basis/GWAS/nmo/nmo_empirical_variances_13_traits_0919.RDS')
ssimp.id <- which(grepl("_ssimp$",all.DT$trait) & all.DT$category == 'estrada_NMO')
keep <- all.DT[ssimp.id,]
keep <- merge(keep,nmo.ssimp.var[,.(trait,variable,new.variance=value)],by=c('trait','variable'))
keep[,variance:=new.variance]
keep$new.variance <- NULL
all.DT <- rbind(all.DT[-ssimp.id,],keep)


## add in control loading
all.DT[,Z:=(value-control.loading)/sqrt(variance)]
all.DT[,Z.seb:=(value-control.loading)/sqrt(seb.var)]
all.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
all.DT[,p.value.seb:=pnorm(abs(Z.seb),lower.tail=FALSE) * 2]
all.DT[,delta:=value-control.loading]
all.DT <- all.DT[!trait %in% c('cousminer_lada','IgA_nephropathy'),]
## correct imputed variances
saveRDS(all.DT,'/home/ob219/share/as_basis/GWAS/RESULTS/10_10_13_traits_0919_summary_results_with_seb_variance.RDS')
## obtain a summary
all.DT[,.(trait,category),by=c('trait','category')][,list(count=.N),by='category'][order(count),]
if(FALSE){
  rm.categories <- c("bb_medications","ad-pid","tachmazidou_osteo",
      "mahajan_t2d","Ig.titre","rhodes_pah",'astle','wong_aav','geneatlas_srd',
      'geneatlas_cancer')
primary.DT <- all.DT[!category %in% rm.categories]
primary.DT[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']
primary.DT[,p.adj.seb:=p.adjust(p.value.seb,method="fdr"),by='variable']
saveRDS(primary.DT,'/home/ob219/share/as_basis/GWAS/RESULTS/13_traits_0919_primary_results.RDS')
}
