## overall file to collate all projection results

## load in UKBB results

BB_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/0619'
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


## read in ICD projections
## Neale ICD projections are currently broken therefore we use GeneATLAS instead
if(FALSE){
  ICD_BB_DIR <- '/home/ob219/share/as_basis/GWAS/bb_projections/0619_ICD'
  icd <- lapply(list.files(path=ICD_BB_DIR,pattern="*.RDS",full.names=TRUE),readRDS) %>% rbindlist
  icd <- melt(icd,id.var='trait')
  icd[,trait:=paste('ICD10',trait,sep=':')]
  icd.traits <- icd$trait %>% unique
  anno.icd <- anno[phe %in% icd.traits]
  icd <- merge(icd,anno.icd[,.(phe,n1,n0,category='bb_icd10')],by.x='trait',by.y='phe')
}

## read in roslin geneatlas results

GENEATLAS <- '/home/ob219/share/as_basis/GWAS/geneatlas/0619'
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


##remove medication for the time being
#ukbb <- ukbb[!is.na(category),]

## load in jia sub types

jia <- readRDS("/home/ob219/share/as_basis/GWAS/jia_projections/summary/jia_0619.RDS")
sample.DT <- fread('/home/ob219/share/Data/GWAS/jia-mar-2019/summary-stats-samplecount-mar2019.csv')
setnames(sample.DT,'n','n1')
sub.DT <- data.table(idx=0:9,subtype=c('case','sys','PO','EO','RFneg','RFpos','ERA','PsA','undiff','missing'))
samp.DT <- merge(sub.DT,sample.DT[,.(ilar_pheno,n1)],by.x='idx',by.y='ilar_pheno')
samp.DT <- rbind(samp.DT,data.table(idx=0,subtype='case',n1=sum(samp.DT$n1)))
samp.DT[,n0:=9196]
jia <- melt(jia,id.var='trait')
samp.DT[,subtype:=sprintf("jia_%s_19",subtype)]
jia <- merge(jia,samp.DT[,.(subtype,n1,n0)],by.x='trait',by.y='subtype')
jia[,category:='bowes_jia_2019']

## do methotrexate study
mtx <- readRDS("/home/ob219/share/as_basis/GWAS/mtx/mtx_matura_0619.RDS")
mtx <- data.table(trait=rownames(mtx),mtx)
mtx<- melt(mtx,id.var='trait')
mtx[,c('n0','n1'):=list(1424,0)]
mtx[,category:='taylor_mtx']
mtx[,sdy:=0.9215145]

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

tian <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/tian_infectious_disease_0619.RDS")
## get sample counts and merge
tian_samples <- readRDS("/home/ob219/share/as_basis/GWAS/tian_projections/sample_size.RDS")
tian <- melt(tian,id.var='trait')
tian <- merge(tian,tian_samples,by='trait')
tian[,category:='tian_infectious_disease']

## load in ferreira_asthma

ferriera <- readRDS("/home/ob219/share/as_basis/GWAS/ferreira_projections/ferreira_asthma_0619.RDS")
ferreira_samples <- readRDS("/home/ob219/share/as_basis/GWAS/ferreira_projections/sample_size.RDS")
ferriera <- melt(ferriera,id.var='trait')
ferriera <- merge(ferriera,ferreira_samples,by='trait')
ferriera[,category:='ferreira_asthma']


## load in astle data

## first get sample sizes - note that because of assumptions in converting linear coefficient
## to OR the cases and control size are equal to the median/mean of the sample size i.e round(ss/2)
ASTLE_DIR <- '/home/ob219/share/as_basis/GWAS/astle/0619/'
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

## load in ig titres
## this is unpublished and based on old method for converting to OR scale so leave for now
if(FALSE){
  eff <- readRDS('/home/ob219/share/as_basis/GWAS/effrosyni_ig/effrosyni_ig_0619.RDS')
  eff <- melt(eff,id.var='trait')
  eff[,c('n0','n1'):=list(3969,3969)]
  eff[,category:='Ig.titre']
}
## load in CD prognosis
cd.prog <- readRDS('/home/ob219/share/as_basis/GWAS/cd_prognosis/cd_prognosis_0619.RDS')
cd.prog <- data.table(trait=rownames(cd.prog),cd.prog)
cd.prog <- melt(cd.prog,id.var='trait')
cd.prog[,c('n0','n1'):=list(389 + 583,669 + 1093)]
cd.prog[,category:='lee_CD_prognosis']

## load in bp,adhd, and scz

files <- list.files(path="/home/ob219/share/as_basis/GWAS/psych/",pattern="*0619.RDS",full.names=TRUE)

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
myogen <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis_0619.RDS')
myogen <- data.table(trait=rownames(myogen),myogen)
myogen <- melt(myogen,id.var='trait')
myogen[,c('n0','n1'):=list(4724,1711)]
myogen[,category:='myogen']



## jdm
#not imputed

n_controls <- 4724
cases <- list(dm=705,jdm=473,pm=533,myogen=1711)
all.myo <- lapply(names(cases),function(trait){
  in.file <- file.path('/home/ob219/share/as_basis/GWAS/myogen_myositis/',sprintf("%s_myositis_0619.RDS",trait))
  myogen_ss <- readRDS(in.file)
  myogen_ss <- data.table(trait=rownames(myogen_ss),myogen_ss)
  myogen_ss <- melt(myogen_ss,id.var='trait')
  n_cases <- cases[[trait]]
  myogen_ss[,c('n0','n1'):=list(n_controls,n_cases)]
  myogen_ss[,category:='myogen']
}) %>% rbindlist

## imputed
n_controls <- 4724
cases <- list(dm=705,jdm=473,pm=533,dmjdmpm=1711)
all.myo.ssimp <- lapply(names(cases),function(trait){
  in.file <- file.path('/home/ob219/share/as_basis/GWAS/myogen_myositis/',sprintf("gwas_%s_new_ssimp_0619.RDS",trait))
  myogen_ss <- readRDS(in.file)
  myogen_ss <- data.table(trait=rownames(myogen_ss),myogen_ss)
  myogen_ss <- melt(myogen_ss,id.var='trait')
  n_cases <- cases[[trait]]
  myogen_ss[,c('n0','n1'):=list(n_controls,n_cases)]
  myogen_ss[,category:='myogen']
}) %>% rbindlist

myogen <- rbindlist(list(all.myo,all.myo.ssimp))


## myogen_myositis
# myogen.flip <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis_flip.RDS')
# myogen.flip <- data.table(trait=rownames(myogen.flip),myogen.flip)
# myogen.flip <- melt(myogen.flip,id.var='trait')
# myogen.flip[,c('n0','n1'):=list(4726,1818)]
# myogen.flip[,category:='myogen']

## NMO - Neuromyelitis Optica
n_controls <- 1244
cases <- list(IgPos=132,IgNeg=83,combined=215)
all.nmo <- lapply(names(cases),function(trait){
  in.file <- file.path('/home/ob219/share/as_basis/GWAS/nmo',sprintf("NMO_%s_0619.RDS",trait))
  nmo_ni <- readRDS(in.file)
  nmo_ni <- data.table(trait=rownames(nmo_ni),nmo_ni)
  nmo_ni <- melt(nmo_ni,id.var='trait')
  n_cases <- cases[[trait]]
  nmo_ni[,c('n0','n1'):=list(n_controls,n_cases)]
  nmo_ni[,category:='estrada_NMO']
}) %>% rbindlist


all.nmo.ssimp <- lapply(names(cases),function(trait){
  in.file <- file.path('/home/ob219/share/as_basis/GWAS/nmo',sprintf("NMO_%s_ssimp_0619.RDS",trait))
  nmo_ni <- readRDS(in.file)
  nmo_ni <- data.table(trait=rownames(nmo_ni),nmo_ni)
  nmo_ni <- melt(nmo_ni,id.var='trait')
  n_cases <- cases[[trait]]
  nmo_ni[,c('n0','n1'):=list(n_controls,n_cases)]
  nmo_ni[,category:='estrada_NMO']
}) %>% rbindlist

nmo <- rbind(all.nmo,all.nmo.ssimp)


# files <- list.files(path="/home/ob219/share/as_basis/GWAS/nmo/",pattern="*0619.RDS",full.names=TRUE)
# nmo <- lapply(files,readRDS) %>% do.call('rbind',.)
# #nmo <- readRDS('/home/ob219/share/as_basis/GWAS/nmo/nmo.RDS')
# nmo <- data.table(trait=rownames(nmo),nmo)
# nmo <- melt(nmo,id.var='trait')
# nmo[trait=='NMO_combined',c('n0','n1'):=list(1244,215)]
# nmo[trait=='NMO_IgGPos',c('n0','n1'):=list(1244,66+66)]
# nmo[trait=='NMO_IgGNeg',c('n0','n1'):=list(1244,20+63)]
# nmo[,category:='estrada_NMO']



## load in egpa vasc

files <- list.files(path="/home/ob219/share/as_basis/GWAS/lyons_egpa/projections",pattern="*0619.RDS",full.names=TRUE)
egpa <- lapply(files,readRDS) %>% do.call('rbind',.)
egpa <- data.table(trait=rownames(egpa),egpa)
egpa <- melt(egpa,id.var='trait')
egpa.samp <- readRDS('/home/ob219/share/as_basis/GWAS/lyons_egpa/sample.RDS')
egpa <- merge(egpa,egpa.samp,by='trait')
egpa[,category:='lyons_egpa']

## MyastheniaGravis_Renton_JAMA_Neurol_2015

mg <- readRDS("/home/ob219/share/as_basis/GWAS/renton_mg/projections/renton_mg_0619.RDS")
mg <- data.table(trait=rownames(mg),mg)
mg <- melt(mg,id.var='trait')
mg[trait=='renton_mg',c('n0','n1'):=c(1977,972)]
mg[trait=='renton_mg_late',c('n0','n1'):=c(1977,737)]
mg[trait=='renton_mg_early',c('n0','n1'):=c(1977,235)]
mg[trait=='renton_mg',trait:='renton_mg_combined']
mg[,category:='renton_mg']

files <- list.files(path="/home/ob219/share/as_basis/GWAS/lyons_egpa/projections",pattern="*0619.RDS",full.names=TRUE)
egpa <- lapply(files,readRDS) %>% do.call('rbind',.)
egpa <- data.table(trait=rownames(egpa),egpa)
egpa <- melt(egpa,id.var='trait')
egpa.samp <- readRDS('/home/ob219/share/as_basis/GWAS/lyons_egpa/sample.RDS')
egpa <- merge(egpa,egpa.samp,by='trait')
egpa[,category:='lyons_egpa']


## load in aav (vasculitis) from Limy Wong

aav <- readRDS("/home/ob219/share/as_basis/GWAS/wong_aav/projections/aav_0619.RDS")
aav <- melt(aav,id.var='trait')
aav.samp <- fread("/home/ob219/share/Data/GWAS-summary/aav_limy_wong/aav_sample_size.txt")
aav <- merge(aav,aav.samp[,.(label,n0,n1)],by.x='trait',by.y='label')
aav[,category:='wong_aav']

## load in osteo data

ost <- readRDS("/home/ob219/share/as_basis/GWAS/tachmazidou_osteo/projections/oa_0619.RDS")
ost <- melt(ost,id.var='trait')
ost.samp <- fread("/home/ob219/share/Data/GWAS-summary/tachmazidou_osteo/osteo_sample_size.txt")
ost <- merge(ost,ost.samp[,.(label,n0,n1)],by.x='trait',by.y='label')
ost[,category:='tachmazidou_osteo']

## load in abdef

abdef <- readRDS('/home/ob219/share/as_basis/GWAS/abdef/abdef_0619.RDS')
abdef <- data.table(trait=rownames(abdef),abdef)
abdef <- melt(abdef,id.var='trait')
abdef[,c('n0','n1'):=list(9225,733)]
abdef[,category:='ad-pid']

## LEAVE OUT FOR NOW
if(FALSE){
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
}

## unpublished psa bowes
psa <- readRDS("/home/ob219/share/as_basis/GWAS/psa_projections/summary/bowes_psa_0619.RDS")
psa <- melt(psa,id.var='trait')
psa[,c('n0','n1','category'):=list(4596,1805,'bowes_psa')]

## psa aterido one north american and one spanish cohort

n_controls <- 1454
cases <- list(span=744,na=1430)
all.psa_aterido.ssimp <- lapply(names(cases),function(trait){
  in.file <- file.path('/home/ob219/share/as_basis/GWAS/psa_aterido',sprintf("%s_psa_ssimp_0619.RDS",trait))
  psa_art <- readRDS(in.file)
  psa_art <- data.table(trait=rownames(psa_art),psa_art)
  psa_art <- melt(psa_art,id.var='trait')
  n_cases <- cases[[trait]]
  psa_art[,c('n0','n1'):=list(n_controls,n_cases)]
  psa_art[,category:='psa_aterido']
}) %>% rbindlist
all.psa_aterido.ssimp[trait=='na_psa_ssimp',n0:=1417]


psa_aterido <- readRDS('/home/ob219/share/as_basis/GWAS/psa_aterido/psa_aterido_0619.RDS')
psa_aterido <- melt(psa_aterido,id.var='trait') %>% data.table
setnames(psa_aterido,c('trait','variable','value'))
psa_aterido[trait=='na_psa',c('n0','n1','category'):=list(1417,1430,'psa_aterido')]
psa_aterido[trait=='span_psa',c('n0','n1','category'):=list(1454,744,'psa_aterido')]

psa_aterido <- rbind(all.psa_aterido.ssimp,psa_aterido)


## this is from Mark Toshner - I was unable to find evidence of immune-mediated component
## therefore leave
if(FALSE){
  pah <- readRDS("/home/ob219/share/as_basis/GWAS/liley_pah/projections/pah_0619.RDS")
  pah <- melt(pah,id.var='trait') %>% data.table
  setnames(pah,c('trait','variable','value'))
  pah[,c('n0','n1','category'):=list(5045,847,'rhodes_pah')]
  pah.np <- readRDS("/home/ob219/share/as_basis/GWAS/liley_pah/projections/pah_no_pid_0619.RDS")
  pah.np <- melt(pah.np,id.var='trait') %>% data.table
  setnames(pah.np,c('trait','variable','value'))
  pah.np[,c('n0','n1','category'):=list(4243,848,'rhodes_pah')]
  pah <- rbind(pah,pah.np)
}


#IgA_nephropathy
iga <- readRDS("/home/ob219/share/as_basis/GWAS/IgA_nephropathy/IgA_nephropathy_0619.RDS")
iga <- melt(iga,id.var='trait') %>% data.table
setnames(iga,c('trait','variable','value'))
iga[,c('n0','n1','category'):=list(3952,2747,'kiryluk_iga_neph')]

## ankylosing_spondylitis
as <- readRDS("/home/ob219/share/as_basis/GWAS/ank_spond/ank_spond_0619.RDS")
as <- melt(as,id.var='trait') %>% data.table
setnames(as,c('trait','variable','value'))
as[,c('n0','n1','category'):=list(1644,4880,'brown_as')]

## Iranain study
li_as <- readRDS('/home/ob219/share/as_basis/GWAS/li_ankspond/li_ankspond_0619.RDS')
li_as <- melt(li_as,id.var='trait') %>% data.table
setnames(li_as,c('trait','variable','value'))
li_as[,c('n0','n1','category'):=list(1841,1480,'li_as')]

## birdshot_retinopathy
bs <- readRDS("/home/ob219/share/as_basis/GWAS/birdshot_retinopathy/birdshot_retinopathy_0619.RDS")
bs <- melt(bs,id.var='trait') %>% data.table
setnames(bs,c('trait','variable','value'))
bs[,c('n0','n1','category'):=list(693,117,'kuiper_bs')]

## Latent autoimmune diabetes in adults
lada <- readRDS("/home/ob219/share/as_basis/GWAS/cousminer_lada/cousminer_lada_0619.RDS")
lada <- melt(lada,id.var='trait') %>% data.table
setnames(lada,c('trait','variable','value'))
lada[,c('n0','n1','category'):=list(5947,2634,'cousminer_lada')]

## T2D to compare with LADA and T1D
t2d_mahajan <- readRDS('/home/ob219/share/as_basis/GWAS/mahajan_t2d/mahajan_t2d_0619.RDS')
t2d_mahajan <- melt(t2d_mahajan,id.var='trait') %>% data.table
setnames(t2d_mahajan,c('trait','variable','value'))
t2d_mahajan[,c('n0','n1','category'):=list(824006,74124,'mahajan_t2d')]

## jia uveitis

hasnoot_uveitis <- readRDS('/home/ob219/share/as_basis/GWAS/hasnoot_uveitis/hasnoot_uveitis_0619.RDS')
hasnoot_uveitis <- melt(hasnoot_uveitis,id.var='trait') %>% data.table
setnames(hasnoot_uveitis,c('trait','variable','value'))
hasnoot_uveitis[,c('n0','n1','category'):=list(330,192,'hasnoot_uveitis_jia')]

## cytokines

cyt <- readRDS('/home/ob219/share/as_basis/GWAS/ahola-olli_cytokine/projections/ck_0619.RDS')
cyt <- melt(cyt,id.var='trait') %>% data.table
setnames(cyt,c('trait','variable','value'))
## read in sample sizes from file
ss.cyt <- fread("/home/ob219/share/Data/GWAS-summary/cytokine_gwas/sample_size.txt")
ss.cyt[,label:=paste('CK',label,sep=':')]
cyt<-merge(cyt,ss.cyt[,.(trait=label,n0=n,n1=0,sdy=1)],by='trait',all.x=TRUE)
cyt[,category:='ahola-olli_cytokine']

## roederer immunophenotype data

ROE_DIR <- '/home/ob219/share/as_basis/GWAS/roederer/'
ROE_SS_FILE <- '/home/ob219/rds/hpc-work/roederer/twinr-ftp.kcl.ac.uk/ImmuneCellScience/2-GWASResults/sample_size.txt'
ss.roe <- fread(ROE_SS_FILE)
setnames(ss.roe,c('trait','mean','max'))
roe.files <- list.files(path=ROE_DIR,pattern="*.RDS",full.names=TRUE)
roe <- lapply(roe.files,function(f){
  dat <- readRDS(f)
  data.table(trait=rownames(dat),dat)
}) %>% rbindlist
roe <- melt(roe,id.var='trait') %>% data.table
roe<-merge(roe,ss.roe[,.(trait,n0=max,n1=0)],by='trait',all.x=TRUE)
roe.sdy <- readRDS('/home/ob219/rds/hpc-work/roederer/twinr-ftp.kcl.ac.uk/ImmuneCellScience/2-GWASResults/sdy.estimate.RDS')
roe <- merge(roe,roe.sdy[,.(trait,sdy=sdY)],by='trait')
roe[,category:='roederer_immunophenotypes']

## for some reason we swapped from n1 to n0 halfway through to n0 n1
## we need to fix otherwise everything gets swapped and case
## sizes become control sizes

ncolorder<-c('trait','variable','value','n0','n1')
setcolorder(ferriera,ncolorder)
setcolorder(tian,ncolorder)
setcolorder(jia,ncolorder)
setcolorder(ukbb,ncolorder)

## ommit for the time being from analysis
if(FALSE){
  ## add eqtlgen
  EQTLGEN_DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/eqtlgen_projections_significant_only_0619/'
  fs <- list.files(path=EQTLGEN_DIR,pattern="*.RDS",full.names=TRUE)
  res.DT <- lapply(fs,readRDS) %>% rbindlist
  ## define a Z score for loading to see if any are significant across pc's
  eqtlgen_fdr <- melt(res.DT,id.vars='trait')
  eqtlgen_fdr[,c('n0','n1','sdy','category'):=list(31684,0,1,'eqtlgen_fdr0.05')]

  PQTL_DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/sun_pqtl/fdr_0.05_by_chr_filtered/'
  fs <- list.files(path=PQTL_DIR,pattern="*.RDS",full.names=TRUE)
  res.DT <- lapply(fs,function(x){
    message(x)
	   mat<-readRDS(x)
	    data.table(trait=rownames(mat),mat)
  }) %>% rbindlist
  ## define a Z score for loading to see if any are significant across pc's
  pqtl_fdr <- melt(res.DT,id.vars='trait')
  pqtl_fdr[,c('n0','n1','sdy','category'):=list(3301,0,1,'sun_pqtl_fdr0.05')]
}

all.proj <- list(
  ferriera=ferriera,
  tian=tian,
  jia=jia,
  ukbb=ukbb,
  #icd=icd, - did not work not sure why !
  astle=astle,
  #eff=eff,
  cd.prog=cd.prog,
  psy=psy,
  myogen=myogen,
  nmo=nmo,
  egpa=egpa,
  abdef=abdef,
  ost=ost,
  #liley=t1d,
  aav=aav,
  psa=psa,
  #pah=pah,
  mtx=mtx,
  iga=iga,
  as=as,
  bs=bs,
  lada=lada,
  li_as=li_as,
  ga=ga,
  t2d_mahajan=t2d_mahajan,
  psa_aterido=psa_aterido,
  hasnoot_uveitis=hasnoot_uveitis,
  cyt,
  roe,
  mg=mg
  #pqtl_fdr=pqtl_fdr,
  #eqtlgen_fdr=eqtlgen_fdr
) %>% rbindlist(.,fill=TRUE)
all.proj[,n:=n1+n0]

## load in basis and variance
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_av_0619.RDS'
var.DT <- readRDS(VARIANCE_FILE)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
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
all.DT[is.na(sdy),variance:=((log(n)-(log(n1) + log(n-n1)))+ log(mfactor)) %>% exp]
all.DT[!is.na(sdy),variance:=(log(sdy^2/n) + log(mfactor)) %>% exp]

## add in the imputation variance for myogen
myogen.ssimp.var <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_empirical_variances_0619.RDS')
ssimp.id <- which(grepl("_ssimp$",all.DT$trait) & all.DT$category == 'myogen')
keep <- all.DT[ssimp.id,]
keep <- merge(keep,myogen.ssimp.var[,.(trait,variable,new.variance=value)],by=c('trait','variable'))
keep[,variance:=new.variance]
keep$new.variance <- NULL
all.DT <- rbind(all.DT[-ssimp.id,],keep)

## add in the imputation variance for arterido
psa_arterido.ssimp.var <- readRDS('/home/ob219/share/as_basis/GWAS/psa_aterido/psa_aterido_empirical_variances_0619.RDS')
ssimp.id <- which(grepl("_ssimp$",all.DT$trait) & all.DT$category == 'psa_aterido')
keep <- all.DT[ssimp.id,]
keep <- merge(keep,psa_arterido.ssimp.var[,.(trait,variable,new.variance=value)],by=c('trait','variable'))
keep[,variance:=new.variance]
keep$new.variance <- NULL
all.DT <- rbind(all.DT[-ssimp.id,],keep)

## add in the imputation variance for nmo
nmo.ssimp.var <- readRDS('/home/ob219/share/as_basis/GWAS/nmo/nmo_empirical_variances_0619.RDS')
ssimp.id <- which(grepl("_ssimp$",all.DT$trait) & all.DT$category == 'estrada_NMO')
keep <- all.DT[ssimp.id,]
keep <- merge(keep,nmo.ssimp.var[,.(trait,variable,new.variance=value)],by=c('trait','variable'))
keep[,variance:=new.variance]
keep$new.variance <- NULL
all.DT <- rbind(all.DT[-ssimp.id,],keep)


## add in control loading
all.DT[,Z:=(value-control.loading)/sqrt(variance)]
all.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
all.DT[,delta:=value-control.loading]
## correct imputed variances
saveRDS(all.DT,'/home/ob219/share/as_basis/GWAS/RESULTS/04_09_19_0619_summary_results.RDS')

## primary results
rm.categories <- c("bb_medications","ad-pid","tachmazidou_osteo",
"mahajan_t2d","Ig.titre","rhodes_pah",'astle','wong_aav','geneatlas_srd','geneatlas_cancer')

primary.DT <- all.DT[!category %in% rm.categories]
primary.DT[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']
saveRDS(primary.DT,'/home/ob219/share/as_basis/GWAS/RESULTS/04_09_19_0619_primary_results.RDS')


if(FALSE){
  all.DT <- readRDS("~/share/as_basis/GWAS/RESULTS/23_07_19_0619_summary_results.RDS")
  qqnorm(all.DT[category=='eqtlgen_fdr0.05',]$Z,main='eqtlgen_fdr0.05')
  qqline(all.DT[category=='eqtlgen_fdr0.05',]$Z,col='red',lty=2)
  qqnorm(all.DT[category=='sun_pqtl_fdr0.05',]$Z,main='sun_pqtl_fdr0.05')
  qqline(all.DT[category=='sun_pqtl_fdr0.05',]$Z,col='red',lty=2)
  qqnorm(all.DT[category=='bb_disease',]$Z,main='Neale SRD')
  qqline(all.DT[category=='bb_disease',]$Z,col='red',lty=2)
  qqnorm(all.DT[category=='astle_blood',]$Z,main='Neale SRD')
  qqline(all.DT[category=='astle_blood',]$Z,col='red',lty=2)
  foo <- all.DT[category=='sun_pqtl_fdr0.05',][,list(mvalue=mean(value),control=mean(control.loading)),by='variable']
  ggplot(all.DT[category=='sun_pqtl_fdr0.05',],aes(x=variable,y=value)) + geom_boxplot()
  all.DT[category=='sun_pqtl_fdr0.05' & variable=='PC1',][order(Z),] %>% head
}
