library(xlsx)
OUT_DIR <- '/home/ob219/share/as_basis/supp_tables/'

### this code creates supplementary tables

## supp table 1 - basis traits
TRAIT_MANIFEST <- '/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
trait.manifest <- fread(TRAIT_MANIFEST)[basis_trait==1,]
trait.manifest[,c('dtr','first_author'):=tstrsplit(disease,'_')]
trait.manifest <- trait.manifest[order(cases+controls),]
basis.out <- trait.manifest[,.(Trait=dtr,`First Author`=first_author,Reference=pmid,N0=controls,N1=cases)]
write.xlsx(basis.out,file=file.path(OUT_DIR,'supplementary_tables_1_4.xlsx'),sheet='Basis',row.names=FALSE)


## supp table 2 - Neale traits

RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/02_08_19_0619_primary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)[,.(trait,category,n0,n1,sdy)] %>% unique
keep.cat <- c('bb_disease','bb_medications','bb_cancer')
neale.out <- res.DT[category %in% keep.cat,]
neale.out[,c('pre','suf'):=tstrsplit(trait,':')]
neale.out[,trait:=suf]
neale.out <- neale.out[order(n1),.(Trait=trait,N0=n0,N1=n1),by=category]

write.xlsx(neale.out,file=file.path(OUT_DIR,'supplementary_tables_1_4.xlsx'),sheet='UKBB Neale',row.names=FALSE,append=TRUE)

## supp table 3 - GeneAtlas traits - these are the ICDs at the moment but could add in srd if required for paper

ga.out <- res.DT[category=='geneatlas_icd',]
ga.out[,c('pre','suf'):=tstrsplit(trait,':')]
ga.out[,trait:=suf]
ga.out <- ga.out[order(n1),.(Trait=trait,N0=n0,N1=n1)]
write.xlsx(ga.out,file=file.path(OUT_DIR,'supplementary_tables_1_4.xlsx'),sheet='UKBB HESS GeneATLAS',row.names=FALSE,append=TRUE)

## supp table 4 All other projected traits and studies

rare.out <- res.DT[!category %in% c(keep.cat,'geneatlas_icd'),]

fs <- list.files(path='/home/ob219/share/as_basis/GWAS/for_fdr',pattern='*.RDS',full.names=TRUE)
miss.dt <- lapply(fs,function(f){
  missing <- readRDS(f)[is.na(p.value),] %>% nrow
  data.table(trait=basename(f) %>% gsub("\\_source.RDS","",.),missing)
}) %>% rbindlist

miss.dt[trait=='li_ankspond',trait:='li_as']

rare.out <- merge(rare.out,miss.dt,by.x='trait',by.y='trait',all.x=TRUE)
## tian and ferreira and astle have no SNPs missing I checked
rare.out[category %in% c('astle_blood','tian_infectious_disease'),missing:=0]



rare.out[category=='ahola-olli_cytokine',trait:=gsub("CK:","",trait)]
rare.out[,c('fa','tr'):=tstrsplit(category,'_')]
rare.out <- rare.out[,.(Trait=trait,`First Author`=fa,Reference='unpublished',N0=n0,N1=n1,sdY=sdy,missing)]


rare.out[Trait=='jia_ERA_19',Trait:='enthesitis_related_jia']
rare.out[Trait=='jia_EO_19',Trait:='extended_oligo_jia']
rare.out[Trait=='jia_PO_19',Trait:='persistent_oligo_jia']
rare.out[Trait=='jia_PsA_19',Trait:='psoriatic_jia']
rare.out[Trait=='jia_RFneg_19',Trait:='polyoligo_rf-_jia']
rare.out[Trait=='jia_RFpos_19',Trait:='polyoligo_rf+_jia']
rare.out[Trait=='jia_undiff_19',Trait:='undifferentiated_jia']
rare.out[Trait=='jia_sys_19',Trait:='systemic_jia']
rare.out[Trait=='jia_case_19',Trait:='combined_jia']
rare.out[grep('jia$',Trait),Reference:='unpublished']
rare.out[Trait=='CD_prognosis',Trait:='crohns_disease_prognosis']
rare.out[Trait=='crohns_disease_prognosis',Reference:='28067912']
rare.out[Trait=='ADHD',Trait:='attention_deficit_hyperactivity_disorder']
rare.out[Trait=='attention_deficit_hyperactivity_disorder',Reference:='29325848']
rare.out[Reference=='29325848',`First Author`:='martin']
rare.out[Trait=='BIP',Trait:='bipolar_disorder']
rare.out[Trait=='bipolar_disorder',Reference:='29906448']
rare.out[Trait=='SCZ',Trait:='schizophrenia']
rare.out[Trait=='schizophrenia',Reference:='29906448']
rare.out[Reference=='29906448',`First Author`:='psyc_consortium']
rare.out[Trait=='myositis_myogen',Trait:='combined_myositis']
rare.out[Trait=='jdm_myogen',Trait:='juvenile_dermatomyositis']
rare.out[Trait=='dm_myogen',Trait:='dermatomyositis']
rare.out[Trait=='pm_myogen',Trait:='polymyositis']
rare.out[Trait=='NMO_combined',Trait:='combined_neuromyelitis_optica']
rare.out[Trait=='NMO_IgGNeg',Trait:='igg-_neuromyelitis_optica']
rare.out[Trait=='NMO_IgGPos',Trait:='igg+_neuromyelitis_optica']
rare.out[grep('neuromyelitis_optica$',Trait),Reference:='29769526']
rare.out[Trait=='anca_Neg',Trait:='anca-_egpa']
rare.out[Trait=='egpa',Trait:='combined_egpa']
rare.out[Trait=='mpo',Trait:='mpo+_egpa']
rare.out[grep('egpa$',Trait),Reference:='https://doi.org/10.1101/491837']
rare.out[Trait=='mpo_Pos',Trait:='mpo+_aav']
rare.out[Trait=='mpo+_aav',Reference:='22808956']
## get rid of abdef
rare.out <- rare.out[Trait!='ABDEF',]
rare.out[Trait=='alloa',Trait:='combined_osteoarthritis']
rare.out[Trait=='hipkneeoa',Trait:='hip_knee_osteoarthritis']
rare.out[Trait=='hipoa',Trait:='hip_osteoarthritis']
rare.out[Trait=='kneeoa',Trait:='knee_osteoarthritis']
rare.out[grep('osteoarthritis$',Trait),Reference:='29559693']
rare.out[Trait=='bowes_psa',Trait:='psoriatic_arthritis_unpublished']
rare.out[Trait=='IgA_nephropathy',Trait:='iga_nephropathy']
rare.out[Trait=='iga_nephropathy',Reference:='25305756']
rare.out[Trait=='cousminer_lada',Trait:='latent_autoimmune_diabetes_in_adults']
rare.out[Trait=='latent_autoimmune_diabetes_in_adults',Reference:='30254083']
rare.out[Trait=='li_as',Trait:='ankylosing_spondylitis_li']
rare.out[Trait=='ankylosing_spondylitis_li',Reference:='30946743']
rare.out[Trait=='mahajan_t2d',Trait:='type_2_diabetes']
rare.out[Trait=='type_2_diabetes',Reference:='30297969']
rare.out[Trait=='na_psa',Trait:='psoriatic_arthritis_north_america']
rare.out[Trait=='psoriatic_arthritis_north_america',Reference:='30552173']
rare.out[Trait=='span_psa',Trait:='psoriatic_arthritis_spanish']
rare.out[Trait=='psoriatic_arthritis_spanish',Reference:='30552173']
rare.out[Trait=='hasnoot_uveitis_jia',Trait:='uveitis_jia']
rare.out[Trait=='uveitis_jia',Reference:='29513936']
rare.out[`First Author`=='ahola-olli',Reference:='27989323']
rare.out[`First Author`=='kuiper',Reference:='24957906']
rare.out[`First Author`=='taylor',Reference:='29795407']
rare.out[`First Author`=='ferreira',Reference:='29083406']
rare.out[`First Author`=='tian',Reference:='28928442']
rare.out[`First Author`=='astle',Reference:='27863252']

## how many snps are dropped from each trait when projecting

rare.out[Trait=='psoriatic_arthritis_spanish',`First Author`:='aterido']
rare.out[Trait=='psoriatic_arthritis_north_america',`First Author`:='aterido']


rare.out <- rare.out[order(N1+N0),.(Trait,Reference,N0,N1,`Missing SNPs`=missing),by=`First Author`]
write.xlsx(rare.out,file=file.path(OUT_DIR,'supplementary_tables_1_4.xlsx'),sheet='Other',row.names=FALSE,append=TRUE)
