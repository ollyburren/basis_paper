## this create input files for aligning alleles
## each of the input files is different so just a list of transformations.

library(data.table)
raw.dir<-'/scratch/ob219/as_basis/gwas_stats/raw/'
process.dir<-'/scratch/ob219/as_basis/gwas_stats/processed/'
col.names<-c('pid','a1','a2','or','p.value')
## by convention chrX - 23 and chrY - 24

##RA - OR with respect to A1

ra<-fread(file.path(raw.dir,'RA_GWASmeta_European_v2.txt'))
setnames(ra,make.names(names(ra)))
ra[,pid:=paste(Chr,Position.hg19.,sep=':')]
ra<-ra[,.(pid,A1,A2,OR.A1.,P.val)]
setnames(ra,col.names)
write.table(ra,file=file.path(process.dir,'ra_okada.tab'),quote=FALSE,sep="\t",row.names=FALSE)

##SLE
## this comes as an archive separated by chromosome this little script merges all together
if(FALSE){
  sle.files<-list.files(path='/Users/oliver/DATA/RAW_GWAS/sle_bennett_2016',pattern='*.txt',full.names=TRUE)
  all.sle.chr <- lapply(sle.files,fread)
  chrs<-gsub("^Chr([0-9XY]+).txt","\\1",basename(sle.files))
  chrs[chrs=="X"]<-23
  chrs[chrs=="XY"]<-24
  chrs<-as.numeric(chrs)
  bennett<-rbindlist(lapply(seq_along(chrs),function(i){
    chrom<-chrs[i]
    tmp.DT<-copy(all.sle.chr[[i]])
    tmp.DT[,Chr:=chrom]
    tmp.DT
  }))
  options(scipen=999)
  write.table(bennett,file="/Users/oliver/DATA/RAW_GWAS/sle_bennett_2016/all_chr.tab",sep="\t",quote=FALSE,row.names=FALSE)
  options(scipen=0)
}
sle <- fread(file.path(raw.dir,'sle_bentham_2016.tab'))
options(scipen=999)
sle[,pid:=paste(Chr,Distance,sep=':')]
options(scipen=0)
sle <- sle[,.(pid,AlleleA,AlleleB,OddR,PValueAdd)]
setnames(sle,col.names)
write.table(sle,file=file.path(process.dir,'sle_bennett.tab'),quote=FALSE,sep="\t",row.names=FALSE)

##CELIAC with respect to A2 (minor allele in original dataset)
## No longer have the raw data will have to use processed ImmunoBase data
cel<-fread(paste0(raw.dir,'hg19_gwas_cel_dubois_4_19_1.tab'))
setnames(cel,make.names(names(cel)))
cel[,c('a1','a2'):=tstrsplit(Alleles.Maj.Min.,'>')]
cel[Chr=='X',Chr:=23]
cel[Chr=='Y',Chr:=24]
options(scipen=999)
cel[,pid:=paste(Chr,Position,sep=':')]
options(scipen=0)
cel <- cel[,.(pid,a1,a2,OR.MinAllele.,PValue)]
setnames(cel,col.names)
write.table(cel,file=file.path(process.dir,'cel_dubois.tab'),quote=FALSE,sep="\t",row.names=FALSE)

##MS with respect to A1 (The test allele) but needs double checking in immunobase

ms <- fread(file.path(raw.dir,"MSAllDataTableREDUCED.txt"))
## these need transferring across to build 37 as currently 36 !
ms.pos.lu <- fread(file.path(raw.dir,"hg19_gwas_ms_imsgc_4_19_2.tab"))[,.(Marker,Chr,Position)]
options(scipen=999)
ms.pos.lu[,pid:=paste(Chr,Position,sep=':')]
options(scipen=0)
setkey(ms.pos.lu,Marker)
setkey(ms,rsid)
ms<-ms[ms.pos.lu[,.(Marker,pid)]]
ms<-ms[,.(pid,TestAllele,OtherAllele,or=exp(combined.beta),combined.pval)]
setnames(ms,col.names)
write.table(ms,file=file.path(process.dir,'ms_imsgc.tab'),quote=FALSE,sep="\t",row.names=FALSE)


processIBD<-function(DT){
	DT$or<-exp(DT$Effect)
	DT<-DT[,c('MarkerName','Allele1','Allele2','P.value','or'),with=FALSE]
	DT$chr<-sub("([^:]+):.*","\\1",DT$MarkerName)
	DT$position<-sub("[^:]+:([0-9]+).*","\\1",DT$MarkerName)
	DT$position<-as.numeric(sub("[^:]+:([0-9]+).*","\\1",DT$MarkerName))
	##convert alleles to upper case
	for(n in c('Allele1','Allele2'))
	        DT[[n]]<-toupper(DT[[n]])
	## note that we switch a1 and a2 as alleles are wrt to a2
	setnames(DT,c('id','a1','a2','p.val','or','chr','position'))
  options(scipen=999)
  DT[,pid:=paste(chr,position,sep=':')]
  options(scipen=0)
  DT<-DT[,.(pid,a1,a2,or,p.val)]
  setnames(DT,col.names)
  DT
}

## CRO disease with respect to allele 2
cd<-fread(paste0(raw.dir,'cd_build37_40266_20161107.txt'))
cdp<-processIBD(copy(cd))
write.table(cdp,file=file.path(process.dir,'cd_delaange.tab'),quote=FALSE,sep="\t",row.names=FALSE)

## UC disease with respect to allele 2
uc<-fread(paste0(raw.dir,'uc_build37_45975_20161107.txt'))
ucp<-processIBD(copy(uc))
write.table(ucp,file=file.path(process.dir,'uc_delaange.tab'),quote=FALSE,sep="\t",row.names=FALSE)

# PSC - all wrt to risk allele which we is allele2
psc<-fread(paste0(raw.dir,'ipscsg2016.result.combined.full.with_header.txt'))
setnames(psc,make.names(names(psc)))
psc[X.chr=='X',X.chr:=23]
psc[X.chr=='Y',X.chr:=24]
options(scipen=999)
psc[,pid:=paste(X.chr,pos,sep=':')]
options(scipen=0)
psc <- psc[,.(pid,allele_0,allele_1,or,p)]
setnames(psc,col.names)
write.table(psc,file=file.path(process.dir,'psc_ji.tab'),quote=FALSE,sep="\t",row.names=FALSE)

## T1D - all wrt to allele 2
t1d <- fread(paste0(raw.dir,'t1d_cooper_2017.txt'))[!is.na(chromosome),]
options(scipen=999)
t1d[,pid:=paste(chromosome,position,sep=':')]
options(scipen=0)
t1d<-t1d[!is.na(p.meta),.(pid,a0,a1,OR.meta,p.meta)]
setnames(t1d,col.names)
write.table(t1d,file=file.path(process.dir,'t1d_cooper.tab'),quote=FALSE,sep="\t",row.names=FALSE)

## PBC
## Cannot find originals
pbc<-fread(paste0(raw.dir,'hg19_gwas_pbc_cordell_4_20_0.tab'))
setnames(pbc,make.names(names(pbc)))
pbc[,c('a1','a2'):=tstrsplit(Alleles.Maj.Min.,'>')]
pbc[Chr=='X',Chr:=23]
pbc[Chr=='Y',Chr:=24]
options(scipen=999)
pbc[,pid:=paste(Chr,Position,sep=':')]
options(scipen=0)
pbc <- pbc[,.(pid,a1,a2,OR.MinAllele.,PValue)]
setnames(pbc,col.names)
write.table(pbc,file=file.path(process.dir,'pbc_cordell.tab'),quote=FALSE,sep="\t",row.names=FALSE)

## JIA
## I computed these

jia.data.dir<-'/home/ob219/scratch/jia/by.trait/'
jia.fs<-list.files(path=jia.data.dir,pattern='*.RDS',full.names=TRUE)
all.jia<-lapply(jia.fs,function(f){
  message(sprintf("Processing %s",basename(f)))
	tmp<-readRDS(f)[,.(chr,position,a1,a2,beta,p.val)]
  options(scipen=999)
  tmp[,pid:=paste(chr,position,sep=':')]
  options(scipen=0)
  tmp<-tmp[,.(pid,a1,a2,or=exp(beta),p.val)]
  setnames(tmp,col.names)
  st<-tolower(gsub("\\.RDS","",basename(f)))
  ofile<-file.path(process.dir,sprintf("jia%s_unpub.tab",st))
  write.table(tmp,file=ofile,quote=FALSE,sep="\t",row.names=FALSE)
})

## process T2D results - late addition on 17/06/2019

library(annotSnpStats)

t2d.DT <- fread("zcat /home/ob219/share/Data/GWAS-summary/T2D-DIAGRAM/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt.gz")
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
man.DT <- fread(SNP_MANIFEST)
M <- merge(t2d.DT[,.(pid=paste(Chr,Pos,sep=':'),a1=EA,a2=NEA,or=exp(Beta),p.value=Pvalue,neff=Neff)],man.DT,by='pid')
alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
#alleles <- alleles[!duplicated(pid),]
#alleles <- M[,list(al.x=paste(uk10_A1,uk10_A2,sep='/'),al.y=paste(a1,a2,sep='/')),by='pid']
## to make quick
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$pid
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
alleles[,g.class:=align.class]
idx<-which(alleles$g.class=='impossible')
if(length(idx) >0){
  M <- M[-idx,]
  alleles <- alleles[-idx,]
}

## check direction which is the effect allele ? It appears that a1 is the effect allele
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
M <- M[!duplicated(pid),]
## so here alleles match we need to flip as we want wrt to a2
M <- M[g.class=='match',or:=1/or]
## note there are a few p.values in here that are less than 10e-300 that R allows
## thus they end up as 0 in the output file - these need to be changed to 1e-300 manually


#M <- fread("/home/ob219/share/as_basis/GWAS/sum_stats/t2d_mahajan.tab")
over.pid<-M[as.numeric(p.value)<=1e-300,]$pid
M[,p.value:=as.numeric(p.value)]
M[pid=='10:114808902',p.value:=3.3e-299]
M[pid=='10:114758349',p.value:=3.3e-300]
M[,p.value:=as.numeric(p.value)]
write.table(M[,.(pid,a1=ref_a1,a2=ref_a2,or,p.value)],file='/home/ob219/share/as_basis/GWAS/sum_stats/t2d_mahajan.tab',sep="\t",quote=FALSE,row.names=FALSE)

## process vitiligo data ## assume effect allele is a1 ?

files <- list.files(path='/home/ob219/share/Data/GWAS-summary/vitiligo-jin/raw',pattern='*.gz',full.names=TRUE)

vit.DT <- lapply(files,function(f){
  message(f)
  tmp.DT<-sprintf("zcat %s",f) %>% fread
  tmp.DT[,.(pid=paste(CHR,BP,sep=':'),A1,MAF,A2,CHISQ,P,ORX,SE,L95,U95)]
}) %>% rbindlist

SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
man.DT <- fread(SNP_MANIFEST)
M <- merge(vit.DT[,.(pid,a1=A1,a2=A2,or=ORX,p.value=P)],man.DT,by='pid')
alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
#alleles <- alleles[!duplicated(pid),]
#alleles <- M[,list(al.x=paste(uk10_A1,uk10_A2,sep='/'),al.y=paste(a1,a2,sep='/')),by='pid']
## to make quick
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$pid
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
alleles[,g.class:=align.class]
idx<-which(alleles$g.class=='impossible')
if(length(idx) >0){
  M <- M[-idx,]
  alleles <- alleles[-idx,]
}

## check direction which is the effect allele ? It appears that a1 is the effect allele
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
M <- M[!duplicated(pid),]
## so here alleles match we need to flip as we want wrt to a2
M <- M[g.class=='match',or:=1/or]
write.table(M[,.(pid,a1=ref_a1,a2=ref_a2,or,p.value)],file='/home/ob219/share/as_basis/GWAS/sum_stats/vitiligo_jin.tab',sep="\t",quote=FALSE,row.names=FALSE)
## new manifest
write.table(man.DT[pid %in% M$pid,],file="/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab",sep="\t",quote=FALSE,row.names=FALSE)
