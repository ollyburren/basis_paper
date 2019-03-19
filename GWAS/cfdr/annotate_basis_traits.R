
if(FALSE){
  library(rtracklayer)
  ## code to prepare for annotating snps discovered by cfdr technique
  ## for each basis trait get variants above a certain threshold - note this generates bcf files and only
  ## needs to be run once
  PTHRESH <- 1e-5

  ##1UC

  dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/uc_build37_45975_20161107.txt")
  uc <- dt[P.value<PTHRESH,.(trait='UC',pid=gsub("^([^\\_]+)\\_.*","\\1",MarkerName),a1=Allele1 %>% toupper,a2=Allele2 %>% toupper,p.value=P.value,beta=Effect,se.beta=StdErr)]
  ##2CD
  dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/cd_build37_40266_20161107.txt")
  cd <- dt[P.value<PTHRESH,.(trait='CD',pid=gsub("^([^\\_]+)\\_.*","\\1",MarkerName),a1=Allele1 %>% toupper,a2=Allele2 %>% toupper,p.value=P.value,beta=Effect,se.beta=StdErr)]
  ##3RA
  dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/RA_GWASmeta_European_v2.txt")
  setnames(dt,make.names(names(dt)))
  ra <- dt[P.val<PTHRESH,.(trait='RA',pid=paste(Chr,Position.hg19.,sep=':'),a1=A1,a2=A2,p.value=P.val,beta=log(OR.A1.),se.beta=abs(log(OR_95.CIlow))/1.96)]
  ##4T1D
  dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/t1d_cooper_2017.txt")
  t1d <- dt[p.meta<PTHRESH,.(trait='T1D',pid=paste(chromosome,position,sep=':'),a1=a0,a2=a1,p.value=p.meta,beta=beta.meta,se.beta=se.meta)]
  ##5MS
  dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/MSAllDataTableREDUCED.txt")
  ms <- dt[combined.pval<PTHRESH,.(trait='MS',chr,p36=position,a1=TestAllele,a2=OtherAllele,p.value=combined.pval,beta=combined.beta,se.beta=combined.se.beta,id=1:.N)]
  ms.36.gr <- with(ms,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=p36,width=1L),id=id))
  c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
  ms.37.gr<-unlist(liftOver(ms.36.gr,c))
  DT.37 <- data.table(id=ms.37.gr$id,p37=start(ms.37.gr))
  ms <- merge(ms,DT.37,by.x='id',by.y='id',all.x=TRUE)[,.(trait='ms',pid=paste(chr,p37,sep=':'),a1,a2,p.value,beta,se.beta)]
  ##6PBC
  dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/PBC_GCMETA_fixedeffects")
  pbc <- dt[P_value<PTHRESH,.(trait='PBC',chr,p36=pos,a1=allele_A,a2=allele_B,p.value=P_value,beta,se.beta=se,id=1:.N)]
  pbc.36.gr <- with(pbc,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=p36,width=1L),id=id))
  c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
  pbc.37.gr<-unlist(liftOver(pbc.36.gr,c))
  DT.37 <- data.table(id=pbc.37.gr$id,p37=start(pbc.37.gr))
  pbc <- merge(pbc,DT.37,by.x='id',by.y='id')[,.(trait='PBC',pid=paste(chr,p37,sep=':'),a1,a2,p.value,beta,se.beta)]
  ##7CEL
  dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/hg19_gwas_cel_dubois_4_19_1.tab")
  setnames(dt,make.names(names(dt)))
  cel <- dt[PValue<PTHRESH,]
  cel[,c('a1','a2'):=tstrsplit(Alleles.Maj.Min.,'>')]
  cel <- cel[,.(trait='CEL',pid=paste(Chr,Position,sep=':'),a1,a2,p.value=PValue,beta=log(OR.MinAllele.),se.beta=(abs(log(OR.MinAllele.))/(qnorm(PValue,lower.tail=FALSE)/2)))]
  ##8SLE
  dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/sle_benett_2016.tab")
  sle <- dt[PValueAdd<PTHRESH,.(trait='SLE',pid=paste(Chr,Distance,sep=':'),a1=AlleleA,a2=AlleleB,p.value=PValueAdd,beta=Beta1Add,se.beta=Se1Add)]
  ##9PSC
  ## note that to me the se.beta look off but this is probably due to number of decimal places we have
  dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/ipscsg2016.result.combined.full.with_header.txt")
  setnames(dt,make.names(names(dt)))
  psc <- dt[p<PTHRESH,.(trait='PSC',pid=paste(X.chr,pos,sep=':'),a1=allele_0,a2=allele_1,p.value=p,beta=log(or),se.beta=se)]
  ## this is code to see that se and p.vals make sense which they don't appear to
  ## probably due to the accuracy of or and se reported ?
  if(FALSE){
    psc[,check.se.beta:=abs(beta)/(qnorm(p.value/2,lower.tail=FALSE))]
    psc[,check.p.value:=pnorm(abs(beta/se.beta),lower.tail=FALSE) * 2]
    par(mfrow=c(1,2))
    with(psc,plot(se.beta,check.se.beta))
    abline(a=0,b=1,col='red')
    with(psc,plot(-log10(p.value),-log10(check.p.value)))
    abline(a=0,b=1,col='red')
  }
  ##ASTHMA
  dt <- fread("zcat /home/ob219/share/Data/GWAS-summary/asthma/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv.gz")
  #dt <- fread("/home/ob219/rds/hpc-work/as_basis/gwas_stats/raw/gabriel.csv")
  #ast <- dt[P_fix<PTHRESH,.(Chr,p36=position,a1=Allele_1,at2=Allele_2,p.value=P_fix,beta=log(OR_fix),se.beta=abs(log(ORl_fix))/1.96,id=1:.N)]
  ast <- dt[European_ancestry_pval_fix<PTHRESH,.(trait='AST',pid=paste(chr,position,sep=':'),a1=reference_allele,a2=alternate_allele,p.value=European_ancestry_pval_fix,beta=European_ancestry_beta_fix,se.beta=European_ancestry_se_fix)]
  #ast.36.gr <- with(ast,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=p36,width=1L),id=id))
  #c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
  #ast.37.gr<-unlist(liftOver(ast.36.gr,c))
  #DT.37 <- data.table(id=ast.37.gr$id,p37=start(ast.37.gr))
  #ast <- merge(ast,DT.37,by.x='id',by.y='id',all.x=TRUE)[,.(trait='AST',pid=paste(chr,p37,sep=':'),a1,a2,p.value,beta,se.beta)]

  all <- list(uc,cd,ra,t1d,ms,pbc,cel,sle,psc,ast) %>% rbindlist
  ## we don't deal with the X chromosome
  all <- all[!grepl("^X:",pid),]
  ## get a list of unique positions to make querying of bcf files much quicker
  all[,c('chr','start'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
  ## remove x and y variants that are 23 and 24 respectively
  all <- all[!chr %in% c(23,24),]
  saveRDS(all,file="/home/ob219/share/as_basis/GWAS/support/basis_trait_suggestive_snps.RDS")


  all.pos<-all[!duplicated(pid),.(chr,start,pid)]

  #all.pos.gr <- with(unique(all.pos),GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,width=1L)),pid=pid)
  ## add in basis snps
  basis <- fread("/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab")
  basis[,c('chr','start'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
  full <- rbind(all.pos,basis[,.(chr,start,pid)])[!duplicated(pid),]
  full[,chr:=paste('chr',chr,sep='')]
  OUTDIR <- '/home/ob219/share/as_basis/GWAS/snp_manifest/lookup'
  options('scipen'=999)
  lapply(split(full,full$chr),function(dt){
    chrname<-dt[1,]$chr
    fname <- file.path(OUTDIR,sprintf("%s.txt",chrname))
    dt <- dt[order(start),]
    write.table(dt[,.(chr,start)],file=fname,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  })
  options('scipen'=0)
  ## use bcftools to dump bcf files with only the variants we care about in on a chr basis
  UK10KDIR <- '/home/ob219/share/Data/reference/UK10K/BCF/'
  BCFOUTDIR <- '/home/ob219/share/as_basis/GWAS/snp_manifest/lookup/uk10k'
  for(chrname in paste('chr',1:22,sep='')){
    message(chrname)
    ukf <- file.path(UK10KDIR,sprintf("%s.bcf.gz",chrname))
    rfile <- file.path(OUTDIR,sprintf("%s.txt",chrname))
    ofile <- file.path(BCFOUTDIR,sprintf("%s.bcf",chrname))
    cmd <- sprintf("bcftools view -R %s -Ob %s > %s",rfile,ukf,ofile)
    write(cmd,file="~/tmp/qstuff/filter.txt",append=TRUE)
    #system(cmd)
    ## run these outside of R so as to preserve session from watchdog
  }
}

## this code allows one to lookup a basis SNP and get all variants with 'interval'
## with an r2 above min.r2 the curated set includes variants from each basis trait
## above a p.value threshold (PTHRESH)
library(snpStats)
library(rtracklayer)
BCFOUTDIR <- '/home/ob219/share/as_basis/GWAS/snp_manifest/lookup/uk10k'
bcf2snpmatrix <- function(pid,quiet=FALSE,interval=1e6,min.r2=0.1){
  message(pid)
  ## convert pid into a region based on interval
  chr <- sub("^([^:]+):.*","\\1",pid); start <- sub("^[^:]+:(.*)","\\1",pid) %>% as.numeric
  region <- sprintf('chr%s:%d-%d',chr,pmax(start-interval,1),start+interval)
  bcf.file <- file.path(BCFOUTDIR,sprintf("chr%s.bcf",chr))
  header_cmd <- sprintf("bcftools view -h %s",bcf.file)
  if(!quiet)
    message(header_cmd)
  my.pipe<-pipe(header_cmd)
  header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
  close(my.pipe)
  cnames<-unlist(strsplit(header,"\t"))
  bcftools_cmd<-sprintf("bcftools view -r %s -O v %s | grep -v '^#'",region,bcf.file)
  if(!quiet)
    message(bcftools_cmd)
  #tryCatch(as.data.frame(fread(cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE)),error=function(e){print(sprintf("Error=%s fread CMD=%s",e,cmd));return(NA)})
  tmp<-tryCatch(fread(bcftools_cmd),error=function(e){print(sprintf("Error=%s fread CMD=%s",e,bcftools_cmd));return(data.table(pid=pid,error=TRUE))})
  #print(tmp)
  #message(grepl("error",names(tmp)))
  if(any(grepl("error",names(tmp))))
    return(tmp)
  setnames(tmp,cnames)
  gt<-tmp[,10:ncol(tmp),with=FALSE]
  if(nrow(gt)==0)
    return(NA)
  info<-tmp[,1:9,with=FALSE]
  setnames(info,'#CHROM','CHROM')
  if(!quiet)
    message("Creating snpMatrix obj")
  sm<-apply(gt,1,function(x) sub("0\\|0","1",x))
  sm<-apply(sm,1,function(x) sub("(0\\|1)|(1\\|0)","2",x))
  sm<-apply(sm,1,function(x) sub("1\\|1","3",x))
  ## set anything else to a missing value
  sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\|[0-9]","0",x))))
  info[,pid:=paste(CHROM,POS,sep=':') %>% gsub("^chr","",.)]
  colnames(sm)<-info$pid
  rownames(sm)<-colnames(gt)
  ## compute r2 for selected pid with everything else
  sm <- new("SnpMatrix", sm)
  idx <- which(info$pid==pid)
  r2 <- ld(sm[,idx],sm,stats="R.squared") %>% t
  dt <- data.table(pid=rownames(r2),r2=r2[,1])[r2>min.r2,]
  dt
}

#bcf2snpmatrix(test.pid)

#pid <- '4:102737250'
#bcf2snpmatrix(pid)

bcf2snpmatrix('10:9062856')

## this file contains annotated variants across basis traits that have p<PVALTHRESH
all <- readRDS('/home/ob219/share/as_basis/GWAS/support/basis_trait_suggestive_snps.RDS')
## add ld block designations
ld <- fread('/home/ob219/share/as_basis/GWAS/support/all.1cM.tab')
ld.gr <- with(ld,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end)),id=1:nrow(ld))
all.gr <- with(all,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,width=1L)))
ol <- findOverlaps(all.gr,ld.gr) %>% as.matrix
all[ol[,1],ld.block:=ol[,2]]


### this is chris' code

## basis p values
## install_github('ollyburren/cupcake')
library(cupcake)
TRAIT_MANIFEST <- '~/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
GWAS_DATA_DIR <- '~/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'~/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)

RESULTS.FILE <- '~/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)
res.DT[,trait:=as.character(trait)]
res.DT[,trait:=sub("bb_","",trait)]
res.DT[,trait:=make.names(trait)]

## annotate by genes
library(GenomicRanges)
if(!file.exists("~/genes.RData")) {
    library(randomFunctions)
    library(biomaRt)
    mart <- useMart('ENSEMBL_MART_ENSEMBL',host="grch37.ensembl.org")
    ensembl.gene <- useDataset("hsapiens_gene_ensembl",mart=mart)

    genes<-getBM(
        ## filters= c("chromosome_name","start","end"),
        attributes= c('ensembl_gene_id','chromosome_name','start_position','end_position','external_gene_name','strand','gene_biotype'),
        ## values= ss,
        mart= ensembl.gene)
   genes <- subset(genes, !grepl("HG|HS|GL",chromosome_name))
    save(genes,file="~/genes.RData")
} else {
    load("~/genes.RData")
    genes <- subset(genes, !grepl("HG|HS|GL",chromosome_name))
    genes <- subset(genes, gene_biotype=="protein_coding")
    ggr <- GRanges(Rle(genes$chromosome_name),
                   IRanges(start=genes$start_position,
                           end=genes$end_position))
    mcols(ggr) <- genes[,c("ensembl_gene_id","external_gene_name")]
}

library(GenomicRanges)
annot <- function(pid) {
    ss <- strsplit(pid,":")
    chr <- sapply(ss,"[[",1)
    bp <- sapply(ss,"[[",2)  %>% as.numeric()
    left <- bp-1e+5
    right <- bp + 1e+5
    gr0 <- GRanges(seqnames=Rle(chr,rep(1,length(chr))),
                  IRanges(start=left,end=right),pid=pid)
    mol <- mergeByOverlaps(gr0,ggr)
    mol <- split(mol$external_gene_name,mol$pid) %>% sapply(.,function(pg){unique(pg) %>% paste(.,collapse=',')})
    data.table(pid=names(mol),genes=mol)
}

annotate <- function(pid){
  ## look up basis traits and r2
  tmp <- lapply(unique(pid),function(ipid){
    r2m <- bcf2snpmatrix(ipid)
    r2m[,qpid:=ipid]
    merge(r2m,all,by='pid')
  }) %>% rbindlist(.,fill=TRUE)
  if(any(grepl("error",names(tmp)))){
    tmp<-tmp[is.na(error),]
    tmp$error<-NULL
  }
  #tmp <- tmp[tmp[, .I[p.value == min(p.value)], by=c('ld.block','trait')]$V1]

  tmp <- tmp[,.(basis.trait=trait,basis.trait.pid=pid,qpid,r2=signif(r2,digits=2),p.value,ld.block)]
  tmp <- merge(tmp,annot(pid),by.x='qpid',by.y='pid',all.x=TRUE)
  setnames(tmp,'qpid','pid')
  tmp
}



################################################################################

## myositis

library(data.table)
library(magrittr)
library(cowplot)

res <- readRDS("~/share/as_basis/GWAS/myogen_myositis/myogen_cut_down_projection_results_18_02_2019.RDS")
ggplot(res,aes(x=delta,xmin=delta-1.96*sqrt(variance),col=p.adj<0.01,xmax=delta+1.96*sqrt(variance),y=variable)) + geom_errorbarh() + geom_point() + geom_vline(xintercept=0) + facet_wrap(~trait)

a <- 5e-8/11 # bonferroni of gwsig 11 components
a <- 0.05/11 # bonferroni of 0.05 11 components
## todo <- res[grepl("myogen",trait) & p.adj < 0.01,]
todo <- res[grepl("myogen",trait) & p.value < a ,]
todo

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_myositis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
## get the rotations
rot.dt <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
## rot.dt <- melt(rot.dt,id.vars='pid')
m <- melt(rot.dt,"pid")

files <- list.files("/home/ob219/share/as_basis/GWAS/myogen_myositis/",
                    pattern="cut_down_source.RDS",full=TRUE)
files
x <- lapply(files,readRDS)  %>% lapply(., function(tmp) tmp[,.(trait,pid,P,shrink=ws_emp_shrinkage)]) %>% rbindlist()
x<- merge(x,m,by='pid',allow.cartesian = TRUE)
## filter
x <- merge(x,todo[,.(trait,variable,overall.fdr=p.adj)],by=c("trait","variable"))

## calc fdr
x[,fdr:=as.numeric(NA)]
x[,unq:=abs(value)]
x[,use:=unq>quantile(unq,0.998) ,by=variable]
x[,fdr0:=p.adjust(P,method="BH"),by=variable]
x[use==TRUE & variable=="PC1",fdr:=p.adjust(P,method="BH"),by=variable]
myo.results <- x[fdr<0.05,.(pid,trait,raw.p=P,pc=variable,basis.fdr=fdr,fdr=fdr0)]

myo.anot <- merge(myo.results,annotate(myo.results$pid),by='pid',all.x=TRUE)



##x[fdr0<0.05,]
##              trait variable         pid         P     shrink       value
## 1: myositis_myogen      PC1 2:198929806 4.561e-06 0.03286767  0.01709288
## 2:       pm_myogen      PC1 1:114377568 2.585e-05 0.10905678 -0.26730173
##     overall.fdr          unq  use      fdr0        fdr
## 1: 1.057988e-12 0.0005618031 TRUE 0.5424459 0.01008437
## 2: 1.435375e-06 0.0291510653 TRUE 0.6562213 0.02857718

##isnps <- x[fdr<0.05,]$pid  %>% unique()
#x[pid %in% isnps,]
#x[pid=="1:114377568",]
#basis.DT[pid %in% isnps,as.list(summary(p.value)),by=pid]
#basis.DT[pid %in% isnps & p.value < 1e-6,]

##annot(isnps) # PTPN22, PLCL1

###############################################################################

## egpa / vasc

library(data.table)
library(magrittr)

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
## get the rotations
rot.dt <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
## rot.dt <- melt(rot.dt,id.vars='pid')
m <- melt(rot.dt,"pid")

## read in vasc snp data
files <- list.files("/home/ob219/share/as_basis/GWAS/vasc/projections/",
                    pattern="_source.RDS",full=TRUE)
x <- lapply(files,readRDS)  %>% rbindlist()
x<- merge(x[,.(trait,pid,P=p.value,shrink=ws_emp_shrinkage)],m,by='pid',allow.cartesian = TRUE)

## filter by significance
todo <- res.DT[grepl("egpa|vasc|mpo|anca",trait,ignore.case=TRUE) & p.adj<0.05,]
todo
x <- merge(x,todo[,.(variable,trait,fdr.overall=p.adj)],by=c("variable","trait"))

x[,fdr:=as.numeric(NA)]
x[,unq:=abs(value)]
##x[,unq:=abs(value)]
x[,use:=unq>quantile(unq,0.998),by=c("trait","variable")]
x[,fdr0:=p.adjust(P,method="BH"),by=c("trait","variable")]
x[use==TRUE,fdr:=p.adjust(P,method="BH"),by=c("trait","variable")]

vasc.results <- x[fdr<0.05,.(pid,trait,raw.p=P,pc=variable,basis.fdr=fdr,fdr=fdr0)]
vasc.anot <- merge(vasc.results,annotate(vasc.results$pid),by='pid',all.x=TRUE,allow.cartesian=TRUE)

summary(x[,.(fdr,fdr0)])

x[fdr<0.01,]
x[fdr0<0.01,]
unique(x[fdr<0.01,.(n=.N),by=c("trait","pid")])
x[fdr<0.01,.(n=.N,traits=paste(sort(unique(trait)),collapse="/")),by=c("pid")]
##             pid n        traits
##  1: 1:114377568 3           mpo
##  2: 7:128594183 3           mpo
##  3: 7:128617466 3           mpo
##  4:  10:9062856 5 anca_Neg/egpa
##  5: 15:67464013 5 anca_Neg/egpa
##  6: 1:117263868 3 anca_Neg/egpa
##  7: 5:131784393 2      anca_Neg
##  8: 5:131796922 5 anca_Neg/egpa
##  9: 5:131819921 5 anca_Neg/egpa
## 10:  6:90814199 5 anca_Neg/egpa
## 11:  6:90996769 5 anca_Neg/egpa
## 12:  7:92236164 2 anca_Neg/egpa
## 13: 5:110567598 2          egpa
## 14: 13:42989660 2          egpa
## 15: 2:111913056 3 anca_Neg/egpa

isnps <- x[fdr<0.01,]$pid  %>% unique()
basis.DT[pid %in% isnps,as.list(summary(p.value)),by=pid]
with(basis.DT[pid %in% isnps & p.value < 1e-6,],table(pid,trait))

tmp <- as.data.frame(annot(isnps))[,c("seqnames","start","end","ensembl_gene_id","external_gene_name")]
tmp[order(tmp$seqnames,tmp$start),]
##    seqnames     start       end ensembl_gene_id external_gene_name
## 7        15  67356101  67487533 ENSG00000166949              SMAD3
## 8        15  67435072  67439277 ENSG00000259202      RP11-342M21.2
## 9        15  67493371  67547533 ENSG00000103591              AAGAB
## 10       15  67547138  67794598 ENSG00000103599               IQCH
## 11       13  42846289  42897396 ENSG00000023516             AKAP11
## 15        7  92116334  92157845 ENSG00000127980               PEX1
## 17        7  92158087  92167319 ENSG00000127993              RBM48
## 18        7  92190107  92219708 ENSG00000234545            FAM133B
## 21        7  92234235  92465908 ENSG00000105810               CDK6
## 16        7 128470431 128499328 ENSG00000128591               FLNC
## 19        7 128502880 128505898 ENSG00000128524            ATP6V1F
## 20        7 128577666 128590089 ENSG00000128604               IRF5
## 22        7 128594948 128695198 ENSG00000064419              TNPO3
## 12        2 111490150 111875799 ENSG00000153093              ACOXL
## 13        2 111876955 111926024 ENSG00000153094            BCL2L11
## 14        6  90636248  91006627 ENSG00000112182              BACH2
## 4         5 110559351 110830584 ENSG00000152495              CAMK4
## 1         5 131705444 131731306 ENSG00000197375            SLC22A5
## 2         5 131746328 131811736 ENSG00000197536            C5orf56
## 3         5 131817301 131826490 ENSG00000125347               IRF1
## 5         5 131877136 131892530 ENSG00000113525                IL5
## 6         5 131891711 131980313 ENSG00000113522              RAD50
## 23        1 114239453 114302111 ENSG00000116793              PHTF1
## 24        1 114304454 114355098 ENSG00000081019              RSBN1
## 25        1 114356433 114414381 ENSG00000134242             PTPN22
## 26        1 114420790 114430169 ENSG00000188761            BCL2L15
## 27        1 114437370 114447823 ENSG00000134262              AP4B1
## 28        1 114447763 114456708 ENSG00000118655            DCLRE1B
## 29        1 114471814 114520426 ENSG00000163349              HIPK1
## 30        1 117117031 117210375 ENSG00000143061              IGSF3
## 31        1 117236734 117249225 ENSG00000203864           C1orf137
## 32        1 117297007 117311850 ENSG00000116824                CD2
################################################################################
###############################################################################

## JIA

library(data.table)
library(magrittr)

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
## get the rotations
rot.dt <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
## rot.dt <- melt(rot.dt,id.vars='pid')
m <- melt(rot.dt,"pid")

## read in vasc snp data
files <- list.files("/home/cew54/share/as_basis/GWAS/sum_stats",
                    pattern="jia",full=TRUE)
x <- lapply(files,function(f) {
    tmp <- fread(f)
    tmp$trait <- sub("_unpub.tab","",basename(f))
    tmp
})%>% rbindlist()
head(x)
x<- merge(x[,.(trait,pid,P=p.value)],m,by='pid',allow.cartesian = TRUE)

## filter by significance
todo <- res.DT[grepl("jia",trait,ignore.case=TRUE) & p.adj<0.05,]
todo
todo$trait <- tolower(sub("_","",todo$trait))
table(todo$trait)
unique(x$trait)
all(unique(todo$trait) %in% unique(x$trait))
x <- merge(x,todo[,.(variable,trait,fdr.overall=p.adj)],by=c("variable","trait"))

x[,fdr:=as.numeric(NA)]
x[,unq:=abs(value)]
##x[,unq:=abs(value)]
x[,use:=unq>quantile(unq,0.998),by=c("trait","variable")]
x[,fdr0:=p.adjust(P,method="BH"),by=c("trait","variable")]
x[use==TRUE,fdr:=p.adjust(P,method="BH"),by=c("trait","variable")]
jia.results <- x[fdr<0.05,.(pid,trait,raw.p=P,pc=variable,basis.fdr=fdr,fdr=fdr0)]
jia.anot <- merge(jia.results,annotate(jia.results$pid),by='pid',all.x=TRUE,allow.cartesian=TRUE)
all.annot.results <- rbindlist(list(myo.anot,vasc.anot,jia.anot))
saveRDS(all.annot.results,file="~/share/as_basis/GWAS/cfdr_results.RDS")
aF <- all.annot.results[r2>0.6,]

## which variants do we miss

all.annot.results[!pid %in% aF$pid,]

all.annot.results[,.SD[1,],by=c('ld.block','trait')]

af[pid %in% ]

aF[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
aF <- aF[order(trait,chr,pos),]
## really want only the max p.values for each disease and trait
## write out as an xls file so we can manually annotate and curate
aF <- aF[aF[,.I[p.value==min(p.value)],by=c('ld.block','trait','basis.trait')]$V1,]
aF[,]

aF <- aF[order(trait,gsub("PC","",pc) %>% as.numeric,chr,pos)]

library(xlsx)
write.xlsx(aF, file="/home/ob219/tmp/fdr_results.xlsx", sheetName="Sheet1",col.names=TRUE, row.names=FALSE, append=FALSE)
## I took this into excel and then began to do manual curation to identify in JIA which hits from cfdr might be worth
## reporting

jia <- read.xlsx(file="/home/ob219/tmp/fdr_results.xlsx", sheetName="jia") %>% data.table
jia <- jia[gw.sig.basis.trait=='Y',.(pid,trait,raw.p,pc,basis.fdr,basis.trait,basis.trait.pid,r2,basis.trait.p.value=p.value,ld.block,genes)]
save(jia,file="~/share/as_basis/GWAS/jia_fdr_results.RDS")

jia[,.(pid,trait)] %>% unique %>% table


results <- aF[aF[,.I[p.value==min(p.value)],by=c('ld.block','trait','basis.trait')]$V1,]

## what is the min p.value if we filter by r2>0.6

all.annot.results[r2>0.6,]


summary(x[,.(fdr,fdr0)])

x[fdr<0.05,]
x[fdr0<0.05,]
unique(x[fdr<0.05,.(n=.N),by=c("trait","pid")])
x[fdr<0.05,.(n=.N,traits=paste(sort(unique(trait)),collapse="/")),by=c("pid")]
isnps <- x[fdr<0.05,]$pid  %>% unique()
basis.DT[pid %in% isnps,as.list(summary(p.value)),by=pid]
with(basis.DT[pid %in% isnps & p.value < 1e-6,],table(pid,trait))

tmp <- as.data.frame(annot(isnps))[,c("seqnames","start","end","ensembl_gene_id","external_gene_name")]
tmp[order(tmp$seqnames,tmp$start),]
