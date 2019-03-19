library(rtracklayer)
## code to prepare for annotating snps discovered by cfdr technique
## for each basis trait get variants above a certain threshold - note this generates bcf files and only
## needs to be run once
PTHRESH <- 1

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


ms[,trait:='MS']
all <- list(uc,cd,ra,t1d,ms,pbc,cel,sle,psc,ast) %>% rbindlist

## plot the number of SNPs by cohort
by.cohort <- all[,list(snp.count=.N),by=trait]
by.snp <- all[,list(snp.count=.N),by=pid]
## remove duplicates
by.snp <- by.snp[!is.na(snp.count),]
by.snp <-by.snp[snp.count<11,]
by.snp <- by.snp[,list(bs=.N),by=snp.count]
by.snp[,snp.count:=factor(snp.count,levels=1:10)]
library(cowplot)

 by.cohort[,trait:=factor(trait,levels=by.cohort[order(snp.count,decreasing=TRUE),]$trait)]

pp1 <-ggplot(by.cohort,aes(x=trait,y=snp.count)) + geom_bar(stat="identity") + ylab("# SNPs") + xlab("Trait") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pp2 <- ggplot(by.snp,aes(x=snp.count,y=bs)) + geom_bar(stat="identity") + ylab("# Unique SNPs") + xlab("# Studies SNP shared across")

save_plot(plot_grid(pp1,pp2,labels=c('a','b')),file="~/tmp/snp_selection.pdf",base_width=8)

## what is responsible for peak at 6 ?

bs <- all[,list(snp.count=.N),by=pid]
gsu<-all[pid %in% bs[snp.count==1,]$pid,]
gs<-all[pid %in% bs[snp.count==6,]$pid,]

## we don't deal with the X chromosome
all <- all[!grepl("^X:",pid),]
## get a list of unique positions to make querying of bcf files much quicker
all[,c('chr','start'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
## remove x and y variants that are 23 and 24 respectively
all <- all[!chr %in% c(23,24),]
