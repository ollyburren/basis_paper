library(data.table)
library(magrittr)

## Process Westra data


DT <- fread('zcat /home/ob219/rds/rds-cew54-wallace-share/Data/expr/Westra-blood-2012/Westra-2012-12-21-CisAssociationsProbeLevelFDR0.5.txt.gz')

snps <- DT[,.(id=SNPName,chr=SNPChr,pos=SNPChrPos,alleles=SNPType,eallele=AlleleAssessed)]
snps[,pid:=paste(chr,pos,sep=':')]
## remove non biallelic SNPs by position
snps<-snps[snps[,.I[unique(eallele) %>% length==1],by=pid]$V1,] %>% unique
## there are still snps duplicated  - remove
snps <- snps[!pid %in% snps[duplicated(pid),],]
## this means that all effect alleles are the same across genes as we might expect.

## these need to be lifted over as on the wrong build
library(rtracklayer)

westra.36.gr <- with(snps,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=pos,width=1L),id=id))
c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
westra.37.gr<-unlist(liftOver(westra.36.gr,c))
DT.37 <- data.table(id=westra.37.gr$id,position.37=start(westra.37.gr))
snps <- merge(snps,DT.37,by.x='id',by.y='id',all.x=TRUE)
## 173 snps don't match after coord conversion
snps <- snps[!is.na(position.37),]
snps[,c('tmp.a1','tmp.a2'):=tstrsplit(alleles,'/')]
snps[eallele==tmp.a1,c('a1','a2'):=list(tmp.a2,tmp.a1)]
snps[eallele==tmp.a2,c('a1','a2'):=list(tmp.a1,tmp.a2)]
snps <- snps[,.(id,chr,position=position.37,a1,a2)]
snps[,pid:=paste(chr,position,sep=':')]
## next add in UK10K stuff so we can align
snp.DT <- fread('/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab')
m.DT <- merge(snp.DT,snps,by.x='pid',by.y='pid')
m.DT[ref_a1==a1 & ref_a2==a2,flip:=FALSE]
m.DT[ref_a1==a2 & ref_a2==a1,flip:=TRUE]

## get the list of variants where we need to flip the effect !

flip.snps <- m.DT[flip==TRUE,]$id
exp.DT <- DT[,.(id=SNPName,Z=OverallZScore,P=PValue,Gene=HUGO,probe=ProbeName)]
## add in extra rows for where genes are covered by the same probes
sgenes <- strsplit(exp.DT$Gene,',')
## identify cases where these are comma delimited but are the same hugo e.g. COG8
idx<-which(sapply(sgenes,function(x) duplicated(x) %>% any))
exp.DT[idx,Gene:=sapply(sgenes[idx],unique)]
sgenes <- strsplit(exp.DT$Gene,',')


rf <- rep(1:length(sgenes),sapply(sgenes,length))
exp.DT <- exp.DT[rf,]
exp.DT[,c('old.Gene','Gene'):=list(Gene,unlist(sgenes))]
f.DT <- merge(m.DT,exp.DT,by.x='id',by.y='id')
## split out genes that are annotated to more than one gene
## compute betas
f.DT[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
f.DT[,beta:=Z/sqrt(maf * (1-maf))]
##
f.DT[flip==TRUE,beta:=beta*-1]
f.DT[,uid:=paste(Gene,probe,sep=':')]

## create_ds_matrix from cupcake expects an odds ratio which it then takes log(OR) we therefore need an
## alternative method to create input for projection


SHRINKAGE_METHOD<-'ws_emp'
SHRINKAGE_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/shrinkage_june10k.RDS'
BASIS_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/basis_june10k.RDS'

shrink.DT <- readRDS(SHRINKAGE_FILE)
shrink.DT<-shrink.DT[,c('pid',shrink=sprintf("%s_shrinkage",SHRINKAGE_METHOD)),with=FALSE]
setkey(shrink.DT,'pid')
pc.emp <- readRDS(BASIS_FILE)
all.DT <- f.DT[,.(pid,uid,beta)]
## add shrinkage here
all.DT <- merge(all.DT,shrink.DT,by.x='pid',by.y='pid')[,shrunk.beta:=beta * ws_emp_shrinkage][,.(pid,uid,shrunk.beta)]


## next add in a dummy gene that has data for all variants
dummy.DT <- snp.DT[,.(pid,uid='DUMMY:-999',shrunk.beta=0)]
all.DT <- rbind(all.DT,dummy.DT)
all.DT <- melt(all.DT,id.vars=c('pid','uid'),measure.vars='shrunk.beta')
setkey(all.DT,'pid')

### easisest to do in batches of 100 ?
all.genes <- unique(all.DT[uid!='DUMMY:-999',]$uid)
chunk.gene <- split(all.genes, ceiling(seq_along(all.genes)/100))
res.DT <- mclapply(chunk.gene,function(gc){
  message("processing")
  tmp <- c(gc,'DUMMY:-999')
  r.DT <- dcast(all.DT[uid %in% tmp,],pid~uid+variable,fill=0)
  mat <- as.matrix(r.DT[,-1])
  rownames(mat) <- r.DT[[1]]
  bc <- predict(pc.emp,newdata=t(mat))
  data.table(trait = rownames(bc) %>% gsub("_shrunk.beta","",.),bc)
  ## these are now ready for projection but need checking
},mc.cores=8) %>% rbindlist
library(cowplot)

basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x,cat='basis')
res.DT[,cat:='westra']

ggplot(rbind(res.DT,basis.DT),aes(x=PC1,y=PC2,label=trait,col=cat)) + geom_point() + geom_text()

## remapping

library(httr)
library(jsonlite)
library(xml2)

server <- "https://rest.ensembl.org"
ext <- "/regulatory/species/homo_sapiens/microarray?"

r <- GET(paste(server, ext, sep = ""), content_type("application/json")) %>% content %>% toJSON %>% fromJSON %>% data.table


saveRDS(rbind(res.DT,basis.DT),file='/home/ob219/rds/rds-cew54-wallace-share/as_basis/westra/westra.RDS')
