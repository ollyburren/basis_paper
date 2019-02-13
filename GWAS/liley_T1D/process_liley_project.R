library(data.table)
library(magrittr)

## Process Westra data


DT <- fread('zcat /home/ob219/share/Data/GWAS-summary/t1d-aab-age-liley+alleles.txt.gz')

#DT.jamie <- fread('zcat /home/ob219/rds/hpc-work/jamie/summary_stats.txt.gz')

snps <- DT[,.(id=snp,chr=chr,pos=pos36,alleles=consensus.alleles)]
snps[,pid:=paste(chr,pos,sep=':')]
## there are still snps duplicated  - remove
snps <- snps[!pid %in% snps[duplicated(pid),],]
## this means that all effect alleles are the same across genes as we might expect.
## remove snps where we don't have an allele

#snps <- snps[alleles!='',]

## these need to be lifted over as on the wrong build
library(rtracklayer)

liley.36.gr <- with(snps,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=pos,width=1L),id=id))
c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
liley.37.gr <-unlist(liftOver(liley.36.gr ,c))
DT.37 <- data.table(id=liley.37.gr$id,position.37=start(liley.37.gr))
snps <- merge(snps,DT.37,by.x='id',by.y='id',all.x=TRUE)
## 173 snps don't match after coord conversion
snps <- snps[!is.na(position.37),]


library(annotSnpStats)
snps[,c('a1','a2'):=tstrsplit(alleles,'/')]
snps <- snps[,.(id,chr,position=position.37,a1,a2)]
snps[,pid:=paste(chr,position,sep=':')]
## next add in UK10K stuff so we can align
snp.DT <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab')
M <- merge(snp.DT,snps,by.x='pid',by.y='pid')
alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
alleles <- alleles[!duplicated(pid),]
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
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
M<-M[g.class!='impossible']
## get the list of variants where we need to flip the effect !

## merge by id into original data set so that we can estimate the se using the reference


SHRINKAGE_METHOD<-'ws_emp_shrinkage'
## just the one shrinkage file
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'

shrink.DT <- readRDS(SHRINKAGE_FILE)
shrink.DT<-shrink.DT[,c('pid',shrink=SHRINKAGE_METHOD),with=FALSE]
setkey(shrink.DT,'pid')
pc.emp <- readRDS(BASIS_FILE)

types <- list(z_t1d=5908 + 8825,z_tpo=5908,z_age=5908 + 8825,z_gad=3208,z_ia2=3197,z_pca=2240)


all.res <- lapply(seq_along(types),function(i){
  abtype <- names(types)[i]
  Mt<-merge(M,DT[,.(id=snp,Z=get(`abtype`))],by='id')
  Mt[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
  ## estimate standard error
  #N <- 5908 + 8825
  N <- types[[i]]
  Mt[,se.beta:=sqrt(1/2) * sqrt(1/N) * (sqrt((1/maf) + (1/(1-maf))) *2)]
  Mt[,beta:=Z * se.beta]
  flip.snps <- M[g.class %in% c('rev','revcomp'),]$id
  Mt[id %in% flip.snps, beta:=-1 * beta]
  f.DT <- Mt
  f.DT[,uid:=names(types)[i]]
  ## create_ds_matrix from cupcake expects an odds ratio which it then takes log(OR) we therefore need an
  ## alternative method to create input for projection
  all.DT <- f.DT[!duplicated(pid),.(pid,uid,beta)]
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
  res.DT[trait!='DUMMY:-999',]
}) %>% rbindlist

saveRDS(all.res,file="/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/liley_t1d/proj.RDS")

library(cowplot)
library(ggrepel)

basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x,cat='basis')
all.res[,cat:='liley_t1d']

ggplot(rbind(all.res,basis.DT),aes(x=PC1,y=PC3,label=trait,col=cat)) + geom_point() + geom_text_repel()

## remapping

library(httr)
library(jsonlite)
library(xml2)

server <- "https://rest.ensembl.org"
ext <- "/regulatory/species/homo_sapiens/microarray?"

r <- GET(paste(server, ext, sep = ""), content_type("application/json")) %>% content %>% toJSON %>% fromJSON %>% data.table


saveRDS(rbind(res.DT,basis.DT),file='/home/ob219/rds/rds-cew54-wallace-share/as_basis/westra/westra.RDS')
