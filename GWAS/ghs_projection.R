library(RSQLite)
con <- dbConnect(RSQLite::SQLite(), "/home/ob219/rds/rds-cew54-wallace-share/Data/expr/GHS-summary/GHS_Express140110.sqlite")
dbListTables(con)
ghs <- dbReadTable(con,"ghs_snp_on_expression_autosomes") %>% data.table
snp.annotations <- dbReadTable(con,"Affymetrix_6_0_snp_annotations") %>% data.table

## compute overall mean for sets of samples
sum_mean <- function(m,n){
    ((m * n) %>% sum)/sum(n)
}

## compute overall mean for sets of samples
sum_var <- function(s,m,n){
    N <- sum(n)
    S <- ((n-1) * (s^2 + m^2)) %>% sum
    (S - (sum_mean(m,n)^2 * N))/N
}

## compute covariance of x and y
covxy <- function(mx,my,n){
    N <- sum(n)
    sxy <- (n * mx * my) %>% sum
    (sxy/N) - (sum_mean(mx,n)*sum_mean(my,n))
}

comp_beta <- function(mx,my,n){
    N <- sum(n)
    var_x <- sum_var(rep(0,3),mx,n)
    cov_x_y <- covxy(mx,my,n)
    cov_x_y/var_x
}

## this does not work not sure why
comp_se_beta <- function(mx,my,sy,n,beta){
    N <- sum(n)
    var_x <- sum_var(rep(0,3),mx,n)
    var_y <- sum_var(sy,my,n)
    ssx <- var_x * N
    ssy <- var_y * N
    ssr <-  ssy - (beta^2 * ssx)
    (ssr/((N-2)*ssx)) %>% sqrt
}

## this is rederived function that does not rely on computing var_y which
## seems to give negative variances !
comp_se_beta_alt <- function(mx,my,n,beta,r2){
    N <- sum(n)
    var_x <- sum_var(rep(0,3),mx,n)
    cov_x_y <- covxy(mx,my,n)
    ssx <- var_x * N
    ssy <- (cov_x_y * N)^2/(ssx * r2)
    ssr <-  ssy - (beta^2 * ssx)
    (ssr/((N-2)*ssx)) %>% sqrt
}

ghs[,uid:=1:.N]
## slow but need to rewrite functions above to allow for vectorisation speedup
ghs[,beta:=comp_beta(mx=c(0,1,2),my=c(mean_genotype_00,mean_genotype_01,mean_genotype_11),
  n=c(number_genotype_00,number_genotype_01,number_genotype_11)),by=uid]

## this does not work for some reason
if(FALSE){
ghs[,se.beta:=comp_se_beta(mx=c(0,1,2),my=c(mean_genotype_00,mean_genotype_01,mean_genotype_11),
  sy=c(sd_genotype_00,sd_genotype_01,sd_genotype_11),
  n=c(number_genotype_00,number_genotype_01,number_genotype_11),
beta=beta),by=uid]
}

ghs[,se.beta:=comp_se_beta_alt(mx=c(0,1,2),my=c(mean_genotype_00,mean_genotype_01,mean_genotype_11),
  n=c(number_genotype_00,number_genotype_01,number_genotype_11),
beta=beta,r2=association_R_square),by=uid]

ghs[,act.Z:=qnorm(association_p_value/2,lower.tail=FALSE)]
ghs[,calc.Z:=abs(beta/se.beta)]

library(ggplot2)
library(cowplot)

ggplot(ghs,aes(x=act.Z,y=calc.Z)) + geom_point() + geom_abline(col='red')

## this looks OK so next job is to to align and project
## first need to project to build 37

library(rtracklayer)

ghs.36.gr <- with(snp.annotations,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=snp_position,width=1L),id=snp))
c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
ghs.37.gr<-unlist(liftOver(ghs.36.gr,c))
DT.37 <- data.table(id=ghs.37.gr$id,position.37=start(ghs.37.gr))
snps <- merge(snp.annotations,DT.37,by.x='snp',by.y='id',all.x=TRUE)
## 161 snps don't match after coord conversion
snps <- snps[!is.na(position.37),]
snps[,c('tmp.a1','tmp.a2'):=tstrsplit(allele,'/')]
## assume allele 2 is the effect allele.
snps <- snps[,.(id=snp,chr,position=position.37,a1=tmp.a1,a2=tmp.a2)]
snps[,pid:=paste(chr,position,sep=':')]
## next add in UK10K stuff so we can align
snp.DT <- fread('/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab')
m.DT <- merge(snp.DT,snps,by.x='pid',by.y='pid')
m.DT[ref_a1==a1 & ref_a2==a2,flip:=FALSE]
m.DT[ref_a1==a2 & ref_a2==a1,flip:=TRUE]
## these are affy 6 array only 1/4 of these appear to overlap the basis - seems too few !

affy.DT <- fread("~/tmp/GenomeWideSNP_6.na35.annot.csv")
setnames(affy.DT,make.names(names(affy.DT)))
affy.DT <- affy.DT[,.(rs=dbSNP.RS.ID,chr=Chromosome,position=Physical.Position,a1=Allele.A,a2=Allele.B)]
affy.DT[,pid:=paste(chr,position,sep=':')]
## looks like I did things right just very small overlap between platforms.
## one suggestion is to try ImpG

## create a GHS.DT that has info we want
ghs.DT <- ghs[,.(id=snp_name,gene=gene_name,beta,se.beta)]
f.DT <- merge(m.DT,ghs.DT,by.x='id',by.y='id')
f.DT[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
##
f.DT[flip==TRUE,beta:=beta*-1]
f.DT[,uid:=gene]


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
OUT.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/ghs/'
saveRDS(res.DT,file=sprintf("%s%s.RDS",OUT.DIR,'ghs'))

basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x,cat='basis')
res.DT[,cat:='ghs']
mdt<-melt(res.DT,id.vars='trait',measure.vars=sprintf("PC%d",1:11))
mdt[,c('mean','sd'):=list(mean(value),sd(value)),by='variable']
mdt[,Z:=(value-mean)/sd]
mdt[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
mdt[,p.adj:=p.adjust(p.value,method="fdr"),by=variable]
mdt[p.adj<0.01,.(trait,pc=variable,Z,p.value,p.adj)]


ggplot(rbind(res.DT,basis.DT),aes(x=PC1,y=PC2,label=trait,col=cat)) + geom_point() + geom_text() +
coord_cartesian(xlim=c(-0.025,0.025),ylim=c(-0.02,-0.01))
