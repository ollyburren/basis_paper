## Chris Wallace code for Ank Spond GWAS


# Rscript -e 'knitr::spin("gwas-v2.R")'

library(magrittr)
library(ggplot2)
library(snpStats)
library(annotSnpStats)
## this has to be full path otherwise read.plink fails
#if(!file.exists('/home/ob219/share/Data/GWAS/ankspond-feb2019-v2/processed_data/snpStats/ankspond_exclusions_removed.RDS')){
  d <- "/home/ob219/share/Data/GWAS/ankspond-feb2019-v2"
  setwd(d)
  library(data.table)
  ex1 <- fread(file.path(d,"bad.illuminous.clusters.txt"),header=FALSE)
  ex2 <- fread(file.path(d,"exclusions.txt"),header=FALSE)
  x=read.plink(file.path(d,"plink"))
  #str(x)
  x <- annot.plink(x)
  #head(x@snps)
  #dim(x)

  ##' # snp and sample exclusions from TASC
  #head(ex1)
  #head(ex2)
  #tail(ex1)
  #tail(ex2)
  #nrow(ex1)
  #nrow(ex2)
  #table(tolower(ex1$V1) %in% x@snps$snp.name)
  #table(ex2$V1 %in% x@samples$member)

  x <- x[ -which(x@samples$member %in% ex2$V1), -which(x@snps$snp.name %in% tolower(ex1$V1))]
#saveRDS(x,file='/home/ob219/share/Data/GWAS/ankspond-feb2019-v2/processed_data/snpStats/ankspond_exclusions_removed.RDS')
#}else{
#  x <- readRDS('/home/ob219/share/Data/GWAS/ankspond-feb2019-v2/processed_data/snpStats/ankspond_exclusions_removed.RDS')
#}

#dim(x)

##' Quote from TASC paper:
##' > The discovery sample population included 69 Australian, 1,129 British and 983 North American individuals as ankylosing spondylitis cases, and 5,847 controls derived from the 1958 British Birth Cohort genotyped by the Wellcome Trust Case-Control Consortium (n = 1,436), and the Illumina iControlDB database (n = 4,149). HapMap CEU (Utah residents with Northern and Western European ancestry), YRI (Yorubans from Ibadan, Nigeria), JPT (Japanese from Tokyo) and CHB (Chinese from Beijing) samples (n = 262) were added to the controls to use for quality control.
##'
##' We seem to have too few cases. Perhaps because of sample
##' exclusions above?  Should check with Paul @ TASC.

grp <- x@samples$affected
table(grp)

##' generally, alleles are aligned, but a plot of RAF between cases and controls looks odd
XL <- list(x[grp==1,],
           x[grp==2,])
F <- lapply(XL, function(x) col.summary(x)$RAF)  %>% do.call("cbind",.)
cor(F,use="pair")

plot(F,xlab="Cases",ylab="Controls")

##' Drop the really rare SNPs
cs <- col.summary(x)
table(cs$MAF<0.001)
x <- x[,cs$MAF > 0.001] # drops 155039 snps

##' missingness varying by chromosome is no longer a problem :)
summ <- tapply(1:ncol(x), x@snps$chromosome, function(i) col.summary(x[,i]))
m <- sapply(summ,function(a) sum(a$MAF>0))
n <- sapply(summ,function(a) nrow(a))
cbind(levels(x@snps$chromosome),total=n,non.missing=m,non.missing.prop=m/n)

##' # look at sample heterozygosity.
rs <- row.summary(x)
head(rs)
plot(rs$Call.rate,rs$Heterozygosity)

group <- ifelse(rs$Call.rate<0.8,1,2)
##' there are clearly two groups, which correspond to cases and controls - two chips?
table(group,x@samples$affected) # 1787 cases, 5162 controls
## numbers approx match https://www.ncbi.nlm.nih.gov/pubmed/20062062

##' drop ctl-only snps - those with 0 Call rates in cases
ctl <- which(x@samples$affected==1)
cs <- col.summary(x[-ctl,])
hist(cs$Call.rate)
table(called.in.cases=cs$Call.rate>0.5)
x <- x[, which(cs$Call.rate>0.5)]

##' # now look at qc in controls
cs <- col.summary(x[ctl,])
table(cs$Call.rate>0.5)
##' drop additional snps not called in any controls
drop <- cs$Call.rate<0.5

head(cs)
plot(cs$Call.rate,cs$z.HWE)
abline(h=c(-5,5),col="red")

drop <- drop | abs(cs$z.HWE)>5
with(cs[!drop,],plot(Call.rate,z.HWE))
abline(v=c(0.95,0.97,0.98,0.99),col="red")

drop <- drop | is.na(cs$z.HWE) | cs$Call.rate <0.97 | (cs$Call.rate > 0.97 & cs$Call.rate <0.99 & abs(cs$z.HWE)>2)

##' snps remaining
sum(!drop)
x <- x[,which(!drop)]

##' # look at cases at these snps
cs <- col.summary(x[-ctl,])
plot(cs$Call.rate,cs$z.HWE)
abline(h=c(-6,-5,5,6),col="red")
drop <- abs(cs$z.HWE)>6 # will lose some HLA

with(cs[!drop,],plot(Call.rate,z.HWE))
abline(v=c(0.95,0.97,0.98,0.99),col="red")
drop <- drop | is.na(cs$z.HWE) | cs$Call.rate <0.97 | (cs$Call.rate > 0.97 & cs$Call.rate <0.99 & abs(cs$z.HWE)>2)

##' snps remaining
sum(!drop)
x <- x[,which(!drop)]

##' # look at samples by call rate and heterozygosity
par(mfrow=c(1,2))
rs <- row.summary(x[ctl,])
plot(rs$Call.rate,rs$Heterozygosity)
abline(v=c(0.99,0.995),col="red")
abline(h=c(0.32,0.325,0.335,0.34),col="red")
title(main="Controls")
ctlkeep <- rs$Call.rate>0.995 & rs$Heterozygosity >0.320 & rs$Heterozygosity < 0.340

rs <- row.summary(x[-ctl,])
plot(rs$Call.rate,rs$Heterozygosity,main="Cases")
abline(v=c(0.99,0.995),col="red")
abline(h=c(0.31,0.315,0.335,0.34),col="red")
csekeep <- rs$Call.rate>0.995 & rs$Heterozygosity >0.320 & rs$Heterozygosity < 0.340

x <- rbind2(x[ctl,][ctlkeep,], x[-ctl,][csekeep,])

dim(x)
table(x@samples$affected)

##' > A subset of SNPs common to all chip types was then extracted, leaving 301,866 SNPs (from a maximum of 317,502 available).
dim(x)
## https://www.nature.com/articles/ng.513

saveRDS(x,"/home/ob219/share/Data/GWAS/ankspond-feb2019-v2/processed_data/snpStats/ankspond_exclusions_removed.RDS")

#save(x,file="~/as-v2.RData")

##' # First pass analysis: check lambda
y <- x@samples$affected-1
ss <- single.snp.tests(snp.data=sm(x),phenotype=y)
p <- p.value(ss,1)
x2 <- chi.squared(ss,1)
qq.chisq(x2,1)
## lambda is 1.046

df <- x@snps
df$p <- p

w <- which(p<1e-200)
qq.chisq(x2[-w],1)

library(data.table)
df=as.data.table(df)
df[,chromosome:=as.numeric(chromosome)]
df[,position:=as.numeric(position)]
df <- df[order(chromosome,position),]
df[,diff:=c(10000,diff(position))]
df[diff<0,diff:=10000]
df[,x:=cumsum(diff)]

library(ggplot2)
ggplot(df,aes(col=chromosome %% 2, x=x,y=-log10(df$p))) + geom_point()
## looks awful!

##' limit on y to see towers?
##'
ggplot(df,aes(col=chromosome %% 2, x=x,y=-log10(df$p))) + geom_point() + ylim(0,10)

table(df$chromosome) ## only 2805 snps on chromosome 5 ??

################################################################################

##' # read in 1kg data
source("kk.R")

#PLINK="/home/cew54/localc/bin/plink" # plink binary
BCFTOOLS="/home/cew54/localc/bin/bcftools" # bcftools binary
module load plink/2.00-alpha
PLINK2="/home/cew54/localc/bin/plink2"

library(rtracklayer)
library(liftOver)
library(GenomicRanges)

system("gunzip /home/cew54/share/Data/reference/hg18ToHg19.over.chain.gz")
ch = import.chain("/home/cew54/share/Data/reference/hg18ToHg19.over.chain")
system("gzip /home/cew54/share/Data/reference/hg18ToHg19.over.chain")
library(data.table)
gr18 <- with(x@snps,GRanges(seqnames=Rle(paste0("chr",chromosome)),
                            IRanges(position,width=1)))
seqlevelsStyle(gr18) = "UCSC"  # necessary
gr19 <- liftOver(gr18,ch)
gr19
ln <- lengths(gr19)
table(ln)
x <- x[,-which(ln!=1)]
gr19 <- gr19[which(ln==1)]
gr19 <- unlist(gr19)
x@snps$pos19 <- start(gr19)

## read in genetic data
rfile <- tempfile()
tmp <- tempfile()

G <- vector("list",22)

chr=22
samp <- fread("~/share/Data/reference/1000GP_Phase3/1000GP_Phase3.sample")
    samp <- as.data.frame(samp)
    rownames(samp) <- samp$ID
for(chr in 1:22) {
    message(chr)
    fh <- paste0("~/share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr",chr,".hap.gz")
    fl <- paste0("~/share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr",chr,".legend.gz")

#h <- fread(paste("zcat",fh))
l <- fread(paste("zcat",fl))
head(l)
tofind <- as.data.table(subset(x@snps,chromosome==chr))
tofind[,chr:=paste0("chr",chromosome)]
head(tofind)
m <- match(tofind$pos19,l$position)
summary(m)
##' check these are unique
lm <- l[m[!is.na(m)],]
nosnp <- l$position[ l$TYPE!="Biallelic_SNP" ]
bad <- lm$position %in% nosnp
dups <- duplicated(lm$position)
if(any(dups))
    bad <- bad | lm$position %in% lm$position[dups]
m[!is.na(m)][bad] <- NA
lm <- l[m[!is.na(m)],]
fwrite(lm[,.(id)],rfile,sep="\t",col.names=FALSE)
    message("searching for ",nrow(lm)," variants")

## comm <- paste0(BCFTOOLS," view -R ",rfile," ~/share/Data/reference/1000GP_Phase3/bcf/chr",chr,".bcf.gz",
    comm <- paste0(PLINK2," --haps ~/share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr",chr,".hap.gz --legend ~/share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr",chr,".legend.gz ",chr," --extract ",rfile," --min-alleles 2 --max-alleles 2 --min-af 0.01 --make-bed --out ",tmp)
    system(comm)
    g <- read.plink(tmp)

    ## keep only most discriminatory subset of snps
    gg <- g$genotypes
    message("variants found after maf filter: ",ncol(gg))
    rownames(gg) <- rownames(samp)
    ss <- snp.lhs.tests(snp.data=gg,~1,~GROUP,data=samp)
    p <- p.value(ss)
    np <- min(length(p),5000)
    keep <- sample(which(p <= sort(p)[np]))[1:1000]
    gg <- gg[,keep]

    r2 <- ld(gg,gg,stat="R.squared")
    d <- as.dist(1-r2)
    hc <- hclust(d)
    plot(hc)
    ct <- cutree(hc,h=0.9)
    keep <- names(ct)[ !duplicated(ct) ]

    G[[chr]] <- gg[,keep]

}
save(G,file="~/as-1kg.RData")

################################################################################

##' PCA for ancestry
library(snpStats)
library(GUESSFM)
(load("~/as-1kg.RData"))
sapply(G,ncol)

tg <- tag(G[[1]],0.1) # for speed
GG <- G[[1]][,unique(tg@tags)]
dim(GG)

library(purrr)
GG <- purrr::reduce(G,cbind)
dim(GG)


library(data.table)
samp <- fread("~/share/Data/reference/1000GP_Phase3/1000GP_Phase3.sample")
    samp <- as.data.frame(samp)
    rownames(samp) <- samp$ID

## xxmat <- xxt(GG, correct.for.missingg=FALSE)
## evv <- eigen(xxmat, symmetric=TRUE)
## pcs <- evv$vectors[,1:10]
## evals <- evv$values[1:10]

## library(ggplot2)
## colnames(pcs) <- paste0("pc",1:10)
## df <- cbind(samp,pcs)
## ggplot(df,aes(x=pc1,y=pc2,col=GROUP)) + geom_point()
## ggplot(df,aes(x=pc3,y=pc4,col=GROUP)) + geom_point()
## ggplot(df[df$GROUP=="AMR",],aes(x=pc1,y=pc2,col=POP)) + geom_point()
## ggplot(df[df$GROUP=="EUR",],aes(x=pc1,y=pc2,col=POP)) + geom_point()
## evals

## btr <- snp.pre.multiply(GG, diag(1/sqrt(evals)) %*% t(pcs))

## link ids in x to those in GG
(load("~/as-v2.RData"))
library(rtracklayer)
library(liftOver)
library(GenomicRanges)

system("gunzip /home/cew54/share/Data/reference/hg18ToHg19.over.chain.gz")
ch = import.chain("/home/cew54/share/Data/reference/hg18ToHg19.over.chain")
system("gzip /home/cew54/share/Data/reference/hg18ToHg19.over.chain")
library(data.table)
gr18 <- with(x@snps,GRanges(seqnames=Rle(paste0("chr",chromosome)),
                            IRanges(position,width=1)))
seqlevelsStyle(gr18) = "UCSC"  # necessary
gr19 <- liftOver(gr18,ch)
gr19
ln <- lengths(gr19)
table(ln)
x <- x[,-which(ln!=1)]
gr19 <- gr19[which(ln==1)]
gr19 <- unlist(gr19)
x@snps$pos19 <- start(gr19)

library(magrittr)
head(colnames(GG))
matchG <- function(GG) {
    posG <- colnames(GG)  %>% strsplit(.,":")  %>% sapply(., "[[", 2)  %>% as.numeric()
    match(posG,x@snps$pos19)
}
m <- matchG(GG)
summary(m)

library(annotSnpStats)
ag <- strsplit(colnames(GG),":")  %>% sapply(., function(x) paste(x[3],x[4],sep="/"))
ax <- with(x@snps[m,],paste(allele.1,allele.2,sep="/"))
sw <- g.class(ag,ax)
drop <- sw %in% c("ambig","impossible")
if(any(drop)) {
    m <- m[!drop]
    GG <- GG[,-which(drop)]
    sw <- sw[!drop]
}

g <- as(GG,"numeric")
gm <- colMeans(g,na.rm=TRUE)
gsd <- apply(g,2,sd,na.rm=TRUE)
pcs <- prcomp(g,center=gm,scale=gsd)
m <- matchG(g)
xsc <- sm(x)[,m]  %>%  switch.alleles(.,sw %in% c("rev","revcomp")) %>%
  as(.,"numeric") %>% scale(.,center=gm,scale=gsd)
xsc[is.na(xsc)] <- 0
colnames(xsc) <- colnames(g)
df2 <- predict(pcs,newdata=xsc)
df <- as.data.table(cbind(samp,pcs$x[,1:10]))
df2 <- as.data.table(df2[,1:10])
df2$GROUP <- "GWAS"
df2$POP <- x@samples$affected

library(ggplot2)
ggplot(df,aes(x=PC1,y=PC2,col=GROUP)) + geom_point(data=df2) +
  geom_point(pch="+")
ggplot(df,aes(x=PC3,y=PC4,col=GROUP)) + geom_point()

xsc[is.na(xsc)] <- 0
pcs2 <- prcomp(xsc)
ggplot(as.data.table(pcs2$x),aes(x=PC1,y=PC2,col=x@samples$affected)) + geom_point()
ggplot(as.data.table(pcs2$x),aes(x=PC3,y=PC4,col=x@samples$affected)) + geom_point()

km <- kmeans(pcs2$x[,"PC1"],3)
dt <- as.data.table(pcs2$x)
dt$cl <- km$cluster
ggplot(dt,aes(x=PC1,y=PC2,col=factor(cl))) + geom_point()

## test if this helps


##' # First pass analysis: check lambda
y <- x@samples$affected-1
## ss <- single.snp.tests(snp.data=sm(x),phenotype=y)
PC1 <- df2$PC1
PC2 <- df2$PC2
cl <- factor(km$cluster)
ss <- snp.rhs.tests(y ~ strata(cl), snp.data=sm(x))
p <- p.value(ss)
x2 <- chi.squared(ss)
qq.chisq(x2,1)


w <- which(p<1e-200)
qq.chisq(x2[-w],1)

library(data.table)
df <- x@snps
df$p <- p
df=as.data.table(df)
df[,chromosome:=as.numeric(chromosome)]
df[,position:=as.numeric(position)]
df <- df[order(chromosome,position),]
df[,diff:=c(10000,diff(position))]
df[diff<0,diff:=10000]
df[,x:=cumsum(diff)]

library(ggplot2)
ggplot(df,aes(col=chromosome %% 2, x=x,y=-log10(df$p))) + geom_point()
## looks awful!

##' limit on y to see towers?
##'
ggplot(df,aes(col=chromosome %% 2, x=x,y=-log10(df$p))) + geom_point() + ylim(0,20)

table(df$chromosome) ## only 2805 snps on chromosome 5 ??

## limit to basis snps?
b <- fread("/home/cew54/share/as_basis/support_tab/as_basis_snp_support_feb_2018_w_ms.tab")
x@snps$pid <- paste(x@snps$chromosome,x@snps$pos19,sep=":")
table(b$pid %in% x@snps$pid)

use <- x@snps$pid %in% b$pid
qq.chisq(x2[use],1)
ggplot(df[use,],aes(col=chromosome %% 2, x=x,y=-log10(p))) + geom_point() #+ ylim(0,20)

snp

unlink(tmp)
unlink(rfile)

save(GG,file="~/as-1kg.RData")


ss <- snp.rhs.estimates(y ~ 1, snp.data=sm(x))
