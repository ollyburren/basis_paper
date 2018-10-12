## uk10k data as our reference data set

## THIS IS FOR REFERENCE ONLY !!!! NOT RUN IN THIS LOCATION

## process PBC

## first convert to annotSnpStats
library(annotSnpStats)
library(data.table)
library(magrittr)
library(optparse)

TEST <- TRUE
DEFAULT_TARGET_DIRNAME <- 'split'

option_list = list(
  make_option(c("-f", "--file"), type="character",default='',
              help="File to process", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)
print(args)



uk10 <- readRDS("/home/ob219/rds/hpc-work/DATA/UK10K/UK10K_0.005_MAF.RDS")
setnames(uk10,'CHROM','CHR')
uk10[CHR=='X',CHR:='23']
uk10[,CHR:=as.numeric(CHR)]
uk10 <- uk10[order(CHR,POS),]
uk10m <- uk10[,.(CHR,BP=POS,uk10_A1=REF,uk10_A2=ALT,uk10_A2_AF=AF)]

switcheroo <- function(x,y,mafdiff=0.1,do.plot=FALSE){
  asw <-  getFromNamespace("asw", "annotSnpStats")
  x.alleles <- apply(x@snps[,alleles(x)],1,paste,collapse="/")
  y.alleles <- apply(y[,.(a1,a2)],1,paste,collapse="/")
  names(y.alleles)<-names(x.alleles)
  message("These are the allele codes as currently defined, before any switching:")
  print(tt <- as.matrix(table(x.alleles, y.alleles)))
  ## genotype classes
  sw.class <- g.class(x.alleles,y.alleles)
  any.comp <- any(sw.class %in% c("comp","revcomp"))
  any.ambig <- any(sw.class=="ambig")
  sw <- sw.class %in% c("rev","revcomp")

  if(length(wh <- which(sw.class=="impossible"))) {
    message(length(wh)," pairwise impossible allele labels found. These genotypes will be set to missing.")
    sw[wh] <- NA
  }
  sw[sw.class=="impossible"] <- NA

  if(any.comp & any.ambig) { # there are reverse complements in the distinguishable cases
    ind <- which(sw.class=="ambig")
    message(length(ind)," SNPs have alleles not completely resolvable without strand information, confirming guess by checking allele freqs.")
    x.cs <- col.summary(x[,ind])
    rdiff <- x.cs[,"RAF"] - y[ind,]$raf
    sw2 <- ifelse(abs(x.cs[,"RAF"] - y[ind,]$raf) < abs(1 - x.cs[,"RAF"] - y[ind,]$raf), FALSE, TRUE)
    too.close <- abs(x.cs[,"MAF"]-0.5)<mafdiff
   # if(any(too.close)) {
   #   can.match <- sw.class %in% c("comp","nochange","rev","revcomp")
   #   xsw <- switch.alleles(x[,-ind],which(sw[-ind]))
   #   ysw <- y[,-ind]
   #   ## step through too.close SNPs checking signed correlation
   #   message("using signed correlation for ",sum(too.close)," SNPs too close to 50% MAF")
   #   ldx <- ld(xsw,x[,ind[too.close],drop=FALSE], stats="R")
   #   ldy <- ld(ysw,y[,ind[too.close],drop=FALSE], stats="R")
   #   ldx[abs(ldx)<0.04] <- NA ## drop uncorrelated - have no information
   #   ldy[abs(ldy)<0.04] <- NA ## drop uncorrelated - have no information
   #   cor.sw <- sapply(1:ncol(ldx), function(j) cor(ldx[,j], ldy[,j], use="pair"))
   #   cor.sw[ abs(cor.sw)<0.8 ] <- NA # NA unless correlation is pretty strong
   #   sw2[too.close] <- cor.sw < 0
   #   too.close <- too.close[is.na(cor.sw)]
   # }
    message(sum(is.na(sw2))," SNPs not resolvable (MAF too close to 0.5).")
    sw[ind] <- sw2
 }

 if(!any.comp & any.ambig) { # there are no reverse complements in distinguishable cases
   ind <- which(sw.class=="ambig")
   message(length(ind)," SNPs have alleles not completely resolvable without strand information,\nbut there is no evidence of strand switches amongst SNPs which are resolvable.\nAssuming fixed strand.")
   ind <- which(sw.class=="ambig")
   sw2 <- x.alleles[ind]==g.rev(y.alleles[ind])
   sw[ind] <- sw2
 }
 ## do switch
 if(any(is.na(sw))){
   x@.Data[,is.na(sw)] <- as.raw("00")
 }


 if(length(wh <- which(sw))) {
     xsw <- asw(x@.Data,wh)
     x@.Data=xsw
     ##x <- switch.alleles(x, wh)
     x@snps[wh, alleles(x) ] <- y[wh, .(a1,a2)]
 }

 if(do.plot) {
   x.cs <- col.summary(x)
   plot(x.cs[,"RAF"], y$raf,main="RAF after switching",xlab="x",ylab="y",pch="+")
   abline(0,1,col="red",lty=1)
 }
 return(x)
}

processJIA <- function(f,out.dir='/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/JIA-2017-data/as_basis'){
  (load(f))
  ## add chromosome to snps
  snps <- snps(G) %>% data.table
  snps[,posid:=1:.N]
  snps <- merge(snps,uk10m,by.x=c('chromosome','position'),by.y=c('CHR','BP'))
  snps[,pid:=paste(chromosome,position,sep=':')]
  ## if there are duplicates this indicate non-binary alleles which we remove
  dup.pid <- snps[duplicated(pid),]$pid
  if(length(dup.pid)>0)
    snps <- snps[!pid %in% dup.pid,]
  Y <- copy(G)
  Y <- Y[,snps$posid]
  y<-snps[,.(a1=uk10_A1,a2=uk10_A2,raf=uk10_A2_AF)]
  Y<-switcheroo(Y,y,do.plot=TRUE)
  fname <- file.path(out.dir,gsub("annotsnpstats[-](.*)","chr\\1",basename(f)))
  save(Y,file=fname)
  samples <- cbind(samples(Y) %>% data.table,row.summary(Y) %>% data.table)
  snps <- cbind(snps(Y) %>% data.table,col.summary(Y) %>% data.table)
  snps <- merge(snps,uk10m,by.x=c('CHR','position'),by.y=c('CHR','BP'))
  fname <- file.path(out.dir,sprintf("summ_%s",gsub("annotsnpstats[-](.*)","chr\\1",basename(f))))
  fname <- gsub("RData","RDS",fname)
  saveRDS(list(snps=snps,samples=samples),file=fname)
  fname
}

# code to do allele switch without requiring y to be annot snp stats
## x is annot snp stats
## y is DT of a1 a2 af_wrt2 - that matches snps in X
## shameless rip off of Chris' code on annotSnpStats

out <- processJIA(args$file)

if(FALSE){
  DATA.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Data/GWAS/JIA-2017-data'
  gt.files <- list.files(path=DATA.DIR,pattern="*.RData",full.names=TRUE)
  RSCRIPT <- '/home/ob219/git/as_basis/R/Individual_projection/switch_alleleles_JIA_q.R'
  cmds <- lapply(gt.files,function(f){
    sprintf("Rscript %s -f %s",RSCRIPT,f)
  }) %>% do.call('c',.)
  write(cmds,file="~/tmp/qsub/jia_align.txt")
  ## run using sh script
}

#summ.files<-lapply(gt.files,processJIA)
