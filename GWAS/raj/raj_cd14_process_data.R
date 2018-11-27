library(annotSnpStats)
#library(data.table)
#library(magrittr)

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



DATA_DIR <-  '/home/ob219/share/Projects/twas'
SNP_MANIFEST_FILE <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
GT_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/individual_data/filtered_gt/raj/cd14'
## loads in as X
load(file.path(DATA_DIR,'raj-cd14-genotypes.RData'))
raj.snps <- snps(X) %>% data.table
raj.snps[,snp.idx:=1:nrow(raj.snps)]
## add positions for variants
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
snps<-SNPlocs.Hsapiens.dbSNP144.GRCh37
lu<-snpsById(snps,raj.snps$snp.name,ifnotfound='drop')
lu<-data.table(as.data.frame(lu))
setkey(lu,'RefSNP_id')
setkey(raj.snps,'snp.name')
raj.snps<-raj.snps[lu][order(snp.idx),.(chromosome,snp.name,allele.1,allele.2,pos,strand,alleles_as_ambig,snp.idx)]
raj.snps<-raj.snps[!is.na(pos),pid:=paste(chromosome,pos,sep=':')]
filtered <- X[,raj.snps$snp.idx]
raj.snps <- raj.snps[,.(chromosome,snp.name=pid,allele.1,allele.2,ID=snp.name)]
## add in RAF which is the allele freq wrt to 2
csum <- col.summary(filtered)
raj.snps[,a2.af:=csum$RAF]
snps(filtered) <- as.data.frame(raj.snps)
rownames(snps(filtered)) <- snps(filtered)$ID
## how many variants are missing from this compared to manifest
sm.DT <- fread(SNP_MANIFEST_FILE)
sum(snps(filtered)$snp.name %in% sm.DT$pid)
## some are missing set posterior log or to zero for those
bf <- filtered[,snps(filtered)$snp.name %in% sm.DT$pid]
bf.snps <- snps(bf)
sn <- split(1:nrow(bf.snps),bf.snps$chromosome)

uk10 <- readRDS("/home/ob219/rds/hpc-work/DATA/UK10K/UK10K_0.005_MAF.RDS")
setnames(uk10,'CHROM','CHR')
uk10[CHR=='X',CHR:='23']
uk10[,CHR:=as.numeric(CHR)]
uk10 <- uk10[order(CHR,POS),]
uk10m <- uk10[,.(CHR,BP=POS,uk10_A1=REF,uk10_A2=ALT,uk10_A2_AF=AF)]
uk10m[,pid:=paste(CHR,BP,sep=':')]
ruk10m <- uk10m[pid %in% snps(bf)$snp.name,]

for(chr in names(sn)){
  idx <- sn[[chr]]
  fsnps <- bf[,idx]
  ofile <- sprintf("chr%s.RDS",chr) %>% file.path(GT_OUT_DIR,.)
  sprintf("Writing to %s",ofile) %>% message
  saveRDS(fsnps,file=ofile)
}

files <- list.files(path=GT_OUT_DIR,pattern="*.RDS",full.names=TRUE)



processRAJ <- function(f,out.dir=GT_OUT_DIR){
  message(f)
  G<-readRDS(f)
  ## add chromosome to snps
  snps <- snps(G) %>% data.table
  snps[,posid:=1:.N]
  snps <- merge(snps,ruk10m,by.x='snp.name',by.y='pid')
  snps[,pid:=snp.name]
  ## if there are duplicates this indicate non-binary alleles which we remove
  dup.pid <- snps[duplicated(pid),]$pid
  if(length(dup.pid)>0)
    snps <- snps[!pid %in% dup.pid,]
  Y <- copy(G)
  Y <- Y[,snps$posid]
  y<-snps[,.(a1=uk10_A1,a2=uk10_A2,raf=uk10_A2_AF)]
  Y<-switcheroo(Y,y,do.plot=TRUE)
  ## overwrite with aligned data
  fname <- file.path(out.dir,basename(f))
  saveRDS(Y,file=fname)
  samples <- cbind(samples(Y) %>% data.table,row.summary(Y) %>% data.table)
  snps <- cbind(snps(Y) %>% data.table,col.summary(Y) %>% data.table)
  snps <- merge(snps,ruk10m,by.x='snp.name',by.y='pid')
  fname <- file.path(out.dir,sprintf("summ_%s",basename(f)))
  fname <- gsub("RData","RDS",fname)
  saveRDS(list(snps=snps,samples=samples),file=fname)
  fname
}

for(f in files){
  processRAJ(f)
}
