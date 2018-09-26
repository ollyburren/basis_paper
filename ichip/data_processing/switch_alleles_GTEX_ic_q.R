## uk10k data as our reference data set

## process PBC

## first convert to annotSnpStats
library(annotSnpStats)
library(data.table)
library(magrittr)
library(optparse)

TEST <- FALSE

option_list = list(
  make_option(c("-f", "--file"), type="character",default='',
              help="File to process", metavar="character")
)

if(TEST){
  args <- list(file='/home/ob219/rds/rds-cew54-wallace-share/Data/gtex/eur-Whole_Blood-001/chr6.RData')
}else{
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}
print(args)




DATA.DIR<-'/home/ob219/rds/rds-cew54-wallace-share/Data/gtex/eur-Whole_Blood-001'
# uk10 <- readRDS("/home/ob219/rds/hpc-work/DATA/UK10K/UK10K_0.005_MAF.RDS")
# setnames(uk10,'CHROM','CHR')
# uk10[CHR=='X',CHR:='23']
# uk10[,CHR:=as.numeric(CHR)]
# uk10 <- uk10[order(CHR,POS),]
# uk10m <- uk10[,.(CHR,BP=POS,uk10_A1=REF,uk10_A2=ALT,uk10_A2_AF=AF)]
# uk10m[,pid:=paste(CHR,BP,sep=':')]
man.DT<-fread('/home/ob219/share/as_basis/ichip/snp_manifest/ichip_summary_stats.tab')
man.DT[,c('CHR','BP'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
man.DT[,ref_a2.af:=1-ref_a1.af]



# code to do allele switch without requiring y to be annot snp stats
## x is annot snp stats
## y is DT of a1 a2 af_wrt2 - that matches snps in X
## shameless rip off of Chris' code on annotSnpStats

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




processGTex <- function(f,out.dir='/home/ob219/share/as_basis/ichip/individual_gt'){
  (load(f))
  ## get liftover file
  lf <- file.path(DATA.DIR,gsub("\\.RData",".lifthg19",basename(f))) %>% fread
  snps <- snps(X) %>% data.table
  snps[,posid:=1:.N]
  lf <- lf[seqnames==sprintf("chr%s",unique(snps$chromosome)),]
  snps <- merge(snps,lf[,.(POS_37=start,snp.name)],by.x='snp.name',by.y='snp.name')
  snps <- merge(snps,man.DT,by.x=c('chromosome','POS_37'),by.y=c('CHR','BP'))
  snps[,pid:=paste(chromosome,POS_37,sep=':')]
  ## if there are duplicates this indicate non-binary alleles which we remove
  dup.pid <- snps[duplicated(pid),]$pid
  snps <- snps[!pid %in% dup.pid,]
  Y <- copy(X)
  Y <- Y[,snps$posid]
  y<-snps[,.(a1=ref_a1,a2=ref_a2,raf=ref_a2.af)]
  Y<-switcheroo(Y,y,do.plot=FALSE)
  samples <- cbind(samples(Y) %>% data.table,row.summary(Y) %>% data.table)
  snps <- cbind(snps(Y) %>% data.table,col.summary(Y) %>% data.table)[,uid:=1:.N]
  snps <- merge(snps,lf[,.(POS_37=start,snp.name)],by.x='snp.name',by.y='snp.name')
  snps[,POS_38:=position]
  snps[,position:=POS_37]
  snps[,POS_37:=NULL]
  snps <- merge(snps,man.DT,by.x=c('chromosome','position'),by.y=c('CHR','BP'))
  snps <- snps[order(uid),]
  ## now fix position so b37 for annotSnpStats object as makes filtering downstream
  ## easier
  snps.ass <- snps(Y)
  snps.pos <- snps[,.(snp.name,pid,p37=position,uid)]
  merge.snps <- merge(snps.ass,snps.pos,by.x='snp.name',by.y='snp.name') %>% data.table
  merge.snps <- merge.snps[order(uid),]
  merge.snps[,c('old.snp.name','b38_position','position','snp.name'):=list(snp.name,position,p37,pid)]
  merge.snps <- as.data.frame(merge.snps)
  rownames(merge.snps) <- merge.snps$snp.name
  colnames(Y) <- merge.snps$pid
  snps(Y)<-merge.snps
  fname <- file.path(out.dir,basename(f))
  save(Y,file=fname)
  fname <- file.path(out.dir,sprintf("summ_%s",basename(f)))
  fname <- gsub("RData","RDS",fname)
  saveRDS(list(snps=snps,samples=samples),file=fname)
  fname
}

processGTex(args$file)

if(FALSE){
  DATA.DIR<-'/home/ob219/rds/rds-cew54-wallace-share/Data/gtex/eur-Whole_Blood-001'
  gt.files <- list.files(path=DATA.DIR,pattern="*.RData",full.names=TRUE)
  cmds <- lapply(gt.files,function(f){
    sprintf("Rscript /home/ob219/git/as_basis/R/ic_basis/switch_alleles_GTEX_ic_q.R -f %s",f)
  }) %>% do.call('c',.)
  write(cmds,file="~/tmp/qsub/gtex_ic.txt")
}
