## script to process eqtl data from DICE database from

## to start process one file

library(optparse)

TEST<-TRUE
option_list = list(
        make_option(c("-f", "--fname"), type="character", default=NULL,
              help="index of phenotype to process ", metavar="character")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$fname)){
	   print_help(opt_parser)
	    stop("Supply a vcf file to process", call.=FALSE)
    }
}else{
  args <- list(fname='/home/ob219/share/Data/GWAS-summary/DICE/CD4_STIM.vcf')
}


VCF_FILE <- args$fname

vcf.DT <- fread(VCF_FILE)
info.cols <- c('Gene','GeneSymbol','Pvalue','Beta')
vcf.DT[,eval(`info.cols`):=tstrsplit(INFO,';')]
for(f in info.cols){
  message(f)
  vcf.DT[,eval(`f`):=gsub("^[^=]+=","",get(`f`))]
}

vcf.DT[,INFO:=NULL]
setnames(vcf.DT,'#CHROM',"CHR")
vcf.DT[,CHR:=gsub("chr","",CHR)]
vcf.DT <- vcf.DT[CHR!='X',]
vcf.DT[,pid:=paste(CHR,POS,sep=':')]
SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
snp.DT <- fread(SNP_MANIFEST)
vcf.DT <- vcf.DT[pid %in% snp.DT$pid,]

basis.DT <- vcf.DT[,.(pid,a1=REF,a2=ALT,or=Beta,p.value=Pvalue,ensg=Gene)]

M <- merge(basis.DT,snp.DT,by='pid')
M[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
M <- M[!duplicated(pid),]

# this to get a derived beta that takes into account differing MAF's and sample sizes.


library(annotSnpStats)
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
## for basis a2 is the risk allele this means flipped alleles are OK but
## matched alleles are the wrong way round
#M[g.class=='rev',der.beta:=der.beta*-1]

M.out <- M[,.(pid,a1=ref_a1,a2=ref_a2,or,p.value,ensg)]
by.gene<-split(M.out,M.out$ensg)
ft <- basename(VCF_FILE) %>% sub("\\.vcf","",.)
OUT_DIR <- file.path('/home/ob219/rds/rds-cew54-wallace-share/as_basis/GWAS/sum_stats/DICE/',ft)
dir.create(OUT_DIR)
for(n in names(by.gene)){
  message(n)
  out <- by.gene[[n]][,.(pid,a1,a2,or,p.value)]
  fname <- file.path(OUT_DIR,sprintf("%s.tab",n))
  write.table(out,file=fname,quote=FALSE,row.names=FALSE,sep="\t")
}
