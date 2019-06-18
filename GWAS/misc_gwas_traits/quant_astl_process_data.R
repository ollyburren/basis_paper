### this is like quant_astl_process_and_project.R but just stores the cc results file
### allowing us to create a basis from


library(annotSnpStats)
library(optparse)
library(cupcake)

TEST<-FALSE
option_list = list(
        make_option(c("-f", "--file"), type="character", default=NULL,
              help="Astle data file to process", metavar="character")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$file)){
	   print_help(opt_parser)
	    stop("Supply a tsv.gz file to process", call.=FALSE)
    }
}else{
  args <- list(file='/home/ob219/share/Data/GWAS-summary/blood-ukbiobank-2016-12-12/baso_build37_171846_20161212.tsv.gz')
}

SNP_MANIFEST <- '/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab'
OUT_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'


if(FALSE){
  DATA_DIR <- '/home/ob219/share/Data/GWAS-summary/blood-ukbiobank-2016-12-12'
  RSCRIPT <- '/home/ob219/git/basis_paper/GWAS/misc_gwas_traits/quant_astl_process_data.R'
  files <- list.files(path=DATA_DIR,pattern="*.gz$",full.names=TRUE)
  sprintf("Rscript %s -f %s",RSCRIPT,files) %>% write(.,file="~/tmp/qstuff/astle_quant.txt")
}

#files <- list.files(path=DATA_DIR,pattern="*.gz",full.names=TRUE)
man.DT <- fread(SNP_MANIFEST)
#paste(man.DT$pid,as.numeric(gsub("[0-9]+:","",man.DT$pid))+1,sep=':') %>% gsub(":"," ",.) %>% write(.,"~/tmp/gwas_june_tabix.tab")
header<-c('VARIANT','ID','CHR','BP','REF','ALT','ALT_MINOR','DIRECTION','EFFECT','SE','P','MLOG10P','ALT_FREQ','MA_FREQ')
f <- args$file
cmd <- sprintf("/home/ob219/bin/htslib/tabix -R ~/tmp/gwas_june_tabix.tab %s",f)
b.DT <- fread(cmd)
setnames(b.DT,header)
b.DT[,pid:=paste(CHR,BP,sep=':')]
M.tmp <- merge(b.DT,man.DT,by='pid')
M <- M.tmp[!duplicated(pid),.(pid,a1=REF,a2=ALT,beta=EFFECT,p.value=P,ref_a1,ref_a2,ref_a1.af)]
fname <- basename(f)
trait <- gsub("(.*)\\_build37\\_[0-9]+\\_20161212.tsv.gz","\\1",fname)
ss <- as.numeric(gsub("(.*)\\_build37\\_([0-9]+)\\_20161212.tsv.gz","\\2",fname))
alleles <- data.table(pid=M$pid,al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
#alleles <- alleles[!duplicated(pid),]
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
idx<-which(alleles$g.class=='impossible')
if(length(idx) >0){
  M <- M[-idx,]
  alleles <- alleles[-idx,]
}
M <- merge(M,alleles[,.(pid,g.class)],by='pid',all.x=TRUE)
M <- M[!duplicated(pid),]
M <- M[g.class!='match',beta:=beta*-1]
M[,trait:= basename(f) %>% sub("(*.)_build37.*","\\1",.)]
out <- M[,.(pid,a1,a2,or=exp(beta),p.value)]
setkey(out,pid)
fname <- file.path(OUT_DIR,sprintf("astle_quant_%s.tab",trait))
write.table(out,file=fname,sep="\t",row.names=FALSE,quote=FALSE)
sprintf("Wrote %s",fname) %>% message
