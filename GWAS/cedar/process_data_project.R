library(data.table)
library(annotSnpStats)
library(optparse)

TEST<-TRUE
option_list = list(
        make_option(c("-f", "--fname"), type="character", default=NULL,
              help="CEDAR file to process", metavar="character")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$fname)){
	   print_help(opt_parser)
	    stop("Supply an file containing a list of genes to process", call.=FALSE)
    }
}else{
  args <- list(fname='/home/ob219/share/as_basis/cedar_eqtl/cedar-CD14-chr1.csv.gz')
}




SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_0619.RDS'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_0619.RDS'
OUT_DIR <- "/home/ob219/share/as_basis/GWAS/cedar/projections"
SRC_OUT_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr'
DATA.DIR <- '/home/ob219/share/as_basis/cedar_eqtl'

man.DT <- fread(SNP_MANIFEST)

sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
pc.emp <- readRDS(BASIS_FILE)


#filez <- list.files(path=DATA.DIR,pattern="*.gz",full.names=TRUE)

#f <- filez[1]
f <- args$fname
sprintf("Processing %s",f) %>% message

## grab cell type
ct <- basename(f) %>% gsub("cedar-([^-]+)-.*","\\1",.)
DT <- sprintf("zcat %s",f) %>% fread
DT[,c('chr','pos','a1','a2'):=tstrsplit(pid,'_')]
DT[,pid:=paste(chr,pos,sep=':')]
by.probe <- split(DT[,.(beta,se.beta=sqrt(vbeta),pid,a1,a2)],DT$probe)
## this code if we need to align but from the looks of things they are already aligned therefore can be much simpler
#library(parallel)
res <- lapply(names(by.probe),function(n){
  probe <- by.probe[[n]]
  otrait <- sprintf("CEDAR:%s:%s",ct,n)
  probe[,trait:=otrait]
  #stat.DT <- probe[,c('Z','p.value',trait):=list(beta/se.beta,2* pnorm(abs(beta/se.beta),lower.tail=FALSE))]
  message(otrait)
  #M <- merge(stat.DT[pid %in% man.DT$pid,.(trait=otrait,pid,a1,a2,or=exp(beta),p.value)],man.DT,by='pid',all.y=TRUE)
  #idx <- which(is.na(M$or))
  #sprintf("%d missing",length(idx)) %>% message
  #if(length(idx)!=0)
  #M[idx,c('trait','a1','a2','or'):=list(otrait,ref_a1,ref_a2,or=1)]
  #setkey(M,pid)
  tmp <- merge(probe,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp$ws_emp_shrinkage * tmp$beta
  #pfile <- file.path(SRC_OUT_DIR,sprintf("%s_source.RDS",tra))
  #saveRDS(tmp[,.(pid,or,p.value,ws_emp_shrinkage)],file=pfile)
  tmp[is.na(metric),c('metric','trait'):=list(0,otrait)]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  predict(pc.emp,newdata=mat.emp)
})

res <- do.call('rbind',res)
ofile <- basename(f) %>% gsub(".csv.gz",".RDS",.) %>% file.path(OUT_DIR,.)
saveRDS(res,file=ofile)
