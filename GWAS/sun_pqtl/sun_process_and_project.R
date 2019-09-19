library(annotSnpStats)

## process pQTL data

library(optparse)

TEST<-FALSE
option_list = list(
        make_option(c("-f", "--file"), type="character", default=NULL,
              help="Sun directory listing file to process", metavar="character")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$file)){
	   print_help(opt_parser)
	    stop("Supply a directory listing file to process", call.=FALSE)
    }
}else{
  args <- list(file='/home/ob219/tmp/qstuff/pqtl/chunk1.txt')
}

OUT_DIR <- '/home/ob219/share/as_basis/GWAS/sun_pqtl/13_traits_0919_unfiltered/'

if(FALSE){
  OUT_DIR <- '/home/ob219/share/as_basis/GWAS/sun_pqtl/13_traits_0919_unfiltered/'
  PQTL_DIR <- '/home/ob219/share/as_basis/sun_pqtl/gwas_basis_june10k_pqtl'
  ## remove dirs that we have already processed
  all.dirs <- list.dirs(path=PQTL_DIR,recursive = FALSE)
  done <- list.files(path=OUT_DIR,pattern="*.RDS") %>% gsub("\\.RDS","",.) %>% file.path(PQTL_DIR,.)
  all.dirs <- all.dirs[!all.dirs %in% done]
  ## process 50 at a time
  spdir <- split(all.dirs,ceiling(seq_along(all.dirs)/50))
  for(i in seq_along(spdir)){
    ofile <- sprintf("/home/ob219/tmp/qstuff/pqtl/chunk%s.txt",i)
    write(spdir[[i]],file=ofile)
  }
  cfiles <- sprintf("/home/ob219/tmp/qstuff/pqtl/chunk%s.txt",seq_along(spdir))
  SCRIPT <- "/home/ob219/git/basis_paper/GWAS/sun_pqtl/sun_process_and_project.R"
  sprintf("Rscript %s -f %s",SCRIPT,cfiles) %>% write(.,"~/tmp/qstuff/sun_pqtl.txt")
}

processPQTL <- function(dir){
  sprintf("Processing %s",dir) %>% message
  ofile <- file.path(OUT_DIR,sprintf("%s.RDS",basename(dir)))
  files <- list.files(path=dir,pattern="*.gz",full.names=TRUE)
  p.DT <- lapply(files,function(f){
    sprintf("zcat %s",f) %>%  fread
  }) %>% rbindlist
  setnames(p.DT,c('snpid','chr','pos','a1','a2','effect','se','lp'))
  SNP_MANIFEST <-'/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_13_traits_0919.tab'
  man.DT <- fread(SNP_MANIFEST)
  p.DT[,pid:=paste(chr,pos,sep=':')]
  p.DT <- p.DT[,.(pid,a1=toupper(a1),a2=toupper(a2),or=exp(effect),p.value=exp(lp))]
  M <- merge(p.DT,man.DT,by='pid')
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
  ## note here effect is wrt to allele 1
  M <- M[g.class=='match',or:=1/or]
  M[,trait:= basename(dir)]

  SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  ## where snp is missing make it zero
  tmp[is.na(metric),metric:=0]
  tmp[,trait:= basename(dir)]
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas_13_traits_0919.RDS'
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  ofile <- file.path(OUT_DIR,sprintf("%s.RDS",M$trait[1]))
  saveRDS(all.proj,file=ofile)
}

#dirlist <- scan(args$file,"character") %>% head(.,n=5)
dirlist <- scan(args$file,"character")
for(d in dirlist){
  ofile <- file.path(OUT_DIR,sprintf("%s.RDS",basename(d)))
  if(file.exists(ofile)){
    sprintf("Directory %s already exists skipping",ofile) %>% message
  }else{
    tryCatch({
      processPQTL(d)
    }, error=function(e){sprintf("Error %s with %s",e,ofile) %>% message})
  }
}
