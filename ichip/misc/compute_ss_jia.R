## Code for computing summary statistics for JIA subtypes.

## parse options first in case there is a problem otherwise slow.
library(optparse)
DEBUG=FALSE
option_list = list(
        make_option(c("-o", "--out"), type="character", default=NULL,
              help="output dir", metavar="character"),
        make_option(c("-p", "--phenotype"), type="character", default='cc',
                    help="Phenotype to use, default cc is just JIA vs Controls", metavar="character"),
        make_option(c("-i", "--input"), type="character", default=NULL,
                                help="Input file (must be of type annotSnpStats)", metavar="character")
        )


## valid phenotypes

PHE<-c('ERA','ext_oligo','jia_undefined','JIA_unknown','jPsA','pers_oligo','RFneg_poly','RFpos_poly','systemic')

if(FALSE){
  OUT.DIR <- '/home/ob219/share/as_basis/ichip/sum_stats/JIA/'
  IN.DIR <- '/home/ob219/share/as_basis/ichip/jia_gt/'
  RSCRIPT <- '/home/ob219/git/basis_paper/ichip/misc/compute_ss_jia.R'
  files <- list.files(path=IN.DIR,pattern="*.RData",full.names=TRUE)
  cmds <- lapply(PHE,function(p){
    lapply(files,function(f){
      sprintf("Rscript %s -o %s -p %s -i %s",RSCRIPT,OUT.DIR,p,f)
    }) %>% do.call('c',.)
  }) %>% do.call('c',.)
  write(cmds,file="~/tmp/qstuff/jia_ichip_gwas.txt")
}

if(!DEBUG){
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)
if (is.null(args$out)){
	print_help(opt_parser)
	stop("At least one argument must be supplied (output file)", call.=FALSE)
}
}else{
  args<-list(
    input='/home/ob219/share/as_basis/ichip/jia_gt/annotSnpStats-22.RData',
    out='/home/ob219/share/as_basis/ichip/sum_stats/JIA/',
    phenotype='ERA'
    )
}
if(!args$phenotype %in% PHE){
  print_help(opt_parser)
	stop(sprintf("%s does not appear to be a vaild phenotype",args$phenotype), call.=FALSE)
}

if(FALSE){
  library(annotSnpStats)
  X <- annot.read.plink("/home/ob219/share/Data/GWAS/IChip/JIA/plink_dataset_files/Immunochip_JIA_PhaseIInew_QCgp_All_SNPQC_UKonly")
  info.DT <- readRDS('/home/ob219/share/as_basis/ichip/ind_sample_info/jia_samples_info.RDS')
  ## split by chromosome
  samp.DT <- samples(X) %>% data.table
  samp.DT[,uid:=1:.N]
  samp.DT <- merge(info.DT[,.(sample_id,ilar_final)],samp.DT,by.x='sample_id',by.y='member')[order(uid),]
  samp.DT[,uid:=NULL]
  samp <- as.data.frame(samp.DT)
  rownames(samp) <- samp$pedigree
  samples(X) <- samp
  OUT.DIR <- '/home/ob219/share/as_basis/ichip/jia_gt'
  snps.DT <- snps(X) %>% data.table
  by.chr.idx <- split(1:nrow(snps.DT),snps.DT$chromosome)
  lapply(names(by.chr.idx),function(chr){
    message(chr)
    idx <- by.chr.idx[[chr]]
    Y <- X[,idx]
    out.file <- file.path(OUT.DIR,sprintf("annotSnpStats-%s.RData",chr))
    save(Y,file=out.file)
  })
}

## for chromosome 1 the largest there are 562166 snps it takes approx 5 seconds to fit glm for 100

library(annotSnpStats)
library(data.table)
## get the data
message(sprintf("Getting data from %s",args$input))
as<-get(load(args$input))
message(sprintf("Finished loading data from %s",args$input))
## create our own column for phenotype
message(sprintf("Filter phenotype %s",args$phenotype))
if(args$phenotype != 'cc'){
  ## phenotype is a bit different as need to remove all other cases as these are not controls
  idx<-which(samples(as)$affected==1 | samples(as)$ilar_final==args$phenotype)
  as<-as[idx,]
}


#system.time(snp.rhs.estimates(phenotype~sex+PC1+PC2+PC3, family='binomial',data=samples(as),sets=1:100 ,snp.data=as(as,"SnpMatrix")))
#for testing
if(DEBUG){
  message("DEBUG mode on")
  idx<-1:50
}else{
  idx<-1:nrow(snps(as))
}
message(sprintf("Fitting models for %d snps",max(idx)))
test=snp.rhs.estimates(affected~sex, family='binomial',data=samples(as),sets=idx,snp.data=as(as,"SnpMatrix"))
message("Reformatting and computing summary statistics")
## convert this to a dataframe and add annotation information as well as maf
tmp<-data.frame(do.call('rbind',test))
tmp$rsid<-rownames(tmp)
for(n in colnames(tmp)){
  tmp[[n]]<-unlist(tmp[[n]])
}
test<-data.table(tmp)
test[,Y.var:=args$phenotype]

## compute the summary stats
controls.idx<-which(samples(as)$affected==0)
cases.idx<-which(samples(as)$affected==1)

cas.ss<-col.summary(as[cases.idx,idx])
ctrl.ss<-col.summary(as[controls.idx,idx])

processSummary<-function(df,prefix){
  rsid<-rownames(df)
  DT<-data.table(df)
  DT<-DT[,.(Calls,MAF)]
  setnames(DT,paste(prefix,names(DT),sep='.'))
  DT$rsid<-rsid
  setkey(DT,rsid)
  DT
}

setkey(test,rsid)

test<-test[processSummary(cas.ss,'case')]
test<-test[processSummary(ctrl.ss,'control')]

## add in allele codings

info<-snps(as[,idx])
info<-data.table(info)
## AF here is with respect to ALT allele
info<-info[,.(snp.name,chromosome,position,allele.1,allele.2)]
setnames(info,c('snp.name','chr','position','a1','a2'))
setkey(info,snp.name)
test<-info[test]
## finally compute Z score and p.value

test[,Z:=beta/sqrt(Var.beta)]
test[,p.val:=2*pnorm(abs(Z),lower.tail=FALSE)]
out.file<-sprintf("%d_%s.RDS",unique(test$chr),args$phenotype)
out.file<-file.path(args$out,out.file)
saveRDS(test,file=out.file)
message(sprintf("Wrote file to %s",out.file))

## package up the results

if(FALSE){
  library(rtracklayer)
  files <- list.files(path="/home/ob219/share/as_basis/ichip/sum_stats/JIA/",pattern="*.RDS",full.names=TRUE)
  by.rtype <- split(files,gsub("[0-9]+\\_([^\\.]+)\\.RDS","\\1",basename(files)))
  lapply(seq_along(by.rtype),function(i){
    trait <- names(by.rtype[i])
    message(trait)
    dat.DT<-lapply(by.rtype[[i]],readRDS) %>% rbindlist
    dat.DT[,id:=1:.N]
    dat.36.gr <- with(dat.DT,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L),id=id))
    c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
    dat.37.gr<-unlist(liftOver(dat.36.gr,c))
    DT.37 <- data.table(id=dat.37.gr$id,position.37=start(dat.37.gr))
    dat.DT <- merge(dat.DT,DT.37,by.x='id',by.y='id',all.x=TRUE)
    ## 173 myo.DT don't match after coord conversion
    dat.DT <- dat.DT[!is.na(position.37),]
    dat.DT[,pid.37:=paste(chr,position.37,sep=":")]

    ## read in snp.manifest as might as well align whilst here
    man.DT <- fread('/home/ob219/share/as_basis/ichip/snp_manifest/ichip_all.tab')
    M <- merge(man.DT,dat.DT,by.x='pid',by.y='pid.37')
    alleles <- data.table(al.x = paste(M$ref_a1,M$ref_a2,sep='/'),al.y=paste(M$a1,M$a2,sep='/'))
    ## to make quick
    align.class <- rep('match',nrow(alleles))
    idx<-which(alleles$al.x!=alleles$al.y)
    x.alleles <- alleles[idx,]$al.x
    names(x.alleles)<-alleles[idx,]$pid
    y.alleles <-  alleles[idx,]$al.y
    names(y.alleles)<-names(x.alleles)
    align.class[idx] <- g.class(x.alleles,y.alleles)
    print(table(align.class))
    M[,g.class:=align.class]
    M[g.class=='comp',c('a1','a2'):=list(ref_a1,ref_a2)]
    sw <- align.class %in% c("rev","revcomp")
    M[,or:=exp(beta)]
    M[sw,c('a1','a2','or'):=list(ref_a1,ref_a2,1/or)]
    ## want to check the the alignment is OK
    # gt.files <- list.files(path='/home/ob219/share/as_basis/ichip/jia_gt',pattern="*.RData",full.names=TRUE)
    # RAF.DT <- lapply(gt.files,function(f){
    #   X <- get(load(f))
    #   cs<-col.summary(X[which(samples(X)$affected==1),])
    #   snps.DT <- snps(X) %>% data.table
    #   snps.DT[,RAF:=cs$RAF]
    #   snps.DT[,pid:=paste(chromosome,position,sep=':')][,.(pid,RAF)]
    # }) %>% rbindlist
    # M[,pid.36:=paste(chr,position,sep=':')]
    # Q <- merge(M,RAF.DT,by.x='pid.36',by.y='pid')
    # sw <- Q$g.class %in% c("rev","revcomp")
    # ## the basis assumes that or is wrt a2, this means that when we flip
    # ## the RAF (a2 frequency) and ref_a1.af are aligned
    # Q[sw,c('a1','a2','or'):=list(ref_a1,ref_a2,1/or)]
    # Q[!sw,RAF:=1-RAF]
    # plot(Q$ref_a1.af,Q$RAF)
    out <- M[,.(pid,a1,a2,or,p.value=p.val)]
    fname <- sprintf("jia_%s.tab",trait)
    fname <- file.path('/home/ob219/share/as_basis/ichip/sum_stats',fname)
    message(fname)
    write.table(out,file=fname,quote=FALSE,sep="\t",row.names=FALSE)
  })
}
