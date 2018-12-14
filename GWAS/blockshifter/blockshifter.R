library(data.table)
library(GenomicRanges)
library(parallel)

## this is the analysis for the data freeze


## FUNCTIONS



superblock.gapsize<-5

## function to remove granges that overlap a given grange

remove.ol<-function(gr,r.gr){
	ol<-as.matrix(findOverlaps(gr,r.gr))
	if(nrow(ol)==0)
		return(gr)
	gr[-ol[,1],]
}

## a function to read in GWAS bed files quickly

pmifile2GRanges<-function(fname){
  dt<-fread(fname,sep="\t",header=TRUE,stringsAsFactors=FALSE)
  dt<subset(dt,chr!="")
	with(dt,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start+1,end=end),name=rsid,p=pval,ppi=ppi))
}

## a function to provide all possible rotations of a vector - returns a matrix

rotate<-function(v){
	l=length(v)
	vc<-l
	m<-matrix(0,nrow=l,ncol=l)
	m[,1]<-v
	for(i in 2:l){
		pm<-m[,i-1]
		m[,i]<-c(pm[l],pm[-l])
	}
	return(m)
}





blockshifter <-
  function(contacts,phenpc.gr,test.set,control.set,perm.no = 1e5,test.set.label,control.set.label,super.target.only =
             FALSE,target.only = TRUE) {

  #if(missing(test.set.label))
   # test.set.label=test.set

  #if(missing(control.set.label))
  #  control.set.label=control.set

  ## if user wants only cell type specific interactions across whole dataset
  if(super.target.only){
    ft<-contacts[,16:length(names(contacts)),with=FALSE]
    contacts<-contacts[rowSums(ft)==1,]
  }

  message(paste(test.set,"vs",control.set))
  test.set<-unlist(strsplit(test.set,","))
  control.set<-unlist(strsplit(control.set,","))
  ## define cell types to analyse and convert to a logical field
  ti<-c(test.set,control.set)
  keep.cols<-c(10:13,which(names(contacts) %in% ti))
  cf<-contacts[,keep.cols,with=FALSE]
  ## as we are interested in targets we get duplicates and also depending
  ##on tissues selected we have targets that have no interactions we remove
  ## these
  # no interactions
  cf<-cf[which(rowSums(cf[,ti,with=FALSE])>0),]
  setkey(cf,oeID)
  ## get coord details
  det<-unique(cf)[,1:4,with=FALSE]
  ## this bit creates a table combining duplicate targets and if 1 has overlap with
  ## cell types of interest then combined object has this set to true
  cf<-cf[, lapply(.SD,function(i) sum(i)>=1),by=oeID,.SDcols=which(names(cf) %in% ti)]
  cf<-merge(cf,det,by.x='oeID',by.y='oeID')
  ## create a granges object
  det.gr<-with(cf,GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,end=oeEnd)))
  mcols(det.gr)<-cf[,c('oeID',ti),with=FALSE]
  ## next remove any targets that don't overlap at least one gwas hit
  det.gr<-subsetByOverlaps(det.gr,phenpc.gr)


  ## compute control an test set
  det.gr$control<-rowSums(as.data.frame(mcols(det.gr[,control.set])))>0
  det.gr$test<-rowSums(as.data.frame(mcols(det.gr[,test.set])))>0
  overlap.idx<-which(det.gr$test & det.gr$control)
  if(target.only){
    message(paste("We have ",length(overlap.idx),"instances where test targets overlap control targets removing .."))
    if(length(overlap.idx)>0)
      det.gr<-det.gr[-overlap.idx,]
  }else{
    message(paste("We have ",length(overlap.idx),"instances where test targets overlap control targets reassign to control .."))
    if(length(overlap.idx)>0)
      det.gr[overlap.idx,]$test<-FALSE
  }

  det.gr$test<-!det.gr$control


  ##next we make sure that det.gr is sorted

  det.gr<-det.gr[order(as.character(seqnames(det.gr)),start(det.gr)),]

  ## define super blocks which are blocks of targets separated by no more than one target

  ## split by chromosome
  det.gr.chr<-split(det.gr,as.character(seqnames(det.gr)))

  ## for each chromosome
  res<-unlist(GRangesList(lapply(det.gr.chr,function(ct){
    ## define the difference between adjacent oeID's
    #ct$bb<-c(diff(ct$oeID),3)
    ct$bb<-c(diff(ct$oeID),10 + superblock.gapsize)
    chr<-as.character(seqnames(ct))
    sid<-1
    bc<-1
    res.list<-list()
    ## get a list of all index positions where
    ## adjacent difference is greater than two these are the positions of breaks
    ## we can then create labels accordingly; we use a for loop as we want to
    ## keep track of variables
    for(i in which(ct$bb>superblock.gapsize + 1)){
      res.list[[bc]]<-rep(bc,i-sid)
      sid<-i
      bc<-bc + 1
    }
    ## we then create label sb that is a mix of chromosome and the label
    ## as between lapply calls we will get the same label without chr prefix
    ct$sb<-paste(chr,c(1,do.call("c",res.list)),sep=":")
    ct
  })))



  ## work out which SNPs go with which target region
  ol<-as.matrix(findOverlaps(phenpc.gr,res))
  phenpc.gr$detID<-0
  phenpc.gr[ol[,1],]$detID<-ol[,2]
  phenpc.gr$metric = phenpc.gr$loading

  ## filter phenpc.gr as for some metric we will have NA ie ppi where MAF is not known in 1K genome.

  phenpc.gr<-phenpc.gr[!is.na(phenpc.gr$metric),]

  ## create a data.table object that matches targets to gwas -lop(p)
  dt<-data.table(metric=phenpc.gr$metric,id=phenpc.gr$detID,key="id")

  ## create a summary object that computes slp, the sum of all -log(p) in a target
  ## and the number of snps overlapping that target object
  dt<-as.data.frame(dt[,list(smetric=sum(metric),n=length(metric)),by=id])

  ## remove ranges from res and recast as a data.frame
  res.df<-as.data.frame(mcols(res))
  ## add an index so we can merge
  res.df$id<-1:nrow(res.df)

  ## merge meta data above so we have a merge object that contains
  ## slp for a target and also whether that target is in test or control set
  mer<-merge(res.df,dt,by.x='id',by.y='id',all.x=TRUE)

  ## which targets have no SNPs overlapping ?
  na.idx<-which(is.na(mer$smetric))
  ## set these to 0 rather than NA
  if(length(na.idx)>0){
    ## should not get here !
    mer[na.idx,]$smetric<-0
    mer[na.idx,]$n<-0
  }
  ## work out two sets of super blocks mixed (i.e. contain test and control targets)
  ## nonmixed contain just test or just control targets
  dt<-data.table(mer,key="sb")
  dt<-dt[,list(n=length(control),control=sum(control)),by=sb]
  mer$mixed<-mer$sb %in% dt[-(which(dt$n==dt$control | dt$control==0)),]$sb
  nonmix<-mer[!mer$mixed,]

  ## [BEGIN PROCESSING NON MIXED SUPERBLOCKS]

  ## here if we decide latterly to allow fragments withing superblocks to
  ## permute then we need to change the key to oeID
  ## below

  dt<-data.table(nonmix,key="sb")
  dt<-dt[,list(smetric=sum(smetric),n=sum(n),test=sum(test)>0),by=sb]

  ## next remove all where slp is na

  #dt<-dt[slp!=0,]
  ## probability of a test is
  p.test<-sum(dt$test)/nrow(dt)
  p.test<-c(p.test,1-p.test)
  mon.perm<-dt[,list(perms=list(p1=c(smetric,0,n,0),p2=c(0,smetric,0,n))),by=sb]

  perms<-do.call("rbind",mon.perm$perms)
  ## to get a set of perms we do
  n.perms=nrow(perms)
  n.sets=n.perms/2

  nomix.perms<-mclapply(1:perm.no,function(i){
    offsets<-seq(1,n.perms,by=2)
    rand.idx<-sample(0:1,n.sets,replace=TRUE,prob=p.test) + offsets
    gh<-perms[rand.idx,]
    colSums(gh)
  },mc.cores=8)

  ## [END OF PROCESSING NON MIXED SUPERBLOCKS]

  ## [BEGIN PROCESSING MIXED SUPERBLOCKS]

  mixed<-mer[mer$mixed,]

  ## leave this in in case removing targets with no coverage was bad
  ## there are some superblocks with no coverage what are these ?
  #dt<-data.table(mixed,key="sb")
  #dt<-dt[,list(empty=sum(slp)==0),by=sb]
  #mixed<-mixed[-(which(mixed$sb %in% dt[dt$empty]$sb)),]

  ## remove factor on sb
  mixed$sb<-as.character(mixed$sb)
  block.idx<-split(seq_len(nrow(mixed)),mixed$sb)
  mixed<-cbind(as.numeric(mixed$control),as.numeric(mixed$smetric),as.numeric(mixed$n))

  ## create every possible permutation for mixed blocks

  mixed<-lapply(block.idx,function(i){
    m<-mixed[i,]
    rv<-rotate(m[,1])
    rma<-matrix(0,nrow=nrow(m),ncol=4)
    for(i in 1:nrow(m)){
      t.idx<-which(rv[,i]==0)
      rma[i,1]<-sum(m[t.idx,2])
      rma[i,3]<-sum(m[t.idx,3])
      rma[i,2]<-sum(m[-t.idx,2])
      rma[i,4]<-sum(m[-t.idx,3])
    }
    rma
  })


  mixedt<-do.call("rbind",mixed)

  ## create a permutation index that allows us to rapidly sample and extract perms for each
  ## mixed super block
  offsets<-do.call("c",lapply(seq_along(block.idx),function(i) rep(i,length(block.idx[[i]]))))
  offsets<-which(c(1,diff(offsets))==1)
  pindex<-cbind(offsets,c(diff(offsets)-1,nrow(mixedt)-offsets[length(offsets)]))

  ## quicker but will use a lot of memory !!
  ## we can offset this perhaps by splitting number of permutations up
  pindex<-do.call("rbind",mclapply(pindex[,2],function(y) sample(0:y,perm.no,replace=TRUE),mc.cores=8)) + pindex[,1]
  mix.perms<-mclapply(seq(ncol(pindex)),function(i){
    colSums(mixedt[pindex[,i],])
  },mc.cores=8)


  ## [END PROCESSING MIXED SUPERBLOCKS]

  ## [COMPUTE NULL DELTA'S]


  null.delta<-sapply(seq_along(mix.perms),function(i){
    mixp<-mix.perms[[i]]
    nomixp<-nomix.perms[[i]]
    tmlp<-(mixp[1]+nomixp[1])/(mixp[3]+nomixp[3])
    cmlp<-(mixp[2]+nomixp[2])/(mixp[4]+nomixp[4])
    tmlp-cmlp
  })



  ## what is the actual delta ?


  act.tmlp<-sum(mer[mer$test,]$smetric)/sum(mer[mer$test,]$n)
  act.cmlp<-sum(mer[mer$control,]$smetric)/sum(mer[mer$control,]$n)
  delta<-act.tmlp-act.cmlp

  test.lab<-paste(test.set,sep=",",collapse=",")
  control.lab<-paste(control.set,sep=",",collapse=",")
  if(super.target.only){
    message("Cell type unique targets only considered")
    ptitle<-paste('Cell type unique targets',test.lab,'VS',control.lab)
    prefix<-'cu'
  }
  if(!super.target.only & target.only){
    message("Set unique targets only considered")
    ptitle<-paste('Set unique targets',test.lab,'VS',control.lab)
    prefix<-'su'
  }
  if(!super.target.only & !target.only){
    message("Test set unique targets")
    ptitle<-paste('Test set unique targets',test.lab,'VS',control.lab)
    prefix<-'tu'
  }
  #gwas_trait<-sub("\\.pmi","",basename(gwas.file))
  #sink(output.file)
  #print(paste0("#####",gwas_trait,"#####"))
  print(paste("TEST mlp:",signif(act.tmlp,digits=3),"n.snps:",sum(mer[mer$test,]$n)))
  print(paste("CONTROL mlp:",signif(act.cmlp,digits=3),"n.snps:",sum(mer[mer$control,]$n)))
  print(paste("Delta:",signif(delta,digits=3)))
  pval.emp<-sum(delta<null.delta)/length(null.delta)
  pval.emp.twotail<-sum(abs(delta)<abs(null.delta))/length(null.delta)

  if(pval.emp==0)
    p.val.emp<-1/perm.no

  ## also compute Z score
  z<-(delta-mean(null.delta))/sd(null.delta)
  pval.z<-2*pnorm(abs(z),lower.tail = FALSE)
  output.df<-data.table(test=test.set.label,control=control.set.label,perm=perm.no,p.emp=pval.emp,p.emp.twotail=pval.emp.twotail,z=z,p.val.z=pval.z,delta=delta)
  output.df
}

stop()
## VARIABLES

#super.target.only implies target.only this does the comparison looking at targets that are unique to a cell type
#target.only this does the comparison looking at targets that are unique to test or control sets
#if both are false then this does the comparison looking at test targets that don't overlaps control targets
## if a target falls in both test and control then we reassign the target to be a control

#chic.thresh<-5
#super.target.only=FALSE ## only consider cell type target overlaps that are unique to cell type
#target.only=TRUE ## after subsetting by test and control sets only consider targets uniquely overlap control sets or test.sets

## read in contact information

## contacts_file = path to R data object containing the CHiCAGO data
## test_tissue = tissue to form test set
## control_tissue = tissue to form control set
## gwas_file = path to bed file containing gwas file to integrate (probably in PMI output format i.e. bed)
## perms = number of permutations
## output_file = file to output results to
## metric = can be one of p (raw pvals),pp (posterior prob),norm.p (normalised pvals)
## test_tissue_label - label for test tissue (if blank is name of test_tissue)
## control_tissue_label - label for control tissue (if blank is name of control_tissue)



PCHIC_FILE <- '/home/ob219/share/Projects/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab'
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
PC <- 'PC1'

contacts<-fread(PCHIC_FILE)
chic.thresh=5
## Need to make sure that we only consider unique contacts (i.e. don't double count where a bait overlaps more than one TSS)
## not sure that this is a problem given that we collapse oeID's later but do this just in case.

contacts$uid<-with(contacts,paste(baitID,oeID,sep=":"))
setkey(contacts,uid)
contacts<-unique(contacts)
contacts$uid<-NULL
## blockshifter operates on thresholded interactions
ft<-contacts[,16:length(names(contacts)),with=FALSE]
for(n in names(ft)){
	contacts[[n]]<-contacts[[n]]>=chic.thresh
}

pc.emp <- readRDS(BASIS_FILE)
# pc.DT <- data.table(variant=rownames(pc.emp$rotation),loading=pc.emp$rotation[,PC])
# pc.DT[,c('chr','pos'):=tstrsplit(variant,':') %>% lapply(.,as.numeric)]
# phenpc.gr <- with(pc.DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,width=1L),loading=loading))
#
# mhc.gr<-GRanges(seqnames=Rle('6'),ranges=IRanges(start=25e6,end=35e6))
# phenpc.gr<-remove.ol(phenpc.gr,mhc.gr)
# #save(null.delta,file=sub("txt","RData",output.file))
# tissues<-names(contacts)[16:length(names(contacts))]
#
# bs<-mclapply(tissues,function(ot){
#   lapply(tissues,function(it){
#     if(ot != it){
#       test.set.label=ot
#       control.set.label=it
#       res.DT <- blockshifter(
#         contacts = contacts,phenpc.gr = phenpc.gr,test.set = ot,control.set =it,
#         test.set.label = test.set.label ,control.set.label = control.set.label
#       )
#       res.DT[,pc:=PC]
#      }
#   }) %>% rbindlist
# },mc.cores=8) %>% rbindlist
## gets kicked off the queue ! ignore
#
#
#
# get_hclust_order <- function(DT){
#   deltas <- melt(DT,id.vars=c('test','control'),measure.vars='z') %>% dcast(.,test~control+variable)
#   delta.mat <- as.matrix(deltas[,-1])
#   rownames(delta.mat) <- deltas$test
#   diag(delta.mat) <- 0
#   delta.hc <- dist(delta.mat) %>% hclust
#   delta.hc$labels[delta.hc$order]
# }
#
# #bs.back <- copy(bs)
#
# bs[,test:=factor(test,levels=get_hclust_order(bs))]
# bs[,control:=factor(control,levels=get_hclust_order(bs))]
#
# ggplot(bs,aes(x=test,y=control,fill=z)) + geom_tile() + scale_fill_gradientn("Z",colours=c('green','white','red')) +
# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

library(cowplot)
bs_test_set <- function(test.set,control.set,tlabel,clabel){
  test.set <- paste(test.set,sep=',',collapse=',')
  control.set <- paste(control.set,sep=',',collapse=',')
  mclapply(paste('PC',1:10,sep=''),function(PC){
    pc.DT <- data.table(variant=rownames(pc.emp$rotation),loading=pc.emp$rotation[,PC])
    pc.DT[,c('chr','pos'):=tstrsplit(variant,':') %>% lapply(.,as.numeric)]
    phenpc.gr <- with(pc.DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,width=1L),loading=loading))
    mhc.gr<-GRanges(seqnames=Rle('6'),ranges=IRanges(start=25e6,end=35e6))
    phenpc.gr<-remove.ol(phenpc.gr,mhc.gr)
    bs.ret <- blockshifter(
      contacts = contacts,phenpc.gr = phenpc.gr,test.set = test.set,control.set = control.set,
      test.set.label = tlabel ,control.set.label = clabel
    )
    bs.ret[,pc:=PC]
    bs.ret
  },mc.cores=8) %>% rbindlist
}

## do lymphoid vs myeloid

#test.set <- c('Naive_B','Total_B','Foetal_thymus','Total_CD8','Total_CD4_MF','Total_CD4_Activated','Naive_CD4','Naive_CD8','Total_CD4_NonActivated')
#control.set <- tissues[!tissues %in% lymphoid ]

test.set <- c('Total_CD4_Activated','Total_CD4_NonActivated')
control.set <- c('Megakaryocytes','Erythroblasts')
## what about b cells vs t cells

test1 <- bs_test_set(test.set,control.set,"a_nact_t","megery")

test.set <- c('Total_CD4_Activated')
control.set <- c('Total_CD4_NonActivated')

test2 <- bs_test_set(test.set,control.set,"act_t","nact_t")
#control.set <- c('Total_CD8','Total_CD4_MF','Total_CD4_Activated','Naive_CD4','Naive_CD8','Total_CD4_NonActivated')


plot.dt <- rbind(test1[,.(pc,comp=paste(control,test,sep='_vs_'),z)],
test2[,.(pc,comp=paste(control,test,sep='_vs_'),z)]) %>% melt(.,id.vars=c('pc','comp')) %>% dcast(.,pc~comp)
library(ggrepel)
ggplot(plot.dt,aes(x=megery_vs_a_nact_t,y=nact_t_vs_act_t,label=pc)) + geom_point() + geom_text_repel() +
geom_hline(yintercept=0) + geom_vline(xintercept=0)
