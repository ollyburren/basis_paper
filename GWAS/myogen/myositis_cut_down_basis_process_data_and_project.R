## bipolar prognosis GWAS
library(annotSnpStats)
library(rtracklayer)

myo.DT <- fread("zcat ~/share/Data/GWAS-summary/MYOGEN/Sep2018_summary_meta.txt.gz")
snps.DT <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_myositis.tab')
myo.DT[,pid:=paste(CHR,BP,sep=':')]
myo.DT <- myo.DT[!pid %in% myo.DT[duplicated(pid),],]
myo.DT[,id:=1:.N]
setnames(myo.DT,make.names(names(myo.DT)))

myo.36.gr <- with(myo.DT,GRanges(seqnames=Rle(paste0('chr',CHR)),ranges=IRanges(start=BP,width=1L),id=id))
c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
myo.37.gr<-unlist(liftOver(myo.36.gr,c))
DT.37 <- data.table(id=myo.37.gr$id,position.37=start(myo.37.gr))
myo.DT <- merge(myo.DT,DT.37,by.x='id',by.y='id',all.x=TRUE)
## 173 myo.DT don't match after coord conversion
myo.DT <- myo.DT[!is.na(position.37),]
myo.DT[,pid.37:=paste(CHR,position.37,sep=":")]
myo.DT[,c('a1','a2'):=list(A1_risk,ifelse(A1_risk!=A1_minor,A1_minor,A2_major))]
out <- myo.DT[,.(pid,a1,a2,or=OR,p.value=P)]
M <- merge(snps.DT,myo.DT,by.y='pid.37',by.x='pid')
M <- M[,.(pid,ref_a1,ref_a2,ref_a1.af,A1_minor,MAF_jdm=MAF_dmjdmpm,P,OR,a1,a2)]
M[is.na(OR),c('A1_minor','MAF_jdm','P','OR','a1','a2','missing'):=list(ifelse(ref_a1.af<0.5,ref_a1,ref_a2),
  ref_a1.af,
  0.99,1,a1=ref_a1,a2=ref_a2,TRUE)]

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
# M[sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,1/OR,1-MAF_jdm)]
# M[!sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,OR,MAF_jdm)]

## note that here effect is wrt to a1 and not allele two - the basis assumes
## that effect is wrt to allele2 therefore revcomp and rev are OK
M[!sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,1/OR,1-MAF_jdm)]
M[sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,OR,MAF_jdm)]
M[align.class=='impossible',c('OR','P'):=list(NA,NA)]


SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_myositis_gwas.RDS'
sDT <- readRDS(SHRINKAGE_FILE)
stmp<-sDT[,.(pid,ws_emp_shrinkage)]
setkey(M,pid)
tmp <- merge(M,stmp,by='pid',all.y=TRUE)
tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
## where snp is missing make it zero
tmp[is.na(metric),metric:=0]
#tmp[,trait:= 'myositis_myogen']
tmp[,trait:= 'myositis_myogen']
saveRDS(tmp,file='/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis_cut_down_source.RDS')
B <- dcast(tmp,pid ~ trait,value.var='metric')
snames <- B[,1]$pid
mat.emp <- as.matrix(B[,-1]) %>% t()
colnames(mat.emp) <- snames
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_myositis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
stop("Something wrong basis and projection matrix don't match")
all.proj <- predict(pc.emp,newdata=mat.emp)
#saveRDS(all.proj,file='/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis.RDS')
saveRDS(all.proj,file='/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis_cut_down.RDS')

## process subset data

for (trait in c('jdm','pm','dm')){
  message(trait)
  meta <- fread(sprintf("~/share/Data/GWAS-summary/MYOGEN/Dec2018/gwas_%s.meta",trait))
  fq <- fread(sprintf("~/share/Data/GWAS-summary/MYOGEN/Sep2018/%s.frq",trait))
  myo.DT <- merge(meta[,.(CHR,BP,SNP,A1,P,OR)],fq[,.(SNP,fq.A1=A1,fq.A2=A2,MAF)],by='SNP')
  myo.DT[,A2:=ifelse(A1==fq.A1,fq.A2,fq.A1)]
  snps.DT <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_myositis.tab')
  myo.DT[,pid:=paste(CHR,BP,sep=':')]
  myo.DT <- myo.DT[!pid %in% myo.DT[duplicated(pid),],]
  myo.DT[,id:=1:.N]
  setnames(myo.DT,make.names(names(myo.DT)))
  myo.36.gr <- with(myo.DT,GRanges(seqnames=Rle(paste0('chr',CHR)),ranges=IRanges(start=BP,width=1L),id=id))
  c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
  myo.37.gr<-unlist(liftOver(myo.36.gr,c))
  DT.37 <- data.table(id=myo.37.gr$id,position.37=start(myo.37.gr))
  myo.DT <- merge(myo.DT,DT.37,by.x='id',by.y='id',all.x=TRUE)
  ## 173 myo.DT don't match after coord conversion
  myo.DT <- myo.DT[!is.na(position.37),]
  myo.DT[,pid.37:=paste(CHR,position.37,sep=":")]
  myo.DT[,c('a1','a2'):=list(A1,A2),]
  #out <- myo.DT[,.(pid,a1,a2,or=OR,p.value=P)]
  M <- merge(snps.DT,myo.DT,by.y='pid.37',by.x='pid')
  M <- M[,.(pid,ref_a1,ref_a2,ref_a1.af,MAF,P,OR,a1,a2)]
  #M[is.na(OR),c('MAF_jdm','P','OR','a1','a2','missing'):=list(ifelse(ref_a1.af<0.5,ref_a1,ref_a2),
  #  ref_a1.af,
  #  0.99,1,a1=ref_a1,a2=ref_a2,TRUE)]
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
  # M[sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,1/OR,1-MAF_jdm)]
  # M[!sw,c('a1','a2','or','MAF_jdm'):=list(ref_a1,ref_a2,OR,MAF_jdm)]
  ## note that here effect is wrt to a1 and not allele two - the basis assumes
  ## that effect is wrt to allele2 therefore revcomp and rev are OK
  M[!sw,c('a1','a2','or'):=list(ref_a1,ref_a2,1/OR)]
  M[sw,c('a1','a2','or'):=list(ref_a1,ref_a2,OR)]
  M[align.class=='impossible',c('OR','P'):=list(NA,NA)]
  SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_myositis_gwas.RDS'
  sDT <- readRDS(SHRINKAGE_FILE)
  stmp<-sDT[,.(pid,ws_emp_shrinkage)]
  setkey(M,pid)
  tmp <- merge(M,stmp,by='pid',all.y=TRUE)
  tmp$metric <- tmp[['ws_emp_shrinkage']] * log(tmp$or)
  ## where snp is missing make it zero
  tmp[is.na(metric),metric:=0]
  #tmp[,trait:= 'myositis_myogen']
  tmp[,trait:= sprintf("%s_myogen",trait)]
  of <- sprintf('/home/ob219/share/as_basis/GWAS/myogen_myositis/%s_myositis_cut_down_source.RDS',trait)
  saveRDS(tmp,file=of)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  mat.emp <- as.matrix(B[,-1]) %>% t()
  colnames(mat.emp) <- snames
  BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_myositis_gwas.RDS'
  pc.emp <- readRDS(BASIS_FILE)
  if(!identical(colnames(mat.emp),rownames(pc.emp$rotation)))
  stop("Something wrong basis and projection matrix don't match")
  all.proj <- predict(pc.emp,newdata=mat.emp)
  saveRDS(all.proj,file=sprintf('/home/ob219/share/as_basis/GWAS/myogen_myositis/%s_myositis_cut_down.RDS',trait))
}

## collate and project

myogen <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis_cut_down.RDS')
myogen <- data.table(trait=rownames(myogen),myogen)
myogen <- melt(myogen,id.var='trait')
myogen[,c('n0','n1'):=list(4724,1711)]
myogen[,category:='myogen']

jdm <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/jdm_myositis_cut_down.RDS')
jdm <- data.table(trait=rownames(jdm),jdm)
jdm <- melt(jdm,id.var='trait')
jdm[,c('n0','n1'):=list(4724,473)]
jdm[,category:='myogen']

dm <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/dm_myositis_cut_down.RDS')
dm <- data.table(trait=rownames(dm),dm)
dm <- melt(dm,id.var='trait')
dm[,c('n0','n1'):=list(4724,705)]
dm[,category:='myogen']

pm <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/pm_myositis_cut_down.RDS')
pm <- data.table(trait=rownames(pm),pm)
pm <- melt(pm,id.var='trait')
pm[,c('n0','n1'):=list(4724,533)]
pm[,category:='myogen']

all.proj <- rbindlist(list(myogen,jdm,dm,pm))
all.proj[,n:=n1+n0]

## load in basis and variance
VARIANCE_FILE <- '/home/ob219/share/as_basis/GWAS/support/myositis_ss_av_june.RDS'
var.DT <- readRDS(VARIANCE_FILE)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_myositis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,value:=(value-control.loading)]
basis.DT[,Z:=sign(value)]
basis.DT[,p.adj:=1]
basis.DT[,control.loading:=NULL]

## compute the variance of a projection
all.DT <- merge(all.proj,var.DT,by.x='variable',by.y='pc')
all.DT <- merge(all.DT,control.DT,by.x='variable',by.y='PC')
all.DT[,variance:=((log(n)-(log(n1) + log(n-n1)))+ log(mfactor)) %>% exp]
## add in control loading
all.DT[,Z:=(value-control.loading)/sqrt(variance)]
all.DT[,p.value:=pnorm(abs(Z),lower.tail=FALSE) * 2]
all.DT[,p.adj:=p.adjust(p.value,method="fdr"),by='variable']
all.DT[,delta:=value-control.loading]


## plot

res.DT <- all.DT

all.traits <- traits<-split(res.DT$trait,res.DT$category) %>% lapply(.,unique)

## add in basis traits for comparison

pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
basis.DT[,category:='zzz_basis']




## work out which PC's we want to draw forest plots for
#PCs <- res.DT[trait %in% FOCUS.DISEASES & variable!='PC11' & p.adj<SIG.THRESH,]$variable %>% unique

#pc <- 'PC3'

forest_plot_focal <- function(proj.dat,basis.dat=basis.DT,pc,focal,title,fdr_thresh=0.05,theme=NA){
  dat <- proj.dat[variable==pc & (p.adj<fdr_thresh | trait %in% focal),]
  dat[trait %in% focal,category:='aa_focal_diseases']
  dat[,ci:=1.96 * sqrt(variance)]
  dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
  dat <- rbind(basis.DT[variable==pc & trait!='control',],dat,fill=TRUE)
  dat[,trait:=factor(trait,levels=dat[order(category,delta,decreasing=TRUE),]$trait)]
  #ggplot(dat,aes(x=trait,y=delta,colour=category)) + geom_point(aes(size=log10(n1))) + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  #coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(pc)
  if(is.na(theme)){
    theme <- theme(panel.grid.major.y=element_line(colour='lightgrey',linetype=3))
  }
  ggplot(dat,aes(x=trait,y=delta,colour=category,alpha=p.adj<fdr_thresh)) + geom_point() + geom_errorbar(aes(ymin=lower,ymax=upper)) +
  coord_flip() + geom_hline(yintercept=0,col='red',linetype=2) + ggtitle(title) + theme +
  xlab("Trait") + ylab("Change in basis loading from control") + guides(alpha=FALSE) + scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0.3))
}

## for talk leave out some of the stuff to make figure clearer

forest_plot_focal(res.DT,pc='PC1',focal=all.traits[['myogen']],title="JIA subtypes PC1")


#fdr stuff
myo <- readRDS('/home/ob219/share/as_basis/GWAS/myogen_myositis/myogen_myositis_cut_down_source.RDS')
n.snps <- 1
## compute prior
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_myositis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
## get the rotations
rot.dt <- data.table(pid=rownames(pc.emp$rotation),pc.emp$rotation)
rot.dt <- melt(rot.dt,id.vars='pid')
rot.dt <- merge(myo[,.(pid,P,or,shrink=ws_emp_shrinkage)],rot.dt,by='pid')
##
rot.dt[,abs_s_v:=abs(shrink * value)]
roty.dt<-merge(rot.dt,rot.dt[,list(sabs_s_v=sum(abs_s_v)),by='variable'],by='variable')
roty.dt[,nprior_h1:=abs_s_v/sabs_s_v]
roty.dt[,nprior_h0:=1-nprior_h1]
## prior is proportional to shrinkage * rotation for a given PC
## need to normalise this so prior adds up to one
#rot.dt[,nprior:=n.snps * ((shrink * value)/sum(shrink * value)),by=variable]
rot.dt[,nprior_h1:=abs(shrink * value)/sum(abs(shrink * value)) * n.snps,by=variable]
rot.dt[,nprior_h0:=1-nprior_h1]
## number of p.values in the set that are less than the current one divided by total
#P(p<a) ??
rot.dt[,plta:=rank(P)/.N,by=variable]
## fdr should be (a * pi)/P(p<a) ??
rot.dt[,fdr:=(P * nprior_h0)/plta]
rot.dt <- rot.dt[order(variable,fdr,decreasing=TRUE),]
rot.dt[, head(.SD, 3), by=variable][order(variable),][,.(variable,pid,P,nprior_h0,plta,fdr)]


f <- ecdf(rot.dt$P)
rot.dt[,chris.fdr:=P/f(P),by='variable']
rot.dt[, head(.SD, 3), by=variable][order(variable),][,.(variable,pid,P,nprior_h0,plta,fdr,chris.fdr)]

rot.dt[,tail(.SD,3),by=variable][,.(pid,P,fdr)]

## multiply the element and the weight
library(cowplot)
## compare between basis
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_myositis_gwas.RDS'
pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basis.DT <- merge(basis.DT,control.DT,by.x='variable',by.y='PC')
basis.DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
basis.DT[,category:='zzz_basis']

BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'

pc.emp <- readRDS(BASIS_FILE)
control.DT <- data.table(PC=names( pc.emp$x["control",]),control.loading= pc.emp$x["control",])
basisold..DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x) %>% melt(.,id.vars='trait')
basisold..DT <- merge(basisold..DT,control.DT,by.x='variable',by.y='PC')
basisold..DT[,delta:=(value-control.loading)]
#basis.DT[,c('lower','upper'):=list(delta,delta)]
basisold..DT[,category:='basis_old']

M <- merge(basisold..DT[,.(pc=variable,trait,old.load=value)],basis.DT[,.(pc=variable,trait,myo.load=value)],by=c('pc','trait'))

library(ggrepel)

M[,pc:=factor(pc,levels=paste('PC',1:11,sep=''))]
## need this as a supp figure
pp <- ggplot(M,aes(x=old.load,y=myo.load,label=trait)) + geom_point() + geom_text_repel() + facet_wrap(~pc,scales="free") +
xlab("GWAS basis") + ylab("Myogen SNP only basis")

save_plot(pp,file="~/tmp/myo_pc_scores.pdf",base_height=10)
