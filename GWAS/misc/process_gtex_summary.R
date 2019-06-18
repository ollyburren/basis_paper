library(optparse)

TEST<-FALSE
option_list = list(
        make_option(c("-f", "--file"), type="character", default=NULL,
              help="chunk file to process", metavar="character")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$file)){
	   print_help(opt_parser)
	    stop("Supply an file to process", call.=FALSE)
    }
}else{
  args <- list(file='/home/ob219/tmp/qstuff/gtex/chunk137.txt')
}




## process and project GTEX data
DATA.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/as_basis/gtex'
if(FALSE){
  fs <- list.files(path=DATA.DIR,pattern='*', full.names=TRUE)
  DT <- data.table(file=basename(fs),ensg=sapply(fs,basename,USE.NAMES=FALSE) %>% gsub("\\..*tsv","",.))
  library(biomaRt)
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
  genedesc <- getBM(attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','entrezgene'), filters = 'ensembl_gene_id', values =DT$ensg, mart =ensembl) %>% data.table
  DT <- merge(genedesc,DT,by.x='ensembl_gene_id',by.y='ensg')
  DT <- DT[gene_biotype=='protein_coding' & !is.na(entrezgene),]
  DT <- DT[!duplicated(ensembl_gene_id),]
  DT[,uid:=paste(external_gene_name,ensembl_gene_id,entrezgene,sep=':')]
  saveRDS(DT,file="/home/ob219/rds/hpc-work/as_basis/gtex/gene_dossier.RDS")
  ## create chunk files
  chunk.gene <- split(DT$file, ceiling(seq_along(DT$file)/25))
  cmds <- sapply(seq_along(chunk.gene),function(i){
      ofile <- file.path('/home/ob219/tmp/qstuff/gtex',sprintf("chunk%d.txt",i))
      write(chunk.gene[[i]],file=ofile)
      sprintf("Rscript /home/ob219/git/as_basis/R/Individual_projection/process_gtex_summary.R -f %s",ofile)
  })
  write(cmds,file="~/tmp/qstuff/gtex_proj.txt")
}
DT <- readRDS("/home/ob219/rds/hpc-work/as_basis/gtex/gene_dossier.RDS")

genes <- scan(args$file,"character")

snp.DT <- fread('/home/ob219/rds/hpc-work/as_basis/gwas_stats/processed_new_aligned_uk10k/snp_manifest/june_10k.tab')


gt.DT <- lapply(genes,function(f){
  tmp <- fread(file.path(DATA.DIR,f))[,file:=f]
}) %>% rbindlist

gtex.SNP <- gt.DT[!duplicated(SNP),.(SNP,seqnames,pos_hg19)]


gtex.SNP[,c('chrom','pos_38','a1','a2','build'):=tstrsplit(SNP,'_')][,c('chrom','pos_38','build'):=list(NULL,NULL,NULL)]
gtex.SNP[,pid:=sub("chr","",seqnames) %>% paste(.,pos_hg19,sep=':')]
M <- merge(snp.DT,gtex.SNP,by.x='pid',by.y='pid')
library(annotSnpStats)
asw <-  getFromNamespace("asw", "annotSnpStats")
x.alleles <- apply(M[,.(ref_a1,ref_a2)],1,paste,collapse="/")
y.alleles <- apply(M[,.(a1,a2)],1,paste,collapse="/")
names(y.alleles)<-names(x.alleles)
message("These are the allele codes as currently defined, before any switching:")
#print(tt <- as.matrix(table(x.alleles, y.alleles)))
## genotype classes
M[,sw.class:=g.class(x.alleles,y.alleles)]
M[,flip:=FALSE]
M[sw.class %in% c('rev','revcomp'),flip:=TRUE]
M <- M[,.(SNP,pid,ref_a1,ref_a2,ref_a1.af,ld.block,flip)]

M <- merge(M,gt.DT,by.x='SNP',by.y='SNP',all.x=TRUE)
## fix the flips
M[flip==TRUE,beta:=beta*-1]

all.DT <- M[!duplicated(paste(file,pid,sep=':')),.(pid,file,beta)]
all.DT <- merge(all.DT,DT[,.(file,uid)],by.x='file',by.y='file')


SHRINKAGE_METHOD<-'ws_emp'
SHRINKAGE_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/shrinkage_june10k.RDS'
BASIS_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/basis_june10k.RDS'

shrink.DT <- readRDS(SHRINKAGE_FILE)
shrink.DT<-shrink.DT[,c('pid',shrink=sprintf("%s_shrinkage",SHRINKAGE_METHOD)),with=FALSE]
setkey(shrink.DT,'pid')
pc.emp <- readRDS(BASIS_FILE)
all.DT <- merge(all.DT,shrink.DT,by.x='pid',by.y='pid')[,shrunk.beta:=beta * ws_emp_shrinkage][,.(pid,uid,shrunk.beta)]

dummy.DT <- snp.DT[,.(pid,uid='DUMMY:-999',shrunk.beta=0)]
all.DT <- rbind(all.DT,dummy.DT)
all.DT <- melt(all.DT,id.vars=c('pid','uid'),measure.vars='shrunk.beta')
setkey(all.DT,'pid')
r.DT <- dcast(all.DT,pid~uid+variable,fill=0)
mat <- as.matrix(r.DT[,-1])
rownames(mat) <- r.DT[[1]]
bc <- predict(pc.emp,newdata=t(mat))
res.DT <- data.table(trait =colnames(mat)  %>% gsub("_shrunk.beta","",.),bc)[trait!='DUMMY:-999',]

OUT.DIR <- '/home/ob219/rds/hpc-work/as_basis/gtex/june_10k/'
saveRDS(res.DT,file=sprintf("%s%s.RDS",OUT.DIR,basename(args$file) %>% gsub("\\.txt","",.)))
message(sprintf("%s%s.RDS",OUT.DIR,basename(args$file) %>% gsub("\\.txt","",.)))


if(FALSE){
  library(cowplot)
  OUT.DIR <- '/home/ob219/rds/hpc-work/as_basis/gtex/june_10k/'
  BASIS_FILE <- '/home/ob219/rds/hpc-work/as_basis/support/basis_june10k.RDS'
  fs <- list.files(path=OUT.DIR,pattern="*.RDS",full.names=TRUE)
  res.DT <- lapply(fs,readRDS) %>% rbindlist
  pc.emp <- readRDS(BASIS_FILE)
  basis.DT <- data.table(trait=rownames(pc.emp$x),pc.emp$x,cat='basis')
  res.DT[,cat:='gtex']
  plot.DT <- rbind(res.DT,basis.DT)
  plot.DT[,c('gname','ensg','entrez'):=tstrsplit(trait,':')]
  plot.DT[,label:=gname]
  plot.DT[is.na(label),label:=trait]
  ggplot(plot.DT,aes(x=PC1,y=PC2,label=label,col=cat)) + geom_point() + geom_text()
  mdt<-melt(res.DT,id.vars='trait',measure.vars=sprintf("PC%d",1:11))
  mdt[,c('mean','sd'):=list(mean(value),sd(value)),by='variable']
  mdt[,Z:=(value-mean)/sd]
  mdt[,p.value:=pnorm(abs(Z)) * 2]
  mdt[,p.adj:=p.adjust(p.value,method="fdr"),by=variable]
  mdt[p.adj<0.05,.(trait,pc=variable,Z,p.value,p.adj)]
  saveRDS(mdt,file="/home/ob219/rds/rds-cew54-wallace-share/as_basis/gtex/results/basis_june10k_protein_coding.RDS")
}
