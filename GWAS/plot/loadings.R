library(cowplot)
library(GenomicRanges)
BASIS_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS'
basis <- readRDS(BASIS_FILE)
rot.DT <- data.table(pid=rownames(basis$rotation),basis$rotation)
plot.DT <- melt(rot.DT,id.vars='pid')
##sort out x axis position
pos.DT <- data.table(pid=rownames(basis$rotation))[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
pos.DT <- pos.DT[,list(pid=pid,off.pos=pos-min(pos) + 1),by=chr]
pos.DT <- pos.DT[order(chr,off.pos),]
pos.DT[,off.pos:=cumsum(off.pos)]
plot.DT <- merge(plot.DT,pos.DT[,.(chr,pid,off.pos)],by.x='pid',by.y='pid')
pp<-ggplot(plot.DT[abs(value)>0.0000025],aes(x=off.pos,y=abs(value))) + geom_point() + facet_wrap(~variable,ncol=2,dir="v")
save_plot("~/tmp/loadings_man.pdf",pp)



bs.DT <- melt(rot.DT,id.vars='pid')[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
bs.gr <- with(bs.DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,width=1L),value=value,pc=variable))
## there should be none !
mhc.gr <- GRanges(seqnames=Rle(6),ranges=IRanges(start=25e6,end=40e6))
ol <- findOverlaps(bs.gr,mhc.gr) %>% as.matrix

## look to see if enrichment in Javierre et al.
hic.DT <- fread('/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full.txt')
cell.types <- names(hic.DT)[12:28]

ct <- cell.types[1]
all.results <- mclapply(cell.types,function(ct){
  tmp.DT <- hic.DT[get(`ct`)>5,.(chr=baitChr,start=oeStart,end=oeEnd)]
  tmp.gr <- with(tmp.DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end)))
  ## we need to remove MHC as this is not in the basis
  ol <- findOverlaps(tmp.gr,mhc.gr) %>% as.matrix
  tmp.gr <- tmp.gr[-ol[,1],]
  res <- lapply(split(bs.gr,bs.gr$pc),function(gr){
    pc <- unique(gr$pc)
    ol <- findOverlaps(gr,tmp.gr) %>% as.matrix
    tt <- t.test(gr$value[ol[,1]] %>% abs,gr$value[-ol[,1]] %>% abs)
    data.table(pc=pc,p.value=tt$p.value,t.stat=tt$statistic,cell.type=ct)
  }) %>% rbindlist
},mc.cores=8) %>% rbindlist
all.results[,mlp:=-log10(p.value) * sign(t.stat)]

mat.DT <- melt(all.results,id.vars=c('pc','cell.type'),measure.vars='t.stat') %>% dcast(pc~cell.type)
mat <- mat.DT[,-1] %>% as.matrix
rownames(mat) <- mat.DT$pc
library(pheatmap)
pheatmap(mat)

mat.DT <- melt(all.results,id.vars=c('pc','cell.type'),measure.vars='mlp') %>% dcast(pc~cell.type)
mat <- mat.DT[,-1] %>% as.matrix
rownames(mat) <- mat.DT$pc
pheatmap(mat)
