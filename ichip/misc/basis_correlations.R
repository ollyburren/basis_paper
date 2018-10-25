gwas.basis <- readRDS('/home/ob219/share/as_basis/GWAS/support/ss_basis_gwas.RDS')
gwas.man.DT <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab')
ichip.basis <- readRDS('/home/ob219/share/as_basis/ichip/support/basis_ic.RDS')
ichip.man.DT <- fread('/home/ob219/share/as_basis/ichip/snp_manifest/ichip_all.tab')

gwas.rot.DT <- data.table(pid=rownames(gwas.basis$rot),gwas.basis$rot) %>% melt(.,id.vars='pid')
setnames(gwas.rot.DT,'value','gwas.proj')
ichip.rot.DT <- data.table(pid=rownames(ichip.basis$rot),ichip.basis$rot) %>% melt(.,id.vars='pid')
setnames(ichip.rot.DT,'value','ic.proj')



br.DT <- merge(gwas.rot.DT,ichip.rot.DT,by.x=c('pid','variable'),by.y=c('pid','variable'))


pids <- intersect(ichip.rot.DT$pid,gwas.rot.DT$pid)
setnames(gwas.man.DT,paste('gwas',names(gwas.man.DT),sep='.'))
setnames(ichip.man.DT,paste('ic',names(ichip.man.DT),sep='.'))
M <- merge(gwas.man.DT,ichip.man.DT,by.x='gwas.pid',by.y='ic.pid')
library(annotSnpStats)
alleles <- data.table(al.x = paste(M$gwas.ref_a1,M$gwas.ref_a2,sep='/'),al.y=paste(M$ic.ref_a1,M$ic.ref_a2,sep='/'))
## to make quick
align.class <- rep('match',nrow(alleles))
idx<-which(alleles$al.x!=alleles$al.y)
x.alleles <- alleles[idx,]$al.x
names(x.alleles)<-alleles[idx,]$pid
y.alleles <-  alleles[idx,]$al.y
names(y.alleles)<-names(x.alleles)
align.class[idx] <- g.class(x.alleles,y.alleles)
print(table(align.class))
M[,class:=align.class]

library(cowplot)
ggplot(M,aes(x=gwas.ref_a1.af,y=ic.ref_a1.af,color=class)) + geom_point()

sw.pid <- M[class %in% c('rev','revcomp'),]$gwas.pid
br.DT[pid %in% sw.pid,ic.proj:=ic.proj*-1 ]
setnames(br.DT,'variable','pc')

gwas <- gwas.basis$rot[M$gwas.pid,]
ic <- ichip.basis$rot[M$gwas.pid,]
icd <- ic
icd[sw.pid,] <- ic[sw.pid,] * -1

c.load <- cor(gwas,icd)

library(pheatmap)
pheatmap(c.load,cluster_rows=FALSE,cluster_cols=FALSE)

basis.DT <- data.table(trait=rownames(gwas.basis$x),gwas.basis$x) %>% melt(.,id.vars='trait')
control.DT <- basis.DT[trait=='control',]
setnames(control.DT,'value','control.loading')
M <- merge(basis.DT,control.DT[,.(variable,control.loading)],by='variable')
M[,variable:=factor(variable,levels=paste0('PC',1:11))]

ggplot(M,aes(x=variable,y=value-control.loading,label=trait,group=trait,color=trait)) + geom_point() + geom_line() + geom_text() 
