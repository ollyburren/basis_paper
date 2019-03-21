## from projections
pdf("~/tmp/hinks_hla_genetic_correlation.pdf",width=12)
par(mfrow=c(1,2))
## I transposed these from the paper
hinks <- fread("/home/ob219/tmp/hinks_hla_correlation.csv")
hi <- as.matrix(hinks[,-1])
rownames(hi) <- hinks$trait
dist(hi) %>% hclust %>% plot(.,,main="MHC only (Hinks et al.)",sub="",cex=2)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)[grep('^jia_',trait),]
M<-melt(res.DT,id.vars=c('trait','variable'),measure.vars='delta') %>% dcast(.,"trait~variable")
m<-as.matrix(M[,-1])
rownames(m) <- M$trait %>% gsub("^jia_","",.)
dist(m) %>% hclust %>% plot(.,,main="Non-MHC (Basis projection)",sub="",cex=2)
dev.off()
