## plot JIA results for thesis


library(cowplot)
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/25_01_19_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)

main.bc <- list(platelet=c('pdw','mpv','plt'),rbc=c('irf','ret','rdw','hct','mch'),lymphoid=c('mono','baso','eo','neut','lymph'))


ast <- res.DT[trait %in% unlist(main.bc),]
ast[,c('ci.lo','ci.hi'):=list(delta-(sqrt(variance) * 1.96),delta+(sqrt(variance) * 1.96)),]
ast[,variable:=factor(variable,levels=paste0('PC',1:11))]
ast[trait %in% c(main.bc[['lymphoid']],c('plt','ret')),trait:=paste0(trait,'#')]
ast[,trait:=toupper(trait)]
## work out order try hclust
#mat <- melt(ast[,.(pc=variable,trait,delta)],id.vars=c('pc','trait')) %>% dcast(.,trait~pc)
#m <- as.matrix(mat[,-1])
#rownames(m) <- mat$trait
#thc <- dist(m) %>% hclust
#ast[,trait:=factor(trait,levels=thc$labels[thc$order])]
## perhaps better to use main subsets

mbct <- unlist(main.bc,use.names=FALSE)
change <- mbct %in% c(main.bc[['lymphoid']],c('plt','ret'))
mbct[change] <- paste0(mbct[change],'#')
mbct <- toupper(mbct)
ast[,trait:=factor(trait,levels=mbct)]


ast[,is.sig:='']
ast[p.adj<0.05,is.sig:='*']
bbplot <- ggplot(ast[variable!='PC11'],aes(x=variable,y=trait,fill=pmax(delta,-0.007),label=is.sig))  +
geom_tile() + geom_text() + scale_fill_gradient2("Difference\nfrom control") + xlab("Principal Component") +
ylab("Trait") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
save_plot(bbplot,file="~/tmp/ukbbast_13main.pdf",base_height=8,base_aspect_ratio=0.9)
