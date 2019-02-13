library(data.table)
library(cowplot)

dt <- fread("/home/ob219/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab")

plot.dt <- dt[basis_trait==1 | grepl("^jia_",trait),.(trait,cases)]

plot.dt[,glabel:='Other IMD']
plot.dt[grep("^jia_",trait),glabel:='JIA']
plot.dt[,trait:=gsub("^jia_","",trait)]
plot.dt[,trait:=factor(trait,levels=plot.dt[order(cases,decreasing=TRUE),]$trait)]

plot.dt[trait=='ERA',dlabel:='Enthesitis Related Arthritis']
plot.dt[trait=='PO',dlabel:='Persistant Oligoarthritis']
plot.dt[trait=='EO',dlabel:='Extended Oligoarthritis']
plot.dt[trait=='RFneg',dlabel:='Rheumatoid Factor +']
plot.dt[trait=='RFpos',dlabel:='Rheumatoid Factor -']
plot.dt[trait=='PsA',dlabel:='Psoriatic Arthritis']
plot.dt[trait=='sys',dlabel:='Systemic Arthritis']

plot.dt[trait=='asthma',dlabel:='Asthma']
plot.dt[trait=='CD',dlabel:="Crohns' Disease"]
plot.dt[trait=='RA',dlabel:='Rheumatoid Arthritis']
plot.dt[trait=='T1D',dlabel:='Type 1 Diabetes']
plot.dt[trait=='CEL',dlabel:='Coeliac Disease']
plot.dt[trait=='PBC',dlabel:='Primary Biliary Cholangitis']
plot.dt[trait=='PSC',dlabel:='Primary Sclerosing Cholangitis']
plot.dt[trait=='UC',dlabel:='Ulcerative Colitis']
plot.dt[trait=='MS',dlabel:='Mutiple Sclerosis']
plot.dt[trait=='SLE',dlabel:='Systemic Lupus Erythematosis']
plot.dt[,dlabel:=sprintf("%s (%s)",dlabel,trait)]

plot.dt[,dlabel:=factor(dlabel,levels=plot.dt[order(cases,decreasing=TRUE),]$dlabel)]

pp<-ggplot(plot.dt,aes(x=dlabel,y=cases,fill=glabel)) + geom_bar(stat="identity") + coord_flip() + xlab("Immune-mediated disease") + ylab("Study cases") +
scale_fill_manual("Type",values=c("dodgerblue","firebrick1"))
save_plot(pp,file="~/tmp/case_count.pdf",base_width=7)



plot.dt <- dt[basis_trait==1 | grepl("_lyons$",disease),.(trait,cases)]

plot.dt[,glabel:='Other IMD']

plot.dt[trait=="mpo",c('dlabel','glabel'):=list('VASC:ANCA_MPO+','Lyons EGPA')]
plot.dt[trait=="anca_Neg",c('dlabel','glabel'):=list('EGPA:ANCA -','Lyons EGPA')]
plot.dt[trait=="egpa",c('dlabel','glabel'):=list('EGPA:(All)','Lyons EGPA')]
plot.dt[trait=="mpo_Pos",c('dlabel','glabel'):=list('EGPA:ANCA MPO +','Lyons EGPA')]

plot.dt[,trait:=factor(trait,levels=plot.dt[order(cases,decreasing=TRUE),]$trait)]

plot.dt[trait=='asthma',dlabel:='Asthma']
plot.dt[trait=='CD',dlabel:="Crohns' Disease"]
plot.dt[trait=='RA',dlabel:='Rheumatoid Arthritis']
plot.dt[trait=='T1D',dlabel:='Type 1 Diabetes']
plot.dt[trait=='CEL',dlabel:='Coeliac Disease']
plot.dt[trait=='PBC',dlabel:='Primary Biliary Cholangitis']
plot.dt[trait=='PSC',dlabel:='Primary Sclerosing Cholangitis']
plot.dt[trait=='UC',dlabel:='Ulcerative Colitis']
plot.dt[trait=='MS',dlabel:='Mutiple Sclerosis']
plot.dt[trait=='SLE',dlabel:='Systemic Lupus Erythematosis']
#plot.dt[,dlabel:=sprintf("%s (%s)",dlabel,trait)]

plot.dt[,dlabel:=factor(dlabel,levels=plot.dt[order(cases,decreasing=TRUE),]$dlabel)]

pp<-ggplot(plot.dt,aes(x=dlabel,y=cases,fill=glabel)) + geom_bar(stat="identity") + coord_flip() + xlab("Immune-mediated disease") + ylab("Study cases") +
scale_fill_manual("Type",values=c("dodgerblue","firebrick1"))
save_plot(pp,file="~/tmp/egpa_case_count.pdf",base_width=7)

## cut down version no basis traits

pp<-ggplot(plot.dt[glabel=='Lyons EGPA',],aes(x=dlabel,y=cases)) + geom_bar(stat="identity",fill='dodgerblue') + coord_flip() + xlab("Immune-mediated disease") + ylab("Study cases")
save_plot(pp,file="~/tmp/egpa_case_count_cut.pdf",base_width=7)
