
## can we load in the whole shooting match

# this command used to work but now seems to time out !
#all.1kg <- fread("bcftools view -H -q 0.01 /home/ob219/reference_panels/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz",sep="\t")
all.1kg <- fread("/home/ob219/reference_panels/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.0.01MAF.txt",sep="\t")
all.1kg[,pid:=paste(V1,V2,sep=":")]
### this reads in a map file and spits out myogen Z score file for use with ImpG
library(rtracklayer)
library(annotSnpStats)

files <- c('renton_mg' = 'MyastheniaGravis_Renton_JAMA_Neurol_2015.txt.gz',
'renton_mg_late' = 'MyastheniaGravis_LateOnset_Renton_JAMA_Neurol_2015.txt.gz',
'renton_mg_early' = 'MyastheniaGravis_YoungOnset_Renton_JAMA_Neurol_2015.txt.gz')
for(f in names(files)){
  fname <- file.path('/home/ob219/share/Data/GWAS-summary',files[f])
  mg.DT <- sprintf("zcat %s",fname) %>% fread
  mg.DT[,pid:=paste(CHR,BP,sep=':')]
  mg.DT <- mg.DT[!pid %in% mg.DT[duplicated(pid),],]
  mg.DT[,id:=1:.N]
  setnames(mg.DT,make.names(names(mg.DT)))
  n0 <- max(mg.DT$NCONTROLS)
  n1 <- max(mg.DT$NCASES)
  # note that joining by SNP name is worse than joining by position
  #merged <- merge(mg.DT,all.1kg[,.(kg.CHR=V1,kg.BP=V2,SNP=V3,REF=V4,ALT=V5)],by='SNP',all.x=TRUE)
  #missing 1289

  merged.pid <- merge(mg.DT,all.1kg[,.(kg.CHR=V1,kg.BP=V2,kg.SNP=V3,REF=V4,ALT=V5,pid)],by='pid')
  ## missing 598333
  merged.pid<- merged.pid[!is.na(kg.SNP),.(CHR,BP,SNP=MARKER,A1=toupper(ALLELE1),A2=toupper(ALLELE2),P=PVALUE,beta=EFFECT_ALLELE1,kg.SNP,REF,ALT,N=n1+n0)]
  ## for impss it seems as if effect allele is a2
  ## looks as if a1 is the effect allele align so that OR 1kg ref allele is a1 and alternative is a2
  alleles <- data.table(al.x = paste(merged.pid$REF,merged.pid$ALT,sep='/'),al.y=paste(merged.pid$A1,merged.pid$A2,sep='/'))
  align.class <- rep('match',nrow(alleles))
  idx<-which(alleles$al.x!=alleles$al.y)
  x.alleles <- alleles[idx,]$al.x
  names(x.alleles)<-alleles[idx,]$pid
  y.alleles <-  alleles[idx,]$al.y
  names(y.alleles)<-names(x.alleles)
  align.class[idx] <- g.class(x.alleles,y.alleles)
  print(table(align.class))
  merged.pid[,g.class:=align.class]
  ## remove impossible and ambig as cannot align these
  merged.pid <- merged.pid[!g.class %in% c('impossible','ambig')]
  ## effect allele is a1 if we have to flip it means that effect allele is allele 2 which we want
  ## thus we flip aligned alleles to get what we want rather than those that don't match
  sw <- merged.pid$g.class %in% c("comp","match")
  merged.pid[sw,beta:=beta * -1]
  mg.filt <- merged.pid[,.(MarkerName=kg.SNP,Allele1=REF,Allele2=ALT,Z=qnorm(P/2,lower.tail=FALSE) * sign(beta),n=N,P,OR=exp(beta))]
  ofile <- sprintf("/home/ob219/rds/hpc-work/ssimp/%s_with_or.txt",f)
  message(ofile)
  write.table(mg.filt[!is.na(Z),],file=ofile,row.names=FALSE,sep=" ",quote=FALSE)
}

files <- sprintf("%s_with_or.txt",names(files))
all.cmds <- lapply(files,function(f){
  cmd_template <- '/home/ob219/git/ssimp_software/compiled/ssimp-linux-0.5.6 --gwas /home/ob219/rds/hpc-work/ssimp/%s --ref 1KG/EUR --out /home/ob219/rds/hpc-work/ssimp/mg/chr%s_%s_imputed.txt --impute.maf 0.01 --impute.range %s'
  trait <- basename(f) %>% gsub(".txt","",.)
  sapply(1:22,function(chr){
    sprintf(cmd_template,f,chr,trait,chr)
  })
}) %>% do.call('c',.)

## NOTE THAT CHR17 WILL NOT WORK! NEED TO DO THE THING BELOW FOR THIS
write(all.cmds,'~/tmp/qstuff/mg_impute.txt')

## annoyingly chr17 does not work as it ssimp guesses the wrong genome build
## create specific files that remove variants withing the first 1Mb of chr17
## to make ssimp pick the correct genome build!
f<-names(files)[1]
for(f in names(files)){
  fname <- file.path('/home/ob219/share/Data/GWAS-summary',files[f])
  mg.DT <- sprintf("zcat %s",fname) %>% fread
  mg.DT <- mg.DT[CHR==17,]
  ## removing SNPs in the fist 1MB region seems to help
  mg.DT <- mg.DT[BP>1e6,]
  mg.DT[,pid:=paste(CHR,BP,sep=':')]
  mg.DT <- mg.DT[!pid %in% mg.DT[duplicated(pid),],]
  mg.DT[,id:=1:.N]
  setnames(mg.DT,make.names(names(mg.DT)))
  n0 <- max(mg.DT$NCONTROLS)
  n1 <- max(mg.DT$NCASES)
  # note that joining by SNP name is worse than joining by position
  #merged <- merge(mg.DT,all.1kg[,.(kg.CHR=V1,kg.BP=V2,SNP=V3,REF=V4,ALT=V5)],by='SNP',all.x=TRUE)
  #missing 1289

  merged.pid <- merge(mg.DT,all.1kg[,.(kg.CHR=V1,kg.BP=V2,kg.SNP=V3,REF=V4,ALT=V5,pid)],by='pid')
  ## missing 598333
  merged.pid<- merged.pid[!is.na(kg.SNP),.(CHR,BP,SNP=MARKER,A1=toupper(ALLELE1),A2=toupper(ALLELE2),P=PVALUE,beta=EFFECT_ALLELE1,kg.SNP,REF,ALT,N=n1+n0)]
  ## for impss it seems as if effect allele is a2
  ## looks as if a1 is the effect allele align so that OR 1kg ref allele is a1 and alternative is a2
  alleles <- data.table(al.x = paste(merged.pid$REF,merged.pid$ALT,sep='/'),al.y=paste(merged.pid$A1,merged.pid$A2,sep='/'))
  align.class <- rep('match',nrow(alleles))
  idx<-which(alleles$al.x!=alleles$al.y)
  x.alleles <- alleles[idx,]$al.x
  names(x.alleles)<-alleles[idx,]$pid
  y.alleles <-  alleles[idx,]$al.y
  names(y.alleles)<-names(x.alleles)
  align.class[idx] <- g.class(x.alleles,y.alleles)
  print(table(align.class))
  merged.pid[,g.class:=align.class]
  ## remove impossible and ambig as cannot align these
  merged.pid <- merged.pid[!g.class %in% c('impossible','ambig')]
  ## effect allele is a1 if we have to flip it means that effect allele is allele 2 which we want
  ## thus we flip aligned alleles to get what we want rather than those that don't match
  sw <- merged.pid$g.class %in% c("comp","match")
  merged.pid[sw,beta:=beta * -1]
  mg.filt <- merged.pid[,.(MarkerName=kg.SNP,Allele1=REF,Allele2=ALT,Z=qnorm(P/2,lower.tail=FALSE) * sign(beta),n=N,P,OR=exp(beta))]
  ofile <- sprintf("/home/ob219/rds/hpc-work/ssimp/%s_17_with_or.txt",f)
  message(ofile)
  write.table(mg.filt[!is.na(Z),],file=ofile,row.names=FALSE,sep=" ",quote=FALSE)
}

files <- sprintf("%s_17_with_or.txt",names(files))
all.cmds <- lapply(files,function(f){
  cmd_template <- '/home/ob219/git/ssimp_software/compiled/ssimp-linux-0.5.6 --gwas /home/ob219/rds/hpc-work/ssimp/%s --ref 1KG/EUR --out /home/ob219/rds/hpc-work/ssimp/mg/chr17_%s_imputed.txt --impute.maf 0.01'
  trait <- basename(f) %>% gsub(".txt","",.)
  trait <- trait %>% gsub("\\_17","",.)
  sprintf(cmd_template,f,trait)
}) %>% do.call('c',.)
write(all.cmds,'~/tmp/qstuff/mg_17_impute.txt')



#~/git/slurmer/qlines_csd3.rb -t 23:59:00 mg_17_impute.txt

#NEED ALL THIS FOR IT TO WORK
#export LC_ALL=C; unset LANGUAGE
if(FALSE){
#/home/ob219/git/ssimp_software/compiled/ssimp-linux-0.5.6 --gwas /home/ob219/rds/hpc-work/ssimp/mg_IgNeg.txt --ref 1KG/EUR --out /home/ob219/rds/hpc-work/ssimp/mg_IgNeg_imputed.txt --impute.maf 0.01
## check imputed vs input ?
library(parallel)
 all.files <- list.files(path='/home/ob219/rds/hpc-work/ssimp/mg',pattern='*.mg_with_or_imputed.txt',full.names=TRUE)
 imp <- mclapply(all.files,function(f){
   fread(f,select=c('chr','pos','SNP','source','r2.pred','z_imp','Allele1','Allele2','maf','N.imp'))
 },mc.cores=8) %>% rbindlist

 input <- fread("/home/ob219/rds/hpc-work/ssimp/renton_mg_with_or.txt")

 M<-merge(input[,.(SNP=MarkerName,orig.z=Z,in.a1=Allele1,in.a2=Allele2,OR)],imp[,.(SNP,source,r2.pred,z_imp,out.a1=Allele1,out.a2=Allele2,maf,N.imp)],by='SNP')
 alleles <- data.table(al.x = paste(M$out.a1,M$out.a2,sep='/'),al.y=paste(M$in.a1,M$in.a2,sep='/'))
 align.class <- rep('match',nrow(alleles))
 idx<-which(alleles$al.x!=alleles$al.y)
 x.alleles <- alleles[idx,]$al.x
 names(x.alleles)<-alleles[idx,]$pid
 y.alleles <-  alleles[idx,]$al.y
 names(y.alleles)<-names(x.alleles)
 align.class[idx] <- g.class(x.alleles,y.alleles)
 print(table(align.class))
 M[,g.class:=align.class]
 sw <- align.class %in% c("rev","revcomp")
 M[sw,orig.z:=orig.z * -1]
 M[sw,OR:=1/OR]
 library(cowplot)
 ggplot(M,aes(x=orig.z,y=z_imp,color=r2.pred)) + geom_point()
 n1 <- 972
 n0 <- 1977
 ## next compute the OR derived from the imputed score
M[,imp_beta_linear:=z_imp/sqrt(2 * N.imp * maf * (1-maf))]
M[,imp_beta_cc:=z_imp/sqrt(2 * ((n1 * n0)/(n1+n0)) * r2.pred * maf * (1-maf))]
aa <- ggplot(M,aes(x=log(OR),y=imp_beta_linear,color=r2.pred)) + geom_point() + geom_abline()
bb <- ggplot(M,aes(x=log(OR),y=imp_beta_cc,color=r2.pred)) + geom_point() + geom_abline()
plot_grid(aa,bb)
## what is the difference in the estimated stderr between imputed and nonimputed
M[,imp.se:=1/sqrt(2 * ((n1 * n0)/(n1+n0)) * r2.pred * maf * (1-maf))]
M[,act.se:=log(OR)/orig.z]
M[,act.se.maf:=(log(OR)/orig.z) * sqrt(2 * ((n1 * n0)/(n1 + n0)) * r2.pred)]
M[,imp.se.maf:=1/sqrt(maf * (1-maf))]
ggplot(M,aes(x=act.se,y=imp.se,color=r2.pred)) + geom_point() + geom_abline()

## collate all results into single files
OUT_DIR <- '/home/ob219/share/Data/GWAS-summary/MyastheniaGravis_Renton_JAMA_Neurol_2015/ssimp_imputed'

early <- lapply(list.files(path='/home/ob219/rds/hpc-work/ssimp/mg',pattern='*renton_mg_early*',full.names=TRUE),function(f){
  fread(f)
}) %>% rbindlist
early <- early[order(chr,pos),]
write.table(early,file=file.path(OUT_DIR,'MyastheniaGravis_YoungOnset_Renton_JAMA_Neurol_2015_ssimp_imputed.tab'),sep="\t",quote=FALSE,row.names=FALSE)
library(parallel)
late <- mclapply(list.files(path='/home/ob219/rds/hpc-work/ssimp/mg',pattern='*renton_mg_late*',full.names=TRUE),function(f){
  fread(f)
},mc.cores=8) %>% rbindlist
late <- late[order(chr,pos),]
write.table(late,file=file.path(OUT_DIR,'MyastheniaGravis_LateOnset_Renton_JAMA_Neurol_2015_ssimp_imputed.tab'),sep="\t",quote=FALSE,row.names=FALSE)
combined <- mclapply(list.files(path='/home/ob219/rds/hpc-work/ssimp/mg',pattern='*renton_mg_with_or*',full.names=TRUE),function(f){
  fread(f)
},mc.cores=8) %>% rbindlist
combined <- combined[order(chr,pos),]
write.table(combined,file=file.path(OUT_DIR,'MyastheniaGravis_Renton_JAMA_Neurol_2015_ssimp_imputed.tab'),sep="\t",quote=FALSE,row.names=FALSE)

}
