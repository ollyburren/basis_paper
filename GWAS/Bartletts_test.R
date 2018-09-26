library(data.table)
library(magrittr)
## load in projections

ind.proj <- readRDS("/home/ob219/rds/hpc-work/as_basis/support/ind_proj_june10k.RDS")

## we want to test to see if PID which is quite a heterogenous phenotype exhibts a different variance
## in the projections for any components compared to something that we know is phenotypically homogenenous

## step 1 use Bartlett's test to test if variances between PID and JIA_sys are different
## for any component


## pid vs jiasys

test.DT <- ind.proj[trait %in% c('jiasys','pid'),]
M <- melt(test.DT,id.vars=c('trait','individual'),measure.vars=paste0('PC',1:11))
lapply(paste0('PC',1:11),function(pc){
  DT <- M[variable==pc,]
  data.table(pc=pc,p.Bartlett=bartlett.test(value~trait,DT)$p.value,traits=paste(unique(DT$trait),sep=',',collapse=','))
}) %>% rbindlist

## next compare with all of JIA but we need to mean centre for each group


jia.cat <- ind.proj[grepl("^jia",trait),]$trait %>% unique
jia.DT <- ind.proj[trait %in% jia.cat,]
M.jia <- melt(jia.DT,id.vars=c('trait','individual'),measure.vars=paste0('PC',1:11))
Mt.jia <- M.jia[,list(individual=individual,variable=variable,value=value-mean(value)),by='trait'][,trait:='jia']

Mf <- rbind(Mt.jia,M[trait=='pid'])
lapply(paste0('PC',1:11),function(pc){
  DT <- Mf[variable==pc,]
  data.table(pc=pc,p.Bartlett=bartlett.test(value~trait,DT)$p.value,traits=paste(unique(DT$trait),sep=',',collapse=','))
}) %>% rbindlist

## jdm vs jiasys

test.DT <- ind.proj[trait %in% c('jiasys','jdm'),]
M <- melt(test.DT,id.vars=c('trait','individual'),measure.vars=paste0('PC',1:11))
lapply(paste0('PC',1:11),function(pc){
  DT <- M[variable==pc,]
  data.table(pc=pc,p.Bartlett=bartlett.test(value~trait,DT)$p.value,traits=paste(unique(DT$trait),sep=',',collapse=','))
}) %>% rbindlist

## next compare with all of JIA but we need to mean centre for each group

jia.cat <- ind.proj[grepl("^jia",trait),]$trait %>% unique
jia.DT <- ind.proj[trait %in% jia.cat,]
M.jia <- melt(jia.DT,id.vars=c('trait','individual'),measure.vars=paste0('PC',1:11))
Mt.jia <- M.jia[,list(individual=individual,variable=variable,value=value-mean(value)),by='trait'][,trait:='jia']

Mf <- rbind(Mt.jia,M[trait=='jdm'])
lapply(paste0('PC',1:11),function(pc){
  DT <- Mf[variable==pc,]
  data.table(pc=pc,p.Bartlett=bartlett.test(value~trait,DT)$p.value,traits=paste(unique(DT$trait),sep=',',collapse=','))
}) %>% rbindlist

## pid vs jdm

test.DT <- ind.proj[trait %in% c('pid','jdm'),]
M <- melt(test.DT,id.vars=c('trait','individual'),measure.vars=paste0('PC',1:11))
lapply(paste0('PC',1:11),function(pc){
  DT <- M[variable==pc,]
  data.table(pc=pc,p.Bartlett=bartlett.test(value~trait,DT)$p.value,traits=paste(unique(DT$trait),sep=',',collapse=','))
}) %>% rbindlist
