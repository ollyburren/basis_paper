### use pooled variance to do t-test taking into account sharing between control for studies.

#RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/03_07_19_0619_summary_results.RDS'
RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/16_07_19_0619_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)

## pooled variance

## we want to test our groups of things to see if there are differences in particular PC's
## we can do this with a t.test but we need to estimate the pooled mean (easy) and the
## pooled variance (harder)

## For UKBB and JIA the summary statistcs have shared controls (and in bb shared cases).
## Let's consider the UKBB
## we know for each trait the number of cases and controls but between traits we don't
## know what the sample overlap is. If we assume that for two traits,i and j they share
## the min(n_i,n_j) then our calculations will be conservative.

## from li and sullivan we can compute
## here C is the sharing freq between cases for a disease

compute_cor_bb <- function(n_i1,n_j1,n_i0,n_j0,C=100){
  #shared controls
  n_ij0 <- min(n_i0,n_j0)
  #shared cases
  n_ij1 <- min(n_i1,n_j1)/C
  n_i <- n_i1 + n_i0
  n_j <- n_j1 + n_j0
  t1 <- n_ij0 * sqrt( (n_i1 * n_j1)/(n_i0 * n_j0 ) )
  t2 <- n_ij1 * sqrt( (n_i0 * n_j0)/(n_i1 * n_j1 ) )
  #(t1 + t2)/sqrt(n_i * n_j)
  ## note that for BB assume a fixed sample size therefore
  (t1 + t2)/n_i # or (t1 + t2)/n_j as n_i==n_j
}

## jia traits will be different as there are no shared cases only shared n_controls
## so don't need t2

compute_cor_jia <- function(n_i1,n_j1,n_i0,n_j0){
  #shared controls
  n_ij0 <- min(n_i0,n_j0)
  n_i <- n_i1 + n_i0
  n_j <- n_j1 + n_j0
  t1 <- n_ij0 * sqrt( (n_i1 * n_j1)/(n_i0 * n_j0 ) )
  t1/sqrt(n_i * n_j)
}

## under  assumptions that shared controls are similar between studies
## which hold for bb and JIA we can cancel terms
## this will overestimate the correlation so give us a
## conservative estimate

compute_conservative_cor_bb <- function(n_i1,n_j1,n_i0,n_j0){
  n_i <- n_i1 + n_i0
  n_j <- n_j1 + n_j0
  (sqrt(n_i1 * n_j1) + sqrt(n_i0 * n_j0))/sqrt(n_i * n_j)
  ## this is what chris wrote on slack but I think wrong
  #(sqrt(n_i1 + n_j1) + sqrt(n_i0 + n_j0))/sqrt(n_i * n_j)
  ## note that for bb where n is fixed across studies
  #(sqrt(n_i1 * n_j1) + sqrt(n_i0 * n_j0))/n_i ## or n_j
}

compute_cov_noshare <- function(DT){
  M <- matrix(0,nrow=nrow(DT),ncol=nrow(DT),dimnames=list(DT$trait,DT$trait))
  M[upper.tri(M,diag=FALSE)] <- t(M)[upper.tri(M,diag=FALSE)]
  ## note that the diag is 1
  diag(M) <- 1
  #pheatmap(M)
  ## 3 compute covariance matrix by multipying elements by sqrt(v_i * v_j)
  ## note for actual application this will depend on principal axis that we
  ## are using
  (tcrossprod(DT$variance,DT$variance) %>% sqrt) * M
}

compute_cov <- function(DT){
  M <- matrix(0,nrow=nrow(DT),ncol=nrow(DT),dimnames=list(DT$trait,DT$trait))
  for(i in 1:nrow(DT)){
    for(j in 1:nrow(DT)){
      if(i==j){
        break()
      }
      smix <- DT$bb[i] + DT$bb[j]
      if(smix == 2){
        ## both biobank
        #cor <- compute_conservative_cor_bb(
        cor <- compute_cor_bb(
          n_i1=DT$n1[i],
          n_j1=DT$n1[j],
          n_i0=DT$n0[i],
          n_j0=DT$n0[j])
        }else if(smix == 1){
          ## jia and bb assume no correlation
          cor <- 0
        }else{
          ## both jia - NB this is not correct if we test between groups where
          ## the groups have shared controls but there are no shared samples between
          ## groups in that case cor should be 0
          cor <- compute_cor_jia(
            n_i1=DT$n1[i],
            n_j1=DT$n1[j],
            n_i0=DT$n0[i],
            n_j0=DT$n0[j])
          }
          M[i,j] <- cor
        }
      }
      ## create a symmetrical matrix
      M[upper.tri(M,diag=FALSE)] <- t(M)[upper.tri(M,diag=FALSE)]
      ## note that the diag is 1
      diag(M) <- 1
      #pheatmap(M)
      ## 3 compute covariance matrix by multipying elements by sqrt(v_i * v_j)
      ## note for actual application this will depend on principal axis that we
      ## are using
      (tcrossprod(DT$variance,DT$variance) %>% sqrt) * M
}


## inefficient as each time we compute covM for a PC whereas we only need to do this once

compute_t <- function(tg,pc,covM){
  cor.DT <-res.DT[variable==pc,.(trait,value,n1,n0=n-n1,variance)]
  cor.DT <- cor.DT[trait!='unclassifiable',]
  cor.DT[,bb:=grepl("^bb",trait)]
  if(missing(covM)){
    message("CovM missing computing")
    covM <- compute_cov(cor.DT)
  }else{
    message("Using supplied CovM")
  }
  ## with this covariance matrix we can now estimate the pooled variance for any arb group of diseases
  ## code to compute groups - we use complete linkage to start with
  ind <- lapply(tg,function(tl){
    idx <- which(cor.DT$trait %in% tl)
  })
  ## using the covariance matrix between all traits computed above compute the pooled
  ## variance across both groups
  gvar <- lapply(ind,function(i){
    var <- sum(covM[i,i])
  })
  ## compute the pooled mean for both samples
  pmean <- lapply(ind,function(i){
    cor.DT[i,list(mean(value))]$V1
  })
  n <- lapply(ind,function(i){
    length(i)
  })
  ## compute the difference
  diffMean <- pmean[[1]] - pmean[[2]]
  ## t-test assesses the difference in means but these are still correlated as group 1 and group 2 still covary
  ## (mean(X_1),V_1), (mean(X_2),V_2)
  ## d = X_1 - X_2
  ## v(d) = V_1 + V_2 + 2cov(mean(X_1),mean(X_2))
  ## we have to adjust the pooled variance to take this into account
  V_1 <- gvar[[1]] * (1/n[[1]]^2)
  V_2 <- gvar[[2]] * (1/n[[2]]^2) ## for this special case (where there are two groups) V_1 + V_2 = sum(covM) !
  btw_sample_cov <- covM[ind[[1]],ind[[2]]] %>% sum
  ## const is 2/n_1^2,n_2^2 - here n is the number of traits in each group
  const <- 2/(sapply(ind,function(x) length(x)^2) %>% prod)
  var_d <- V_1 + V_2 - (const * btw_sample_cov)
  n <- sapply(ind,length)
  #t.stat <- diffMean/(sqrt(var_d) * sqrt(1/n[[1]] + 1/n[[2]]))
  t.stat <- diffMean/(sqrt(var_d))
  P <- pnorm(abs(t.stat),lower.tail=FALSE) * 2
  #data.table(p.value=P,t.stat=t.stat,df=df,pc=pc,diff.mean=diffMean,var1=V_1,var2=V_2,overall_var=var_d)
  t.DT <- data.table(pc=pc,p.value=P,t.stat=t.stat,set1=paste(tg[[1]],collapse=','),set2=paste(tg[[2]],collapse=','))
  return(list(t=t.DT,covM=covM))
}


compute_t_no_share <- function(tg,pc,covM){
  cor.DT <-res.DT[variable==pc,.(trait,value,n1,n0=n-n1,variance)]
  cor.DT <- cor.DT[trait!='unclassifiable',]
  cor.DT[,bb:=grepl("^bb",trait)]
  covM <- compute_cov_noshare(cor.DT)
  ## with this covariance matrix we can now estimate the pooled variance for any arb group of diseases
  ## code to compute groups - we use complete linkage to start with
  ind <- lapply(tg,function(tl){
    idx <- which(cor.DT$trait %in% tl)
  })
  ## using the covariance matrix between all traits computed above compute the pooled
  ## variance across both groups
  gvar <- lapply(ind,function(i){
    var <- sum(covM[i,i])
  })
  ## compute the pooled mean for both samples
  pmean <- lapply(ind,function(i){
    cor.DT[i,list(mean(value))]$V1
  })
  n <- lapply(ind,function(i){
    length(i)
  })
  ## compute the difference
  diffMean <- pmean[[1]] - pmean[[2]]
  ## t-test assesses the difference in means but these are still correlated as group 1 and group 2 still covary
  ## (mean(X_1),V_1), (mean(X_2),V_2)
  ## d = X_1 - X_2
  ## v(d) = V_1 + V_2 + 2cov(mean(X_1),mean(X_2))
  ## we have to adjust the pooled variance to take this into account
  V_1 <- gvar[[1]] * (1/n[[1]]^2)
  V_2 <- gvar[[2]] * (1/n[[2]]^2) ## for this special case (where there are two groups) V_1 + V_2 = sum(covM) !
  btw_sample_cov <- covM[ind[[1]],ind[[2]]] %>% sum
  ## const is 2/n_1^2,n_2^2 - here n is the number of traits in each group
  const <- 2/(sapply(ind,function(x) length(x)^2) %>% prod)
  var_d <- V_1 + V_2 - (const * btw_sample_cov)
  n <- sapply(ind,length)
  #t.stat <- diffMean/(sqrt(var_d) * sqrt(1/n[[1]] + 1/n[[2]]))
  t.stat <- diffMean/(sqrt(var_d))
  P <- pnorm(abs(t.stat),lower.tail=FALSE) * 2
  #data.table(p.value=P,t.stat=t.stat,df=df,pc=pc,diff.mean=diffMean,var1=V_1,var2=V_2,overall_var=var_d)
  t.DT <- data.table(pc=pc,p.value=P,t.stat=t.stat,set1=paste(tg[[1]],collapse=','),set2=paste(tg[[2]],collapse=','))
  return(list(t=t.DT,covM=covM))
}


## compute the covariance matrix so that we can do use and save a lot of effort
scovM <- compute_t(tg=list(g1='pr3_meta',g2='mpo_meta'),pc='PC1')$covM
pcs <- paste('PC',1:11,sep='')


dat.ladavst1d <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='bb_SRD:type.1.diabetes',g2='cousminer_lada'),pc=pc,scovM)$t
}) %>% rbindlist

dat.ladavst1d[,bonf.p:=p.adjust(p.value,method="bonferroni")]

dat.egpa <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='mpo_Pos',g2='anca_Neg'),pc=pc,scovM)$t
}) %>% rbindlist

dat.aav <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='pr3_meta',g2='mpo_meta'),pc=pc,scovM)$t
}) %>% rbindlist

## for labtalk select differences between subtypes that are bonferroni significant

lt.results <- rbindlist(list(dat.egpa,dat.aav))[,bonf.p:=p.adjust(p.value,method="bonferroni")][bonf.p<0.05,]

## let us plot so we have something to show

lt.results <- rbindlist(list(dat.egpa,dat.aav))

lt.results[set1=='mpo_Pos',disease:='EGPA:mpo+ vs anca-']
lt.results[set1=='pr3_meta',disease:='AAV: mpo+ vs pr3+']
lt.results[,pc:=factor(pc,levels=paste0('PC',1:11))]
sig.level <- qnorm(0.05/44,lower.tail=FALSE)
library(ggplot2)
library(cowplot)
pp2 <- ggplot(lt.results,aes(x=pc,y=abs(t.stat),fill=disease)) + geom_bar(stat="identity",position=position_dodge()) +
geom_hline(yintercept=c(sig.level),color='red',lty=2) + xlab('PC') + ylab('t-statistic') +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

save_plot(file="~/tmp/smith_090719_vasc_barplot.pdf",pp2,base_width=11)

dat.egpampovsaavmpo <- lapply(pcs,function(pc){
  compute_t(tg=list(g1=c('mpo_Pos'),g2=c('mpo_meta')),pc=pc,scovM)$t
}) %>% rbindlist

dat.egpaancavsaavmpo <- lapply(pcs,function(pc){
  compute_t(tg=list(g1=c('anca_Neg'),g2=c('mpo_meta')),pc=pc,scovM)$t
}) %>% rbindlist


lta.results <- rbindlist(list(dat.egpampovsaavmpo,dat.egpaancavsaavmpo))[,bonf.p:=p.adjust(p.value,method="bonferroni")][bonf.p<0.05,]

## let us plot so we have something to show

lta.results <- rbindlist(list(dat.egpampovsaavmpo,dat.egpaancavsaavmpo))

lta.results[set1=='mpo_Pos',disease:='EGPA:mpo+ vs AAV:mpo+']
lta.results[set1=='anca_Neg',disease:='EGPA:ANCA- vs AAV:mpo+']
lta.results[,pc:=factor(pc,levels=paste0('PC',1:11))]
sig.level <- qnorm(0.05/44,lower.tail=FALSE)
library(ggplot2)
library(cowplot)
pp3 <- ggplot(lta.results,aes(x=pc,y=abs(t.stat),fill=disease)) + geom_bar(stat="identity",position=position_dodge()) +
geom_hline(yintercept=c(sig.level),color='red',lty=2) + xlab('PC') + ylab('t-statistic') +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

save_plot(file="~/tmp/smith_090719_egpavsaavmpo_barplot.pdf",pp3,base_width=11)


dat.egpavsaav <- lapply(pcs,function(pc){
  compute_t(tg=list(g1=c('mpo_Pos','anca_Neg'),g2=c('pr3_meta','mpo_meta')),pc=pc,scovM)$t
}) %>% rbindlist



dat.jia.sys <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='jia_sys_19',g2=c('jia_EO_19','jia_PO_19','jia_PsA_19','jia_RFneg_19','jia_RFpos_19','jia_undiff_19')),pc=pc,scovM)$t
}) %>% rbindlist

dat.jia.era <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='jia_ERA_19',g2=c('jia_EO_19','jia_PO_19','jia_PsA_19','jia_RFneg_19','jia_RFpos_19','jia_undiff_19')),pc=pc,scovM)$t
}) %>% rbindlist

dat.nmo <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='NMO_IgGNeg',g2='NMO_IgGPos'),pc=pc,scovM)$t
}) %>% rbindlist

dat.myo <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='pm_myogen',g2=c('dm_myogen','jdm_myogen')),pc=pc,scovM)$t
}) %>% rbindlist


res <- rbindlist(list(
egpa=dat.egpa,
egpavsaav = dat.egpavsaav,
aav = dat.aav,
sys=dat.jia.sys,
era=dat.jia.era,
nmo=dat.nmo,
myo = dat.myo
))

res[,p.adj:=p.adjust(p.value,method="bonferroni")]
saveRDS(res,"~/share/as_basis/GWAS/RESULTS/trait_comparison.RDS")
