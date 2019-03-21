### use pooled variance to do t-test taking into account sharing between control for studies.

RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/19_12_18_summary_results.RDS'
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


## test where one or more is different from control
res.DT[trait %in% c('bb_ankylosing.spondylitis','jia_ERA') & p.adj<0.05,]$variable
era_ankspond <- compute_t(tg=list(g1='jia_ERA',g2='bb_ankylosing.spondylitis'),pc='PC1')

pcs <- res.DT[trait %in% c('bb_psoriatic.arthropathy','bb_psoriasis','jia_PsA') & p.adj<0.05,]$variable
psa_ukbbpsa <- compute_t(tg=list(g1='jia_PsA',g2='bb_psoriatic.arthropathy'),pc='PC3')
#psa_ukbbpso <- compute_t(tg=list(g1='jia_PsA',g2='bb_psoriasis'),pc='PC3',psa_ukbbpsa$covM)
dat <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='jia_PsA',g2='bb_psoriasis'),pc=pc,psa_ukbbpsa$covM)$t
}) %>% rbindlist
dat2 <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='jia_PsA',g2='bb_psoriatic.arthropathy'),pc=pc,psa_ukbbpsa$covM)$t
}) %>% rbindlist

dat3 <- lapply(pcs,function(pc){
  compute_t(tg=list(g1='bb_psoriasis',g2='bb_psoriatic.arthropathy'),pc=pc,psa_ukbbpsa$covM)$t
}) %>% rbindlist

dat <- rbindlist(list(dat,dat2,dat3))
dat[,p.adj:=p.adjust(p.value,method="bonferroni")]



compute_t(tg=list(g1='jia_PsA',g2='bb_psoriasis'),pc='PC1',psa_ukbbpsa$covM)$t

myo.res <- compute_t(tg=list(g1='pm_myogen',g2=c('jdm_myogen','dm_myogen')),pc='PC10')
egpa.res <- compute_t(tg=list(g1='anca_Neg',g2=c('mpo_Pos')),pc='PC6')
jia.res <- compute_t(tg=list(g1=c('jia_sys','jia_ERA'),g2=c('jia_EO','jia_PO','jia_PsA','jia_RFneg','jia_RFpos')),pc='PC3')
jia.res.noshare <- compute_t_no_share(tg=list(g1=c('jia_sys','jia_ERA'),g2=c('jia_EO','jia_PO','jia_PsA','jia_RFneg','jia_RFpos')),pc='PC3')


## for pc10 can we compute the shared variance

library(xtable)
rbindlist(list(myo.res$t,egpa.res$t,jia.res$t)) %>% xtable
