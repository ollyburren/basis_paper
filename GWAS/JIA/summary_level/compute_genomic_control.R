## check the JIA GWAS for inflation etc.

DATA_DIR <- '/home/ob219/share/as_basis/GWAS/sum_stats'
jf <- list.files(path=DATA_DIR,pattern='jia.*',full.names=TRUE)

dt <- lapply(jf,fread)

names(dt) <- gsub("_unpub.tab","",basename(jf))

## compute lambda

lambdas <- lapply(dt,function(X){
  Z <- qnorm(X$p.value/2,lower.tail=FALSE)
  ## lambda is just this divided by the average chi squared
  round(median(Z^2,na.rm = TRUE)/qchisq(0.5,1),3)
})

% > lambdas
% $jiacc
% [1] 1.059

% $jiaeo
% [1] 1.016

% $jiaera
% [1] 0.988

% $jiapo
% [1] 1.009

% $jiapsa
% [1] 1.004

% $jiarfneg
% [1] 1.019

% $jiarfpos
% [1] 1.008

% $jiasys
% [1] 1.015
