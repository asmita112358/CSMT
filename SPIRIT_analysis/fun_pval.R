library(purrr)
library(tidyr)
library(dplyr)
##Function to input X, M and Y and output p-values for maxP, Sobel, Subbel and ABtest
# Function to divide data into K subsamples
divide_into_subsamples <- function(data, K) {
  n <- nrow(data)  # Total number of samples
  base_size <- floor(n / K)  # Size for the first K-1 subsamples
  leftover_size <- n - (base_size * (K - 1))  # Size for the last subsample
  subsamples <- vector("list", K)
  data <- data[sample(nrow(data)), ]
  for (i in 1:(K - 1)) {
    subsamples[[i]] <- data[(1 + (i - 1) * base_size):(i * base_size), ]
  }
  subsamples[[K]] <- data[((K - 1) * base_size + 1):nrow(data), ]
  return(subsamples)
}


fun.pval = function(X, M, Y, K = 4, niter = 2000, other_covariates, omega = 0.05)
{
  tstat = c()
  pval = c()
  M = scale(M)
  Y = scale(Y)
  data = na.omit(data.frame(X,M,Y, other_covariates))
  n = nrow(data)
  X = data$X
  M = data$M
  Y = data$Y
  other_covariates = data[,-(1:3)]
  
  for(iter in 1:niter)
  {
    #Sample partition
    splits <- split(data, data$X)
    
    sub1 = divide_into_subsamples(splits[[1]], K)
    sub2 = divide_into_subsamples(splits[[2]], K)
    sub = Map(rbind,sub1, sub2)
    
    #identical(sort(do.call('c', index)), 1:n) #check the sample partition
    S_K = c()
    for(i in 1:K)
    {
      data_s = sub[[i]]
      X_s = data_s$X
      M_s = data_s$M
      Y_s = data_s$Y
      other_covariates_s = as.matrix(data_s[,-(1:3)])
      obj1 = lm(M_s ~ -1 + cbind(X_s, other_covariates_s))
      obj2 = lm(Y_s ~ -1 + cbind(M_s, X_s, other_covariates_s))
      t1 = summary(obj1)[["coefficients"]][1,3]
      t2 = summary(obj2)[["coefficients"]][1,3]
      
      S_K[i] = (t1*t2)/sqrt(t1^2 + t2^2)
    } 
    tstat[iter] = sqrt(K)*mean(S_K)/sqrt(var(S_K))
    pval[iter] = 2*pt(-abs(tstat[iter]), df = K-1)
  }
  #wts = runif(niter)
  #wts = wts/sum(wts)
  #hist(wts)
  CCT = mean(tan(3.141593*(0.5 - pval))) ##Cauchy combination statistic
  p_cct = 1- pcauchy(CCT)
  
  #ABtest
  
  mat = DBmypackage6two::One_MC_DB_X_noProj_C(cbind(X,M,Y), as.matrix(other_covariates),
                                              n, omega, num_med = 1, lambda1all = 1, lambda2all = 1, N_boot_out = 500, N_boot_in = 0)
  pvalue = mat[1,1]
  
  #fitting full model for remaining methods
  object1 = lm(M ~ -1 + cbind(X, as.matrix(other_covariates)))
  object2 = lm(Y ~ -1 + cbind(M, X, as.matrix(other_covariates)))
  t1 = summary(object1)[["coefficients"]][1,3]
  t2 = summary(object2)[["coefficients"]][1,3]
  p1 = summary(object1)[["coefficients"]][1,4]
  p2 = summary(object2)[["coefficients"]][1,4]
  
  
  
  ##sobel
  T.sobel = t2/sqrt(1+(t2/t1)^2)
  p.sobel = 2*pnorm(-abs(T.sobel),mean = 0, sd=1, lower.tail = TRUE)
  
  return(c(p_cct, pvalue, max(p1,p2), p.sobel))
  
}
