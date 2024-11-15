##sobel -> sample partiioning -> tstatistic p value -> cauchy combination test
rm(list = ls())


library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

sim = function(alpha, beta, gamma, coef_X_alpha, coef_X_beta, 
               n, K, omega = 0.05)
{
  other_covariates = cbind(rep(1,n), runif(n), runif(n))
  ##Generating Data
  Y<- G <- S <- c()
  
  S = runif(n, 1,2)
  G = S*alpha + rowSums(coef_X_alpha*other_covariates) + rnorm(n)
  Y = G*beta + S*gamma + rowSums(coef_X_beta*other_covariates) + rnorm(n)
  
  
  
  ##Subsampling
  
  #K = 10
  niter = n*5
  
  single_iteration <- function(iter) {
    # Sample partitioning
    all <- sample.int(n = n, size = n, replace = FALSE)
    index <- split(all, rep(1:K, length.out = n))
    
    S_K <- numeric(K)
    for (i in 1:K) {
      G_s <- G[index[[i]]]
      S_s <- S[index[[i]]]
      Y_s <- Y[index[[i]]]
      other_covariates_s <- other_covariates[index[[i]], ]
      
      obj1 <- lm(G_s ~ -1 + cbind(S_s, other_covariates_s))
      obj2 <- lm(Y_s ~ -1 + cbind(G_s, S_s, other_covariates_s))
      t1 <- summary(obj1)[["coefficients"]][1, 3]
      t2 <- summary(obj2)[["coefficients"]][1, 3]
      
      S_K[i] <- (t1 * t2) / sqrt(t1^2 + t2^2)
    }
    
    # Compute t-statistic and p-value
    tstat_iter <- sqrt(K) * mean(S_K) / sqrt(var(S_K))
    pval_iter <- 2 * pt(-abs(tstat_iter), df = K - 1)
    
    # Return results as a list
    list(tstat = tstat_iter, pval = pval_iter)
  }
  
  # Parallelize across 1:niter
  results <- mclapply(1:niter, single_iteration, mc.cores = detectCores()-1)
  
  # Extract tstat and pval from results
  tstat <- sapply(results, `[[`, "tstat")
  pval <- sapply(results, `[[`, "pval")
  
  #p_global = combination.test(pval, method ="Levy" ,truncate.threshold = 0.97)
  wts = runif(niter)
  wts = wts/sum(wts)
  #hist(wts)
  CCT = sum(tan(3.141593*(0.5 - pval))*wts) ##Cauchy combination statistic
  p_cct = 1- pcauchy(CCT)
  #ABtest
  
  mat = DBmypackage6two::One_MC_DB_X_noProj_C(cbind(S,G,Y), other_covariates,
                                              n, omega, num_med = 1, lambda1all = 2, lambda2all = 2, N_boot_out = 500, N_boot_in = 0)
  pvalue = mat[1,1]
  
  #fitting full model for remaining methods
  object1 = lm(G ~ -1 + cbind(S, other_covariates))
  object2 = lm(Y ~ -1 + cbind(G, S, other_covariates))
  t1 = summary(object1)[["coefficients"]][1,3]
  t2 = summary(object2)[["coefficients"]][1,3]
  p1 = summary(object1)[["coefficients"]][1,4]
  p2 = summary(object2)[["coefficients"]][1,4]
  ##fix this part. wtf is happening
  
  
  ##sobel
  T.sobel = t2/sqrt(1+(t2/t1)^2)
  p.sobel = 2*pnorm(-abs(T.sobel),mean = 0, sd=1, lower.tail = TRUE)
  
  return(c(p_cct, pvalue, max(p1,p2), p.sobel))
  
}


