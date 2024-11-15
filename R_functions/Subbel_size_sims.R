##Size sims for Subbel
source("~/Downloads/CSMT/R_functions/main_fun.R")
n = 550
K = 11
#Coefficients
alpha = 0.1
beta = 0.1
gamma = 0.1
omega = 0.05
coef_X_alpha = c(1,1,1)
coef_X_beta = c(1,1,1)
##Checking size, qqplot
prop = c(0.4, 0.3, 0.3)
N_rep = 500       ##Number of monte carlo repetitions
mc = matrix(nrow = N_rep, ncol = 4)
samp = sample(1:3, N_rep, replace = T, prob = prop)
for(j in 1: N_rep)
{
  alphaj = 0 + alpha*(samp[j]== 2)
  betaj = 0 + beta*(samp[j] == 3)
  #print(alphaj*betaj)
  mc[j,] = sim(alphaj, betaj,  gamma = 0.1, coef_X_alpha, coef_X_beta, 
               n, K = K, omega = 0.05)
}

size = apply((mc<=omega), 2, mean)
se_size = sqrt(size*(1 - size)/ N_rep)
write.csv(mc, "dense550.csv")
