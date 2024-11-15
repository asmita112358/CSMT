##Subbel power simulations

source("~/Downloads/CSMT/R_functions/main_fun.R")
# Load the parallel library for parallel computing
library(parallel)

# Initialize variables
n <- 300
K <- 8
gamma <- 0.1
omega <- 0.05
coef_X_alpha <- c(1, 1, 1)
coef_X_beta <- c(1, 1, 1)
#div = c(0.02, 0.1, 0.5, 1, 2, 3, 4)
#mul = 0.1
str = seq(0.1, 0.5, length.out = 9)
alpha = str
beta = str


# Prepare matrices to store results
power <- matrix(nrow = length(alpha), ncol = 4)
colnames(power) <-  c("Subbel", "ABtest", "MaxP", "Sobel")

N_rep <- 500  # Number of repetitions
power = matrix(nrow = length(alpha), ncol = 4)

# Define a function to perform a single simulation run
for(i in 1:length(alpha)){
  # Create the result matrix for the repetitions
  mc <- replicate(N_rep, {
    alphaj <-alpha[i]
    betaj <- beta[i]
    
    
    # Run the sim function and return the result
    sim(alphaj, betaj, gamma, coef_X_alpha, coef_X_beta, n, K)
  })
  
  power[i,] = rowMeans(mc <= omega)
}


# Prepare the output data frame and save it
#intensity_factor <- factor(c("low", "medium", "high"))
data <- data.frame(str , power)

print(data)
write.csv(data, "power300.csv")
beep(2)
