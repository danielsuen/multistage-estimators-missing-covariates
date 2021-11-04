library(simsurv)
library(flexsurv)
library(MASS)
library(nnet)
library(bindata)
library(doSNOW)
library(doParallel)

source("functions_coxsim.R")

###### parallelization

# parallel
numCores = detectCores()
cluster = makeCluster(numCores - 2)
registerDoSNOW(cluster)

# progress bar
pb = txtProgressBar(max = totalIter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

################## simulate

M = 50
totalIter = 1000

transformed_estimates = array(NA, c(4, totalIter, 2))
MI_estimates = array(NA, c(4, totalIter, 2))
ipw_estimates = array(NA, c(4, totalIter, 2))

sample_sizes = c(350, 500, 1000, 2000)

n_iter = 1

## sensitivity parameters - no perturbation
zeta = 1/2
rho = 0

for (n in sample_sizes)
{
  writeLines(paste("\nn: ", n, "\n", sep=""))
  output = foreach (ii = 1:totalIter, .combine='comb', .packages=c('MASS','flexsurv','simsurv','nnet','bindata'), .multicombine=TRUE,
                  .init=list(list(), list(), list(), list(), list(), list()), .options.snow = opts) %dopar%
  {
  
    ########################## generate data
      
    dataset = gen_complete_data(n)
    X = dataset$X
    R = dataset$R
    delta = dataset$delta
    Y = dataset$Y
    X_obs = dataset$X_obs
    Y_cc = dataset$Y_cc
    delta_cc = dataset$delta_cc
    X_cc = dataset$X_cc
      
    ########################## EM
      
    
    # estimate logistic
    logistic_coef = estimate_logistic(X_obs, Y, delta, R)
    eta_10_est = logistic_coef$eta_10_est
    eta_01_est = logistic_coef$eta_01_est
      
    # estimate empirical p(x1, x2, R=11)
    cc_empirical_prob = get_empirical_probs(X_obs)
      
    # EM
    EM_output = EM(X_obs, Y, delta, R, eta_10_est, eta_01_est, cc_empirical_prob)
      
    # transformed MLE
    transformed = rep(NA,2)
    transformed[1] = log(EM_output$gamma_est[2]/EM_output$gamma_est[1]) / 2 + (log(EM_output$gamma_est[4]/EM_output$gamma_est[1]) - log(EM_output$gamma_est[3]/EM_output$gamma_est[1])) / 2
    transformed[2] = log(EM_output$gamma_est[3]/EM_output$gamma_est[1]) / 2 + (log(EM_output$gamma_est[4]/EM_output$gamma_est[1]) - log(EM_output$gamma_est[2]/EM_output$gamma_est[1])) / 2
      
    ########################## multiple imputation
      
    M = 50
    stacked_impute = doImputation(X_obs, Y, delta, R, cc_empirical_prob, EM_output, zeta, M)
      
    impute_output = coxph(Surv(Y, delta) ~ X1+X2, data = stacked_impute)$coefficients
      
    ############################ ipw
      
    lambda = get_IPW_weights(X_cc, Y_cc, delta_cc, eta_10_est, eta_01_est, rho)
      
    ipw_output = coxph(Surv(Y_cc, delta_cc) ~ X_cc[,1] + X_cc[,2], weights=lambda)$coefficients
      
    list(transformed[1], transformed[2],
          impute_output[1], impute_output[2],
          ipw_output[1], ipw_output[2])
  }
    
  transformed_estimates[n_iter, , 1] = unlist(output[[1]])
  transformed_estimates[n_iter, , 2] = unlist(output[[2]])
  MI_estimates[n_iter, , 1] = unlist(output[[3]])
  MI_estimates[n_iter, , 2] = unlist(output[[4]])
  ipw_estimates[n_iter, , 1] = unlist(output[[5]])
  ipw_estimates[n_iter, , 2] = unlist(output[[6]])

  n_iter = n_iter + 1
}
save(transformed_estimates, MI_estimates, ipw_estimates, sample_sizes, file=paste("survival_samplesize_simulation.RData", sep=""))
