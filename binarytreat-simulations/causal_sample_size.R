library(simsurv)
library(flexsurv)
library(MASS)
library(nnet)
library(bindata)
library(doSNOW)
library(doParallel)

source("functions_causalsim.R")

################## simulate

M = 50
totalIter = 1000

ipw_R_matrix = array(NA, c(4, totalIter, 3))
reg_R_matrix = array(NA, c(4, totalIter, 3))

sample_sizes = c(1000, 2000, 5000, 10000)

n_iter = 1


# parallel
numCores = detectCores()
cluster = makeCluster(numCores - 2)
registerDoSNOW(cluster)

# progress bar
pb = txtProgressBar(max = totalIter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

for (n in sample_sizes)
{
  writeLines(paste("\nn: ", n, "\n", sep=""))
  output = foreach (ii = 1:totalIter, .combine='comb', .packages=c('MASS','flexsurv','simsurv','nnet','bindata','mvtnorm','matrixStats'), .multicombine=TRUE,
                  .init=list(list(), list(), list(),
                             list(), list(), list(),
                             list()), .options.snow = opts) %dopar%
    {
      
      ########################## generate data
      
      data = gen_complete_data(n)
      
      Xobs = data$Xobs
      X = data$X
      R = data$R
      Y = data$Y
      A = data$A
      
      params = estimate_params(Xobs, R, Y, A)
      empirical_probs = estimate_empirical_probs(R, Y, A)
      
      ######################## ipw
      zeta = 0
      
      ipw_R = estimateIPW(Xobs, R, Y, A, zeta)
      
      ipw_R_ipw_A = ipw_R$ipw_R_ipw_A
      ipw_R_reg_A = ipw_R$ipw_R_reg_A
      ipw_R_dr_A = ipw_R$ipw_R_dr_A
      
      inv_prop_score_R = ipw_R$inv_prop_score_R
      
      ########################## multiple imputation
      
      
      #### perform stacked imputation
      
      M = 50
      kappa = 0
      
      reg_R = estimateReg(Xobs, R, Y, A, M, kappa)
      
      reg_R_ipw_A = reg_R$reg_R_ipw_A
      reg_R_reg_A = reg_R$reg_R_reg_A
      reg_R_dr_A = reg_R$reg_R_dr_A
      
      list(ipw_R_ipw_A, ipw_R_reg_A, ipw_R_dr_A, reg_R_ipw_A, reg_R_reg_A, reg_R_dr_A, max(inv_prop_score_R))
    }
  
  ipw_R_matrix[n_iter, , 1] = unlist(output[[1]])
  ipw_R_matrix[n_iter, , 2] = unlist(output[[2]])
  ipw_R_matrix[n_iter, , 3] = unlist(output[[3]])
  
  reg_R_matrix[n_iter, , 1] = unlist(output[[4]])
  reg_R_matrix[n_iter, , 2] = unlist(output[[5]])
  reg_R_matrix[n_iter, , 3] = unlist(output[[6]])


  n_iter = n_iter + 1
}

save(ipw_R_matrix, reg_R_matrix, sample_sizes, file=paste("causal_samplesize_simulation.RData", sep=""))
