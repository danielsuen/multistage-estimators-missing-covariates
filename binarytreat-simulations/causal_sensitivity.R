library(simsurv)
library(flexsurv)
library(MASS)
library(nnet)
library(bindata)
library(doSNOW)
library(doParallel)

source("functions_causalsim.R")

################## simulate


totalIter = 1000

ipw_R_matrix = array(NA, c(totalIter, 3))
reg_R_matrix = array(NA, c(totalIter, 3))


n_iter = 1

# parallel
numCores = detectCores()
cluster = makeCluster(numCores - 2)
registerDoSNOW(cluster)

# progress bar
pb = txtProgressBar(max = totalIter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


xi = seq(0,4,length.out=101)
rho = seq(-2,1,length.out=101)

n = 1000

for (kk in 1:length(xi))
{
  writeLines(paste("\nn: ", n, ", xi: ", xi[kk], " , rho: ", rho[kk], "\n", sep=""))
  output = foreach (ii = 1:totalIter, .combine='comb', .packages=c('MASS','flexsurv','simsurv','nnet','bindata','mvtnorm','matrixStats'), .multicombine=TRUE,
                  .init=list(list(), list(), list(),
                             list(), list(), list()), .options.snow = opts) %dopar%
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

      ipw_R = estimateIPW(Xobs, R, Y, A, rho[kk])
      
      ipw_R_ipw_A = ipw_R$ipw_R_ipw_A
      ipw_R_reg_A = ipw_R$ipw_R_reg_A
      ipw_R_dr_A = ipw_R$ipw_R_dr_A
      
      ########################## multiple imputation
      
      
      #### perform stacked imputation
      
      M = 50

      reg_R = estimateReg(Xobs, R, Y, A, M, xi[kk])
      
      reg_R_ipw_A = reg_R$reg_R_ipw_A
      reg_R_reg_A = reg_R$reg_R_reg_A
      reg_R_dr_A = reg_R$reg_R_dr_A
      
      list(ipw_R_ipw_A, ipw_R_reg_A, ipw_R_dr_A, reg_R_ipw_A, reg_R_reg_A, reg_R_dr_A)
    }
  
  ipw_R_matrix[, 1] = unlist(output[[1]])
  ipw_R_matrix[, 2] = unlist(output[[2]])
  ipw_R_matrix[, 3] = unlist(output[[3]])
  
  reg_R_matrix[, 1] = unlist(output[[4]])
  reg_R_matrix[, 2] = unlist(output[[5]])
  reg_R_matrix[, 3] = unlist(output[[6]])

  save(reg_R_matrix, file=paste("causal_sim_MI_sensitivity_",xi[kk],".RData", sep=""))
  save(ipw_R_matrix, file=paste("causal_sim_ipw_sensitivity_",rho[kk],".RData", sep=""))

}
