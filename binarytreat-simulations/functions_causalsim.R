library(MASS)
library(nnet)
library(bindata)
library(matrixStats)
library(mvtnorm)
library(latex2exp)

############## output of for loop
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

####### generate data according to CCMV

gen_complete_data = function(n)
{
  index = sample(144,n,replace=TRUE)
  Y = rep(NA, n)
  A = rep(NA, n)
  R = rep(NA, n)
  X = matrix(NA, n, 2)
  
  probs = c(1/16, 1/48, 1/16, 5/48,
            1/12, 1/24, 1/24, 1/6,
            3/32, 1/32, 1/24, 1/12,
            1/16, 1/48, 1/36, 1/18)
  endpoints = cumsum(probs*144)
  
  # R=11, Y=0, A=0
  Y[index <= endpoints[1]] = 0
  A[index <= endpoints[1]] = 0
  R[index <= endpoints[1]] = 3
  
  # R=11, Y=0, A=1
  Y[index > endpoints[1] & index <= endpoints[2]] = 0
  A[index > endpoints[1] & index <= endpoints[2]] = 1
  R[index > endpoints[1] & index <= endpoints[2]] = 3
  
  # R=11, Y=1, A=0
  Y[index > endpoints[2] & index <= endpoints[3]] = 1
  A[index > endpoints[2] & index <= endpoints[3]] = 0
  R[index > endpoints[2] & index <= endpoints[3]] = 3
  
  # R=11, Y=1, A=1
  Y[index > endpoints[3] & index <= endpoints[4]] = 1
  A[index > endpoints[3] & index <= endpoints[4]] = 1
  R[index > endpoints[3] & index <= endpoints[4]] = 3
  
  # R=10, Y=0, A=0
  Y[index > endpoints[4] & index <= endpoints[5]] = 0
  A[index > endpoints[4] & index <= endpoints[5]] = 0
  R[index > endpoints[4] & index <= endpoints[5]] = 2
  
  # R=10, Y=0, A=1
  Y[index > endpoints[5] & index <= endpoints[6]] = 0
  A[index > endpoints[5] & index <= endpoints[6]] = 1
  R[index > endpoints[5] & index <= endpoints[6]] = 2
  
  # R=10, Y=1, A=0
  Y[index > endpoints[6] & index <= endpoints[7]] = 1
  A[index > endpoints[6] & index <= endpoints[7]] = 0
  R[index > endpoints[6] & index <= endpoints[7]] = 2
  
  # R=10, Y=1, A=1
  Y[index > endpoints[7] & index <= endpoints[8]] = 1
  A[index > endpoints[7] & index <= endpoints[8]] = 1
  R[index > endpoints[7] & index <= endpoints[8]] = 2
  
  # R=01, Y=0, A=0
  Y[index > endpoints[8] & index <= endpoints[9]] = 0
  A[index > endpoints[8] & index <= endpoints[9]] = 0
  R[index > endpoints[8] & index <= endpoints[9]] = 1
  
  # R=01, Y=0, A=1
  Y[index > endpoints[9] & index <= endpoints[10]] = 0
  A[index > endpoints[9] & index <= endpoints[10]] = 1
  R[index > endpoints[9] & index <= endpoints[10]] = 1
  
  # R=01, Y=1, A=0
  Y[index > endpoints[10] & index <= endpoints[11]] = 1
  A[index > endpoints[10] & index <= endpoints[11]] = 0
  R[index > endpoints[10] & index <= endpoints[11]] = 1
  
  # R=01, Y=1, A=1
  Y[index > endpoints[11] & index <= endpoints[12]] = 1
  A[index > endpoints[11] & index <= endpoints[12]] = 1
  R[index > endpoints[11] & index <= endpoints[12]] = 1
  
  # R=00, Y=0, A=0
  Y[index > endpoints[12] & index <= endpoints[13]] = 0
  A[index > endpoints[12] & index <= endpoints[13]] = 0
  R[index > endpoints[12] & index <= endpoints[13]] = 0
  
  # R=00, Y=0, A=1
  Y[index > endpoints[13] & index <= endpoints[14]] = 0
  A[index > endpoints[13] & index <= endpoints[14]] = 1
  R[index > endpoints[13] & index <= endpoints[14]] = 0
  
  # R=00, Y=1, A=0
  Y[index > endpoints[14] & index <= endpoints[15]] = 1
  A[index > endpoints[14] & index <= endpoints[15]] = 0
  R[index > endpoints[14] & index <= endpoints[15]] = 0
  
  # R=00, Y=1, A=1
  Y[index > endpoints[15] & index <= endpoints[16]] = 1
  A[index > endpoints[15] & index <= endpoints[16]] = 1
  R[index > endpoints[15] & index <= endpoints[16]] = 0
  
  ind = (R==3 & Y==0 & A==0)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(3,4), Sigma=matrix(c(0.5,0.1,0.1,0.5),2,2))
  ind = (R==3 & Y==0 & A==1)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(4,4), Sigma=matrix(c(0.5,0.2,0.2,0.5),2,2))
  ind = (R==3 & Y==1 & A==0)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(3,2), Sigma=matrix(c(0.4,0.1,0.1,0.4),2,2))
  ind = (R==3 & Y==1 & A==1)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(2,2), Sigma=matrix(c(0.5,0.1,0.1,0.5),2,2))
  
  ind = (R==2 & Y==0 & A==0)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(2,19/5), Sigma=matrix(c(0.5,0.1,0.1,0.5),2,2))
  ind = (R==2 & Y==0 & A==1)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(2,16/5), Sigma=matrix(c(0.5,0.2,0.2,0.5),2,2))
  ind = (R==2 & Y==1 & A==0)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(1,3/2), Sigma=matrix(c(0.4,0.1,0.1,0.4),2,2))
  ind = (R==2 & Y==1 & A==1)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(3,11/5), Sigma=matrix(c(0.5,0.1,0.1,0.5),2,2))
  
  ind = (R==1 & Y==0 & A==0)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(13/5,2), Sigma=matrix(c(0.5,0.1,0.1,0.5),2,2))
  ind = (R==1 & Y==0 & A==1)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(14/5,1), Sigma=matrix(c(0.5,0.2,0.2,0.5),2,2))
  ind = (R==1 & Y==1 & A==0)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(3.125,2.5), Sigma=matrix(c(0.4,0.1,0.1,0.4),2,2))
  ind = (R==1 & Y==1 & A==1)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(1.9,1.5), Sigma=matrix(c(0.5,0.1,0.1,0.5),2,2))
  
  ind = (R==0 & Y==0 & A==0)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(3,4), Sigma=matrix(c(0.5,0.1,0.1,0.5),2,2))
  ind = (R==0 & Y==0 & A==1)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(4,4), Sigma=matrix(c(0.5,0.2,0.2,0.5),2,2))
  ind = (R==0 & Y==1 & A==0)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(3,2), Sigma=matrix(c(0.4,0.1,0.1,0.4),2,2))
  ind = (R==0 & Y==1 & A==1)
  X[ind, ] = mvrnorm(n=sum(ind), mu=c(2,2), Sigma=matrix(c(0.5,0.1,0.1,0.5),2,2))
  
  Xobs = X
  Xobs[R==1, 1] = NA
  Xobs[R==2, 2] = NA
  Xobs[R==0, c(1,2)] = NA
  
  return (list("X" = X, "Xobs" = Xobs, "R" = R, "A" = A, "Y" = Y))
}

######### estimate Gaussian parameters for X | R, Y, A
estimate_params = function(Xobs, R, Y, A)
{
  
  ind = (R==3 & Y==0 & A==0)
  mu_R3_Y0_A0 = colMeans(Xobs[ind, ])
  cov_R3_Y0_A0 = cov(Xobs[ind, ])
  
  ind = (R==3 & Y==0 & A==1)
  mu_R3_Y0_A1 = colMeans(Xobs[ind, ])
  cov_R3_Y0_A1 = cov(Xobs[ind, ])
  
  ind = (R==3 & Y==1 & A==0)
  mu_R3_Y1_A0 = colMeans(Xobs[ind, ])
  cov_R3_Y1_A0 = cov(Xobs[ind, ])
  
  ind = (R==3 & Y==1 & A==1)
  mu_R3_Y1_A1 = colMeans(Xobs[ind, ])
  cov_R3_Y1_A1 = cov(Xobs[ind, ])
  
  ind = (R==2 & Y==0 & A==0)
  mu_R2_Y0_A0 = mean(Xobs[ind, 1])
  var_R2_Y0_A0 = var(Xobs[ind, 1])
  
  ind = (R==2 & Y==0 & A==1)
  mu_R2_Y0_A1 = mean(Xobs[ind, 1])
  var_R2_Y0_A1 = var(Xobs[ind, 1])
  
  ind = (R==2 & Y==1 & A==0)
  mu_R2_Y1_A0 = mean(Xobs[ind, 1])
  var_R2_Y1_A0 = var(Xobs[ind, 1])
  
  ind = (R==2 & Y==1 & A==1)
  mu_R2_Y1_A1 = mean(Xobs[ind, 1])
  var_R2_Y1_A1 = var(Xobs[ind, 1])

  ind = (R==1 & Y==0 & A==0)
  mu_R1_Y0_A0 = mean(Xobs[ind, 2])
  var_R1_Y0_A0 = var(Xobs[ind, 2])
  
  ind = (R==1 & Y==0 & A==1)
  mu_R1_Y0_A1 = mean(Xobs[ind, 2])
  var_R1_Y0_A1 = var(Xobs[ind, 2])
  
  ind = (R==1 & Y==1 & A==0)
  mu_R1_Y1_A0 = mean(Xobs[ind, 2])
  var_R1_Y1_A0 = var(Xobs[ind, 2])
  
  ind = (R==1 & Y==1 & A==1)
  mu_R1_Y1_A1 = mean(Xobs[ind, 2])
  var_R1_Y1_A1 = var(Xobs[ind, 2])
  
  means_R3 = array(NA, dim=c(2,2,2))
  means_R3[1,1,] = mu_R3_Y0_A0
  means_R3[1,2,] = mu_R3_Y0_A1
  means_R3[2,1,] = mu_R3_Y1_A0
  means_R3[2,2,] = mu_R3_Y1_A1
  
  covs_R3 = array(NA, dim=c(2,2,2,2))
  covs_R3[1,1,,] = cov_R3_Y0_A0
  covs_R3[1,2,,] = cov_R3_Y0_A1
  covs_R3[2,1,,] = cov_R3_Y1_A0
  covs_R3[2,2,,] = cov_R3_Y1_A1
  
  means_R2 = array(NA, dim=c(2,2))
  vars_R2 = array(NA, dim=c(2,2))
  means_R2[1,1] = mu_R2_Y0_A0
  means_R2[1,2] = mu_R2_Y0_A1
  means_R2[2,1] = mu_R2_Y1_A0
  means_R2[2,2] = mu_R2_Y1_A1
  vars_R2[1,1] = var_R2_Y0_A0
  vars_R2[1,2] = var_R2_Y0_A1
  vars_R2[2,1] = var_R2_Y1_A0
  vars_R2[2,2] = var_R2_Y1_A1
  
  means_R1 = array(NA, dim=c(2,2))
  vars_R1 = array(NA, dim=c(2,2))
  means_R1[1,1] = mu_R1_Y0_A0
  means_R1[1,2] = mu_R1_Y0_A1
  means_R1[2,1] = mu_R1_Y1_A0
  means_R1[2,2] = mu_R1_Y1_A1
  vars_R1[1,1] = var_R1_Y0_A0
  vars_R1[1,2] = var_R1_Y0_A1
  vars_R1[2,1] = var_R1_Y1_A0
  vars_R1[2,2] = var_R1_Y1_A1
  
  return (list("means_R3" = means_R3, "covs_R3" = covs_R3,
               "means_R2" = means_R2, "vars_R2" = vars_R2,
               "means_R1" = means_R1, "vars_R1" = vars_R1))
}

############# get empirical frequencies p(r,y,a)

estimate_empirical_probs = function(R,Y,A)
{
  total = length(R)

  empirical_probs = array(NA, dim=c(4,2,2))
  
  for (r in c(0,1,2,3))
  {
    for (y in c(0,1))
    {
      for (a in c(0,1))
      {
        empirical_probs[r+1,y+1,a+1] = sum(R==r & Y==y & A==a) / total
      }
    }
  }
  return (empirical_probs)
}

############# calculate complete odds according to CCMV P(R=r|X_r,Y,A) / P(R=1d|X_r,Y,A)

calculate_odds = function(Xobs, r, y, a, params, empirical_probs, rho)
{
  if (r == 0)
  {
    odds = empirical_probs[1,y+1,a+1] / empirical_probs[4,y+1,a+1]
    odds = odds * exp(Xobs[R==3,1]*rho + Xobs[R==3,2]*rho)
  } else if (r == 1)
  {
    odds = dnorm(Xobs[R==3, 2], mean=params$means_R1[y+1,a+1], sd=sqrt(params$vars_R1[y+1,a+1]))
    odds = odds / dnorm(Xobs[R==3, 2], mean=params$means_R3[y+1,a+1,2], sd=sqrt(params$covs_R3[y+1,a+1,2,2]))
    odds = odds * empirical_probs[2,y+1,a+1] / empirical_probs[4,y+1,a+1]
    
    odds = odds * exp(Xobs[R==3,1]*rho)
  } else if (r == 2)
  {
    odds = dnorm(Xobs[R==3, 1], mean=params$means_R2[y+1,a+1], sd=sqrt(params$vars_R2[y+1,a+1]))
    odds = odds / dnorm(Xobs[R==3, 1], mean=params$means_R3[y+1,a+1,1], sd=sqrt(params$covs_R3[y+1,a+1,1,1]))
    odds = odds * empirical_probs[3,y+1,a+1] / empirical_probs[4,y+1,a+1]
    
    odds = odds * exp(Xobs[R==3,2]*rho)
  } else if (r == 3)
  {
    odds = 1
  }
  
  return (odds)
}

###### calculate normal pdf from the observed data according to CCMV

calculate_norm_pdf = function(x1, x2, Xobs, R, r, y, a, X)
{
  if (r == 0)
  {
    norm_pdf = dmvnorm(cbind(x1,x2), mean = colMeans(Xobs[R==3 & Y==y & A==a, c(1,2)]), sigma = cov(Xobs[R==3 & Y==y & A==a, c(1,2)]))
  } else if (r == 1)
  {
    norm_pdf = dnorm(x1,
                     mean = mean(Xobs[R==3 & Y==y & A==a,1]) 
                     + cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2]) * 
                       sd(Xobs[R==3 & Y==y & A==a, 1]) / sd(Xobs[R==3 & Y==y & A==a, 2]) * (x2 - mean(Xobs[R==3 & Y==y & A==a,2])),
                     sd = sqrt((1-cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2])^2) * var(Xobs[R==3 & Y==y & A==a,1])))
    
    norm_pdf = norm_pdf * dnorm(x2, mean = mean(Xobs[R==1 & Y==y & A==a, 2]),
                                sd = sqrt(var(Xobs[R==1 & Y==y & A==a, 2])))
  } else if (r == 2)
  {
    norm_pdf = dnorm(x2,
                     mean = mean(Xobs[R==3 & Y==y & A==a,2]) 
                     + cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2]) * 
                       sd(Xobs[R==3 & Y==y & A==a, 2]) / sd(Xobs[R==3 & Y==y & A==a, 1]) * (x1 - mean(Xobs[R==3 & Y==y & A==a,1])),
                     sd = sqrt((1-cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2])^2) * var(Xobs[R==3 & Y==y & A==a,2])))
    
    norm_pdf = norm_pdf * dnorm(x1, mean = mean(Xobs[R==2 & Y==y & A==a, 1]),
                                sd = sqrt(var(Xobs[R==2 & Y==y & A==a, 1])))
  } else if (r == 3)
  {
    norm_pdf = dmvnorm(cbind(x1,x2), mean = colMeans(Xobs[R==3 & Y==y & A==a, c(1,2)]), sigma = cov(Xobs[R==3 & Y==y & A==a, c(1,2)]))
  }
  return (norm_pdf)
}


########## evaluate the propensity score P(A=1|X)
calc_prop_score = function(x1, x2)
{
  empirical_probs = estimate_empirical_probs(R, Y, A)
  numerator = 0
  denominator = 0
  for (r in c(0,1,2,3))
  {
    for (y in c(0,1))
    {
      numerator = numerator + calculate_norm_pdf(x1, x2, Xobs, R, r, y, 1, X) * empirical_probs[r+1,y+1,2]
    }
  }
  
  for (r in c(0,1,2,3))
  {
    for (y in c(0,1))
    {
      for (a in c(0,1))
      {
        denominator = denominator + calculate_norm_pdf(x1, x2, Xobs, R, r, y, a, X) * empirical_probs[r+1,y+1,a+1]
      }
    }
  }
  
  return (numerator/denominator)
}


########## evaluate the regression function E[Y|A=a,X]
calc_reg_function = function(x1, x2, a)
{
  empirical_probs = estimate_empirical_probs(R, Y, A)
  numerator = 0
  denominator = 0
  for (r in c(0,1,2,3))
  {
    for (y in c(0,1))
    {
      denominator = denominator + calculate_norm_pdf(x1, x2, Xobs, R, r, y, a, X) * empirical_probs[r+1,y+1,a+1]
    }
  }
  
  for (r in c(0,1,2,3))
  {
    numerator = numerator + calculate_norm_pdf(x1, x2, Xobs, R, r, 1, a, X) * empirical_probs[r+1,2,a+1]
  }
  
  return (numerator/denominator)
}

######## evaluate 1 / P(R=1_d|X,Y,A) according to CCMV and rho sensitivity parameter
calculate_inv_prop_score_R = function(Xobs, R, Y, A, r, params, empirical_probs, rho)
{
  inv_prop_score_R_y1_a1 = 0
  inv_prop_score_R_y0_a1 = 0
  inv_prop_score_R_y1_a0 = 0
  inv_prop_score_R_y0_a0 = 0
  for (r in c(0,1,2,3))
  {
    inv_prop_score_R_y1_a1 = inv_prop_score_R_y1_a1 + calculate_odds(Xobs, r, y=1, a=1, params, empirical_probs, rho)
    inv_prop_score_R_y0_a1 = inv_prop_score_R_y0_a1 + calculate_odds(Xobs, r, y=0, a=1, params, empirical_probs, rho)
    inv_prop_score_R_y1_a0 = inv_prop_score_R_y1_a0 + calculate_odds(Xobs, r, y=1, a=0, params, empirical_probs, rho)
    inv_prop_score_R_y0_a0 = inv_prop_score_R_y0_a0 + calculate_odds(Xobs, r, y=0, a=0, params, empirical_probs, rho)
  }
  
  inv_prop_score_R_a1 = inv_prop_score_R_y1_a1 * (Y[R==3]==1) + inv_prop_score_R_y0_a1 * (Y[R==3]==0)
  inv_prop_score_R_a0 = inv_prop_score_R_y1_a0 * (Y[R==3]==1) + inv_prop_score_R_y0_a0 * (Y[R==3]==0)
  
  inv_prop_score_R = inv_prop_score_R_a1 * (A[R==3]==1) + inv_prop_score_R_a0 * (A[R==3]==0)
  inv_prop_score_R
  
  return (inv_prop_score_R)
}

#### calculate all 3 IPW-R estimators, rho is a sensitivity parameter
estimateIPW = function(Xobs, R, Y, A, rho)
{
  params = estimate_params(Xobs, R, Y, A)
  empirical_probs = estimate_empirical_probs(R, Y, A)
  
  ## calculate weights according to CCMV and rho sensitivity parameter
  inv_prop_score_R = calculate_inv_prop_score_R(Xobs, R, Y, A, r, params, empirical_probs, rho)
  
  prop_score_A = calc_prop_score(Xobs[R==3,1], Xobs[R==3,2])
  
  ipw_R_ipw_A = sum(inv_prop_score_R[A[R==3]==1] / prop_score_A[A[R==3]==1] * Y[R==3 & A==1])/n - sum(inv_prop_score_R[A[R==3]==0] / (1-prop_score_A[A[R==3]==0]) * Y[R==3 & A==0])/n
  ipw_R_reg_A = sum( (calc_reg_function(Xobs[R==3,1], Xobs[R==3,2], 1) - calc_reg_function(Xobs[R==3,1], Xobs[R==3,2], 0)) * inv_prop_score_R) / n
  
  ipw_R_dr_A = ipw_R_ipw_A + ipw_R_reg_A
  ipw_R_dr_A = ipw_R_dr_A - sum(inv_prop_score_R[A[R==3]==1] / prop_score_A[A[R==3]==1] * calc_reg_function(Xobs[R==3 & A==1,1], Xobs[R==3 & A==1,2], 1))/n
  ipw_R_dr_A = ipw_R_dr_A + sum(inv_prop_score_R[A[R==3]==0] / (1-prop_score_A[A[R==3]==0]) * calc_reg_function(Xobs[R==3 & A==0,1], Xobs[R==3 & A==0,2], 0))/n
  
  return (list("ipw_R_ipw_A" = ipw_R_ipw_A, "ipw_R_reg_A" = ipw_R_reg_A, "ipw_R_dr_A" = ipw_R_dr_A, "inv_prop_score_R" = inv_prop_score_R))
}

####### calculate all 3 RA-R estimators, xi is a sensitivity parameter
estimateReg = function(Xobs, R, Y, A, M, xi)
{
  stacked = impute(Xobs, R, Y, A, M, xi)
  Y_imputed = stacked[,3]
  A_imputed = stacked[,4]
  
  prop_score_A_imputed = calc_prop_score(stacked[,1], stacked[,2])
  reg_R_ipw_A = ( sum(Y_imputed[A_imputed==1]/prop_score_A_imputed[A_imputed==1]) - sum(Y_imputed[A_imputed==0]/(1-prop_score_A_imputed[A_imputed==0])) ) / dim(stacked)[1]
  
  reg_R_reg_A = sum( calc_reg_function(stacked[,1], stacked[,2], 1) - calc_reg_function(stacked[,1], stacked[,2], 0) ) / dim(stacked)[1]
  
  reg_R_dr_A = reg_R_ipw_A + reg_R_reg_A
  reg_R_dr_A = reg_R_dr_A - (sum(calc_reg_function(stacked[A_imputed==1,1], stacked[A_imputed==1,2], 1)/prop_score_A_imputed[A_imputed==1]) - sum(calc_reg_function(stacked[A_imputed==0,1], stacked[A_imputed==0,2], 0)/(1-prop_score_A_imputed[A_imputed==0]))) / dim(stacked)[1]

  return (list("reg_R_ipw_A" = reg_R_ipw_A, "reg_R_reg_A" = reg_R_reg_A, "reg_R_dr_A" = reg_R_dr_A))  
}

######## perform stacked imputation according to CCMV and using xi sensitivity parameter
impute = function(Xobs, R, Y, A, M, xi)
{
  params = estimate_params(Xobs, R, Y, A)
  
  stacked = matrix(, nrow=0, ncol=4)
  
  for (m in 1:M)
  {
    imputed = cbind(Xobs, Y, A)
    for (y in c(0,1))
    {
      for (a in c(0,1))
      {
        #imputed[A==a & Y==y & R==0, c(1,2)] = rmvnorm(n=sum(A==a & Y==y & R==0), mean = colMeans(Xobs[R==3 & Y==y & A==a, c(1,2)]), sigma = cov(Xobs[R==3 & Y==y & A==a, c(1,2)]))
        
        ## R = 0
        
        mu = mean(Xobs[R==3 & Y==y & A==a, 1])
        sigma = sd(Xobs[R==3 & Y==y & A==a, 1])
        
        mu = mu / (2*sigma^2*xi + 1)
        sigma = sigma / sqrt(2*sigma^2*xi + 1)
        
        imputed[A==a & Y==y & R==0, 1] = rnorm(n=sum(A==a & Y==y & R==0), mean = mu, sd = sigma)
        
        mu = mean(Xobs[R==3 & Y==y & A==a,2]) + 
          cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2]) * sd(Xobs[R==3 & Y==y & A==a, 2]) / sd(Xobs[R==3 & Y==y & A==a, 1]) * (imputed[A==a & Y==y & R==0, 1] - mean(Xobs[R==3 & Y==y & A==a,1]))
        sigma = sqrt((1-cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2])^2) * var(Xobs[R==3 & Y==y & A==a,2]))
        
        mu = mu / (2*sigma^2*xi + 1)
        sigma = sigma / sqrt(2*sigma^2*xi + 1)
        
        imputed[A==a & Y==y & R==0, 2] = rnorm(n=sum(A==a & Y==y & R==0),
                                               mean = mu,
                                               sd = sigma)
        ## R = 1
        
        mu = mean(Xobs[R==3 & Y==y & A==a,1]) +
          cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2]) * sd(Xobs[R==3 & Y==y & A==a, 1]) / sd(Xobs[R==3 & Y==y & A==a, 2]) * (Xobs[A==a & Y==y & R==1, 2] - mean(Xobs[R==3 & Y==y & A==a,2]))
        sigma = sqrt((1-cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2])^2) * var(Xobs[R==3 & Y==y & A==a,1]))
        
        mu = mu / (2*sigma^2*xi + 1)
        sigma = sigma / sqrt(2*sigma^2*xi + 1)
        
        imputed[A==a & Y==y & R==1, 1] = rnorm(n=sum(A==a & Y==y & R==1),
                                               mean = mu,
                                               sd = sigma)
        
        ## R = 2
        
        mu = mean(Xobs[R==3 & Y==y & A==a,2]) +
          cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2]) * sd(Xobs[R==3 & Y==y & A==a, 2]) / sd(Xobs[R==3 & Y==y & A==a, 1]) * (Xobs[A==a & Y==y & R==2, 1] - mean(Xobs[R==3 & Y==y & A==a,1]))
        sigma = sqrt((1-cor(Xobs[R==3 & Y==y & A==a,1], Xobs[R==3 & Y==y & A==a,2])^2) * var(Xobs[R==3 & Y==y & A==a,2]))
        
        mu = mu / (2*sigma^2*xi + 1)
        sigma = sigma / sqrt(2*sigma^2*xi + 1)
        
        imputed[A==a & Y==y & R==2, 2] = rnorm(n=sum(A==a & Y==y & R==2),
                                               mean = mu,
                                               sd = sigma)
      }
    }
    
    stacked = rbind(stacked, imputed)
  }
 return (stacked) 
}
