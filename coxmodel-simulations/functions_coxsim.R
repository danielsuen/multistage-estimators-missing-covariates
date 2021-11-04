library(simsurv)
library(flexsurv)
library(MASS)
library(nnet)
library(bindata)

############## output of for loop
comb = function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

############### generate data

gen_complete_data = function(n)
{
  p1 = 0.5
  p2 = 0.3
  r = 0.3
  
  X = rmvbin(n, c(p1,p2), bincorr=(1-r)*diag(2)+r)
  
  beta = c(-0.5,2)
  
  shape = 1 # shape - exponential
  alpha1 = 1
  alpha2 = 2 # scale
  
  # generate T, C according to exponential; find Y = min(T,C)
  T = (- log(runif(n)) / (alpha1 * exp(X %*% beta)))^(1 / shape)
  C = rweibull(n, shape, 1/alpha2)
  delta = as.numeric(T <= C)
  Y = pmin(Tlat, C)
  
  ## generate missing from CCMV and multinomial logistic
  
  eta_10 = c(-1,0.5,1,1)
  norm_prob_R10 = exp(eta_10[1] + eta_10[2]*X[,1] + eta_10[3]*Y + eta_10[4]*delta)
  
  eta_01 = c(0,-0.5,0.75,0.5)
  norm_prob_R01 = exp(eta_01[1] + eta_01[2]*X[,2] + eta_01[3]*Y + eta_01[4]*delta)
  
  #norm_prob_R10 = exp(-1 + X[,1] + X[,2] + 0.5*Y + 1*delta)
  #norm_prob_R01 = exp(-0.5 + X[,1] + X[,2] + -0.5*Y + 0.5*delta)
  
  cc_prob = 1/(1 + norm_prob_R10 + norm_prob_R01)
  prob_R01 = norm_prob_R01 * cc_prob
  prob_R10 = norm_prob_R10 * cc_prob
  
  prob_R = cbind(prob_R01, prob_R10, cc_prob)
  R = rep(0, n)
  for (ii in 1:n)
  {
    R[ii] = sample(3,1,prob=prob_R[ii,])
  }
  
  X_obs = X
  X_obs[R==1,1] = NA
  X_obs[R==2,2] = NA
  X_cc = X[R==3,]
  Y_cc = Y[R==3]
  delta_cc = delta[R==3]

  X_cc = X_obs[R==3,]
  
  return (list("X" = X, "R" = R, "delta" = delta, "Y" = Y, "eta_10" = eta_10, "eta_01" = eta_01,
               "X_obs" = X_obs, "X_cc" = X_cc, "Y_cc" = Y_cc, "delta_cc" = delta_cc))
}

############### helper functions

prob_y_delta_given_x = function(y, delta, x, gamma, alpha2)
{
  return (gamma^delta * alpha2^(1-delta) * exp(-y*(gamma+alpha2)))
}

exp_linear_R01 = function(x2, y, delta, eta)
{
  exp(eta[1] + eta[2]*x2 + eta[3]*y + eta[4]*delta)
}
exp_linear_R10 = function(x1, y, delta, eta)
{
  exp(eta[1] + eta[2]*x1 + eta[3]*y + eta[4]*delta)
}
selection_R11 = function(x1, x2, y, delta)
{
  return (1 / (1 + exp_linear_R01(x2,y,delta,eta_01_est) + exp_linear_R10(x1,y,delta,eta_10_est)))
}

# integrable function
f = function(y, gamma, alpha2, x, delta)
{
  expression = prob_y_delta_given_x(y, delta, x, gamma, alpha2) #gamma * exp(-(gamma + alpha2)*y)
  expression = expression / (1 + exp_linear_R01(x[2], y, delta, eta_01_est) + exp_linear_R10(x[1], y, delta, eta_10_est))
  expression
}

################## EM algorithm

EM = function(X_obs, Y, delta, R, eta_10_est, eta_01_est, cc_empirical_prob)
{
  beta_init = c(1, 1)
  alpha1_init = 1
  alpha2_init = 3
  gamma_init_11 = alpha1_init * exp(beta_init%*%c(1,1))
  gamma_init_10 = alpha1_init * exp(beta_init%*%c(1,0))
  gamma_init_01 = alpha1_init * exp(beta_init%*%c(0,1))
  gamma_init_00 = alpha1_init * exp(beta_init%*%c(0,0))
  
  weights = matrix(0, nrow=n, ncol=4)
  # 11, 10, 01, 00

  
  alpha2_curr = alpha2_init
  gamma_curr_11 = gamma_init_11
  gamma_curr_10 = gamma_init_10
  gamma_curr_01 = gamma_init_01
  gamma_curr_00 = gamma_init_00
  beta_curr = beta_init
  
  gamma_old = c(0,0,0,0)
  
  for (ii in 1:100)
  {
    
    ######## E-step, calculate weights
    
    ## R = 11
    
    ind_11 = X_obs[,1]%in%1 & X_obs[,2]%in%1
    weights[ind_11, ] = rep(c(1,0,0,0), each=sum(ind_11))
    ind_10 = X_obs[,1]%in%1 & X_obs[,2]%in%0
    weights[ind_10, ] = rep(c(0,1,0,0), each=sum(ind_10))
    ind_01 = X_obs[,1]%in%0 & X_obs[,2]%in%1
    weights[ind_01, ] = rep(c(0,0,1,0), each=sum(ind_01))
    ind_00 = X_obs[,1]%in%0 & X_obs[,2]%in%0
    weights[ind_00, ] = rep(c(0,0,0,1), each=sum(ind_00))
    
    ## R = 01
  
    prob_R_X_11 = integrate(f, lower = 0, upper = Inf, gamma = gamma_curr_11, alpha2 = alpha2_curr, x = c(1,1), delta = 1)$value
    prob_R_X_11 = prob_R_X_11 + integrate(f, lower = 0, upper = Inf, gamma = gamma_curr_11, alpha2 = alpha2_curr, x = c(1,1), delta = 0)$value
    
    prob_R_X_10 = integrate(f, lower = 0, upper = Inf, gamma = gamma_curr_10, alpha2 = alpha2_curr, x = c(1,0), delta = 1)$value
    prob_R_X_10 = prob_R_X_10 + integrate(f, lower = 0, upper = Inf, gamma = gamma_curr_10, alpha2 = alpha2_curr, x = c(1,0), delta = 0)$value
    
    prob_R_X_01 = integrate(f, lower = 0, upper = Inf, gamma = gamma_curr_01, alpha2 = alpha2_curr, x = c(0,1), delta = 1)$value
    prob_R_X_01 = prob_R_X_01 + integrate(f, lower = 0, upper = Inf, gamma = gamma_curr_01, alpha2 = alpha2_curr, x = c(0,1), delta = 0)$value
    
    prob_R_X_00 = integrate(f, lower = 0, upper = Inf, gamma = gamma_curr_00, alpha2 = alpha2_curr, x = c(0,0), delta = 1)$value
    prob_R_X_00 = prob_R_X_00 + integrate(f, lower = 0, upper = Inf, gamma = gamma_curr_00, alpha2 = alpha2_curr, x = c(0,0), delta = 0)$value
    
    prob_R_X_matrix = matrix(NA, 2, 2)
    prob_R_X_matrix[1, 1] = prob_R_X_00
    prob_R_X_matrix[1, 2] = prob_R_X_01
    prob_R_X_matrix[2, 1] = prob_R_X_10
    prob_R_X_matrix[2, 2] = prob_R_X_11
    # 01, X1 missing
    
    index = (R==1) & (X[,2]==1)
    impute_prob_num = cc_empirical_prob[2,2] / prob_R_X_11 * selection_R11(x1 = 1, x2 = X[index,2], Y[index], delta[index]) * prob_y_delta_given_x(Y[index], delta[index], x=c(1,1), gamma_curr_11, alpha2_curr)
    impute_prob_denom = impute_prob_num
    impute_prob_denom = impute_prob_denom + cc_empirical_prob[1,2] / prob_R_X_01 * selection_R11(x1 = 0, x2 = X[index,2], Y[index], delta[index]) * prob_y_delta_given_x(Y[index], delta[index], x=c(0,1), gamma_curr_01, alpha2_curr)
    
    weights[index, c(1,3)] = c(impute_prob_num/impute_prob_denom, 1-impute_prob_num/impute_prob_denom)
    
    
    index = (R==1) & (X[,2]==0)
    impute_prob_num = cc_empirical_prob[2,1] / prob_R_X_10 * selection_R11(x1 = 1, x2 = X[index,2], Y[index], delta[index]) * prob_y_delta_given_x(Y[index], delta[index], x=c(1,0), gamma_curr_10, alpha2_curr)
    impute_prob_denom = impute_prob_num
    impute_prob_denom = impute_prob_denom + cc_empirical_prob[1,1] / prob_R_X_00 * selection_R11(x1 = 0, x2 = X[index,2], Y[index], delta[index]) * prob_y_delta_given_x(Y[index], delta[index], x=c(0,0), gamma_curr_00, alpha2_curr)
    
    weights[index, c(2,4)] = c(impute_prob_num/impute_prob_denom, 1-impute_prob_num/impute_prob_denom)
    
    
    # 10, X2 missing
    
    index = (R==2) & (X[,1]==1)
    impute_prob_num = cc_empirical_prob[2,2] / prob_R_X_11 * selection_R11(x1 = X[index,1], x2 = 1, Y[index], delta[index]) * prob_y_delta_given_x(Y[index], delta[index], x=c(1,1), gamma_curr_11, alpha2_curr)
    impute_prob_denom = impute_prob_num
    impute_prob_denom = impute_prob_denom + cc_empirical_prob[2,1] / prob_R_X_10 * selection_R11(x1 = X[index,1], x2 = 0, Y[index], delta[index]) * prob_y_delta_given_x(Y[index], delta[index], x=c(1,0), gamma_curr_10, alpha2_curr)
    
    weights[index, c(1,2)] = c(impute_prob_num/impute_prob_denom, 1-impute_prob_num/impute_prob_denom)
    
    
    index = (R==2) & (X[,1]==0)
    impute_prob_num = cc_empirical_prob[1,2] / prob_R_X_01 * selection_R11(x1 = X[index,1], x2 = 1, Y[index], delta[index]) * prob_y_delta_given_x(Y[index], delta[index], x=c(0,1), gamma_curr_01, alpha2_curr)
    impute_prob_denom = impute_prob_num
    impute_prob_denom = impute_prob_denom + cc_empirical_prob[1,1] / prob_R_X_00 * selection_R11(x1 = X[index,1], x2 = 0, Y[index], delta[index]) * prob_y_delta_given_x(Y[index], delta[index], x=c(0,0), gamma_curr_00, alpha2_curr)
    
    weights[index, c(3,4)] = c(impute_prob_num/impute_prob_denom, 1-impute_prob_num/impute_prob_denom)
    
    
    ######## M-step
    
    alpha2_curr = (n - sum(delta)) / sum(Y)
    
    gamma_curr_00 = sum(delta*weights[,4]) / sum(Y*weights[,4])
    gamma_curr_01 = sum(delta*weights[,3]) / sum(Y*weights[,3])
    gamma_curr_10 = sum(delta*weights[,2]) / sum(Y*weights[,2])
    gamma_curr_11 = sum(delta*weights[,1]) / sum(Y*weights[,1])
    
    gamma_est = c(gamma_curr_00, gamma_curr_10, gamma_curr_01, gamma_curr_11)
    
    # check convergence
    if (sqrt(sum(gamma_old-gamma_est)^2) < 0.001)
      break
    
    gamma_old = gamma_est
  }
  
  return (list("gamma_est"=gamma_est, "weights"=weights))
  
}

############### estimate logistic

estimate_logistic = function(X_obs, Y, delta, R)
{
  class = as.factor(R[R==2 | R==3] - 1)
  class <- relevel(class, ref = 2)
  model = glm(class ~ X_obs[R==2 | R==3,1] + Y[R==2 | R==3] + delta[R==2 | R==3], family="binomial")
  eta_10_est = coef(model)
  
  class = R[R==1 | R==3]
  class[class==3] = 2
  class = as.factor(class)
  class = relevel(class, ref = 2)
  model = glm(class ~ X_obs[R==1 | R==3,2] + Y[R==1 | R==3] + delta[R==1 | R==3], trace = FALSE, family="binomial")
  eta_01_est = coef(model)
  
  return (list("eta_10_est" = eta_10_est, "eta_01_est" = eta_01_est))
}

############# IPW weights

get_IPW_weights = function(X_cc, Y_cc, delta_cc, eta_10, eta_01, rho)
{
  lambda_1 = (1 + exp(eta_10[1] + eta_10[2]*X_cc[,1] + eta_10[3]*Y_cc + eta_10[4] + rho*X_cc[,2]) + exp(eta_01[1] + eta_01[2]*X_cc[,2] + eta_01[3]*Y_cc + eta_01[4] + rho*X_cc[,1])) * delta_cc
  lambda_0 = (1 + exp(eta_10[1] + eta_10[2]*X_cc[,1] + eta_10[3]*Y_cc + rho*X_cc[,2]) + exp(eta_01[1] + eta_01[2]*X_cc[,2] + eta_01[3]*Y_cc + rho*X_cc[,1])) * (1-delta_cc)
  lambda = lambda_0 + lambda_1
  
  return (lambda)
}

get_empirical_probs = function(X_obs)
{
  ## return p(x1, x2, R=11)
  cc_empirical_prob = matrix(NA, 2, 2)
  # row = X1, col = X2
  
  cc_empirical_prob[1,1] = sum(R==3 & X_obs[,1]%in%0 & X_obs[,2]%in%0) / n
  cc_empirical_prob[1,2] = sum(R==3 & X_obs[,1]%in%0 & X_obs[,2]%in%1) / n
  cc_empirical_prob[2,1] = sum(R==3 & X_obs[,1]%in%1 & X_obs[,2]%in%0) / n
  cc_empirical_prob[2,2] = sum(R==3 & X_obs[,1]%in%1 & X_obs[,2]%in%1) / n
  
  return (cc_empirical_prob)
}

#################### imputation

doImputation = function(X_obs, Y, delta, R, cc_empirical_prob, EM_output, zeta, M)
{

  weights = EM_output$weights
  
  weights = weights * (1-zeta) / (1 - weights*zeta - (1-zeta)*(1-weights))
  
  #### perform stacked imputation
  
  stacked = matrix(, nrow=0, ncol=4)
  
  for (m in 1:M)
  {
    imputed = cbind(X_obs, Y, delta)
    # impute X1, given X2 = 1
    ind = ((R==1) & X[,2]==1)
    imputed[ind, 1] = rbinom(n=sum(ind), size=1, prob=weights[ind,1])
    
    # impute X1, given X2 = 0
    ind = ((R==1) & X[,2]==0)
    imputed[ind, 1] = rbinom(n=sum(ind), size=1, prob=weights[ind,2])
    
    # impute X2, given X1 = 1
    ind = ((R==2) & X[,1]==1)
    imputed[ind, 2] = rbinom(n=sum(ind), size=1, prob=weights[ind,1])
    
    # impute X2, given X1 = 0
    ind = ((R==2) & X[,1]==0)
    imputed[ind, 2] = rbinom(n=sum(ind), size=1, prob=weights[ind,3])
    
    stacked = rbind(stacked, imputed)
  }
  
  stacked = data.frame(stacked)
  colnames(stacked) = c("X1", "X2", "Y", "delta")
  
  return(stacked)
}