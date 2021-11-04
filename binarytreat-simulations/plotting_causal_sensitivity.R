#################### plotting for causal sensitivity analysis simulation

library(tidyverse)
library(latex2exp)
require(plotrix)
require(matrixStats)

xi = seq(0,4,length.out=101)
rho = seq(-2,1,length.out=101)

ipw_means = matrix(NA, 3, 101)
ipw_sds = matrix(NA, 3, 101)

ipw_upper_quantiles = matrix(NA, 3, 101)
ipw_lower_quantiles = matrix(NA, 3, 101)

mi_means = matrix(NA, 3, 101)
mi_sds = matrix(NA, 3, 101)

### load data

for (kk in 1:101)
{
    

  load(paste("causal_sim_MI_sensitivity_",xi[kk],".RData", sep=""))
  load(paste("causal_sim_ipw_sensitivity_",rho[kk],".RData", sep=""))
  
  ipw_means[, kk] = colMeans(ipw_R_matrix)
  ipw_sds[, kk] = colSds(ipw_R_matrix)
  
  mi_means[, kk] = colMeans(reg_R_matrix)
  mi_sds[, kk] = colSds(reg_R_matrix)
  
  ipw_upper_quantiles[, kk] = colQuantiles(ipw_R_matrix, probs = rep(0.975, 4))[,1]
  ipw_lower_quantiles[, kk] = colQuantiles(ipw_R_matrix, probs = rep(0.025, 4))[,1]
      
  
}

for (ii in c(1,2,3))
{
  ### titles
  
  if (ii == 1)
  {
    text1 = paste('$(IPW-R, \\, IPW-A) \\, - \\, ATE \\, Simulation', sep="")
    text2 = paste('$(RA-R, \\, IPW-A) \\, - \\, ATE \\, Simulation', sep="")
  } else if (ii == 2)
  {
    text1 = paste('$(IPW-R, \\,  RA-A) \\, - \\, ATE \\, Simulation', sep="")
    text2 = paste('$(RA-R, \\, RA-A) \\, - \\, ATE \\, Simulation', sep="")
  } else if (ii == 3)
  {
    text1 = paste('$(IPW-R, \\,  DR-A) \\, - \\, ATE \\, Simulation', sep="")
    text2 = paste('$(RA-R, \\, DR-A) \\, - \\, ATE \\, Simulation', sep="")
  }
  
  ystart = -2
  ylast = 2.8
  
  ### plot IPW-R
  
  data = data.frame(xi = seq(0,4,length.out=101), rho = seq(-2,1,length.out=101), y_ipw = ipw_means[ii,], upper_ipw = ipw_means[ii,] + ipw_sds[ii,]*1.96,
        lower_ipw = ipw_means[ii,] - ipw_sds[ii,]*1.96,
        upper_quantiles = ipw_upper_quantiles[ii,],
        lower_quantiles = ipw_lower_quantiles[ii,],
        y_mi = mi_means[ii,], upper_mi = mi_means[ii,] + mi_sds[ii,]*1.96,
        lower_mi = mi_means[ii,] - mi_sds[ii,]*1.96)
  
  pdf(file=paste("sim_sensitivity_causal",ii,"_ipw.pdf", sep=""), width=8, height=4.5)
  print(ggplot(data[1:85,], aes(rho[1:85], y_ipw[1:85])) +
        geom_point() +
        geom_smooth() +
        geom_ribbon(aes(ymin = lower_quantiles[1:85], ymax = upper_quantiles[1:85], x = rho[1:85]), alpha=0.3) +
        geom_hline(yintercept=0.269, linetype='dotted', col = 'red') +
        ylim(ystart, ylast) +
        labs(title = TeX(text1)) +
        xlab(TeX('$\\rho$')) +
        ylab("ATE") +
        theme(plot.title = element_text(size=22), 
              axis.text.x = element_text(size=18), 
              axis.title.x = element_text(size=18),
              axis.text.y = element_text(size=18),
              axis.title.y = element_text(size=18)))
      
  dev.off()
  
  ### plot RA-R
  
  pdf(file=paste("sim_sensitivity_causal",ii,"_ra.pdf", sep=""), width=8, height=4.5)
  print(ggplot(data, aes(xi, y_mi)) +
    geom_point() +
    geom_smooth() +
    geom_ribbon(aes(ymin = lower_mi, ymax = upper_mi, x = xi), alpha=0.3) +
    geom_hline(yintercept=0.269, linetype='dotted', col = 'red') +
    ylim(ystart, ylast) +
    labs(title = TeX(text2)) +
    xlab(TeX('$\\xi$')) +
    ylab("ATE") +
    theme(plot.title = element_text(size=22), 
          axis.text.x = element_text(size=18), 
          axis.title.x = element_text(size=18),
          axis.text.y = element_text(size=18),
          axis.title.y = element_text(size=18)))
      
  dev.off()
 
}