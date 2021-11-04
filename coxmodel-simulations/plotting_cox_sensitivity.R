library(tidyverse)
library(latex2exp)
require(plotrix)
require(matrixStats)

ii = 2
zeta = seq(0,1,length.out=101)
rho = seq(-3,3,length.out=101)

ipw_means = matrix(NA, 2, 101)
ipw_sds = matrix(NA, 2, 101)

mi_means = matrix(NA, 2, 101)
mi_sds = matrix(NA, 2, 101)

for (kk in 1:101)
{
    

  load(paste("survival_sim_MI_sensitivity_",zeta[kk],".RData", sep=""))
  load(paste("survival_sim_ipw_sensitivity_",rho[kk],".RData", sep=""))
  ipw_means[, kk] = colMeans(ipw_estimates[1,,])
  ipw_sds[, kk] = colSds(ipw_estimates[1,,])
  
  mi_means[, kk] = colMeans(MI_estimates[1,,])
  mi_sds[, kk] = colSds(MI_estimates[1,,])
      
  
}
    
text1 = paste('$\\beta_', ii, '\\, Simulation \\, - \\, IPW \\, Estimator', sep="")
text2 = paste('$\\beta_', ii, '\\, Simulation \\, - \\, RA \\, Estimator', sep="")

if (ii==1)
{
  yend = -1.25
  ylast = 0.5
  truth = -0.5
} else
{
  yend = 0
  ylast = 3
  truth = 2
}

data = data.frame(zeta = seq(0,1,length.out=101), rho = seq(-3,3,length.out=101), y_ipw = ipw_means[ii,], upper_ipw = ipw_means[ii,] + ipw_sds[ii,]*1.96,
      lower_ipw = ipw_means[ii,] - ipw_sds[ii,]*1.96,
      y_mi = mi_means[ii,], upper_mi = mi_means[ii,] + mi_sds[ii,]*1.96,
      lower_mi = mi_means[ii,] - mi_sds[ii,]*1.96)

pdf(file=paste("sim_sensitivity_beta",ii,"_ipw.pdf", sep=""), width=8, height=4.5)
ggplot(data, aes(rho, y_ipw)) +
      geom_point() +
      geom_smooth() +
      geom_ribbon(aes(ymin = lower_ipw, ymax = upper_ipw, x = rho), alpha=0.3) +
      geom_hline(yintercept=truth, linetype='dotted', col = 'red') +
      ylim(yend, ylast) +
      labs(title = TeX(text1)) +
      xlab(TeX('$\\rho$')) +
      ylab(TeX(paste('$\\beta_',ii,sep=""))) +
      theme(plot.title = element_text(size=22), 
        axis.text.x = element_text(size=18), 
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=18))
    #labs(title = paste("beta",ii," Age Group ",start_age,"-",start_age+4,sep=""), y = "beta"))
    
dev.off()

pdf(file=paste("sim_sensitivity_beta",ii,"_ra.pdf", sep=""), width=8, height=4.5)
ggplot(data, aes(zeta, y_mi)) +
  geom_point() +
  geom_smooth() +
  geom_ribbon(aes(ymin = lower_mi, ymax = upper_mi, x = zeta), alpha=0.3) +
  geom_hline(yintercept=truth, linetype='dotted', col = 'red') +
  ylim(yend, ylast) +
  labs(title = TeX(text2)) +
  xlab(TeX('$\\zeta$')) +
  ylab(TeX(paste('$\\beta_',ii,sep=""))) +
  theme(plot.title = element_text(size=22), 
        axis.text.x = element_text(size=18), 
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title.y = element_text(size=18))
    
dev.off()
 
