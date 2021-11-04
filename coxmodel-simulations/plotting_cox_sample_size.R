require(plotrix)
require(matrixStats)

load("survival_samplesize_simulation.RData")

#################### plotting

sample_sizes = c(350, 500, 1000, 2000)

mi_beta_means = MI_estimates[,1000,]
#mi_beta_means = matrix(NA, 4, 2)
#for (kk in 1:4)
#  mi_beta_means[kk,] = colMeans(MI_estimates[kk,,])

mi_beta_sd = matrix(NA, 4, 2)
for (kk in 1:4)
  mi_beta_sd[kk,] = colSds(MI_estimates[kk,,])

ipw_beta_means = ipw_estimates[,1000,]
#ipw_beta_means = matrix(NA, 4, 2)
#for (kk in 1:4)
#  ipw_beta_means[kk,] = ipw_estimates[,1,]#colMeans(ipw_estimates[kk,,])

ipw_beta_sd = matrix(NA, 4, 2)
for (kk in 1:4)
  ipw_beta_sd[kk,] = colSds(ipw_estimates[kk,,])


transformed_means = matrix(NA, 4, 2)
for (kk in 1:4)
  transformed_means[kk,] = colMeans(transformed_estimates[kk,,])

transformed_sd = matrix(NA, 4, 2)
for (kk in 1:4)
  transformed_sd[kk,] = colSds(transformed_estimates[kk,,])


################################# plot all 3


##### beta1

pdf(file="beta1_final.pdf", width=5, height=4)
plotCI(c(0.9, 1.9, 2.9, 3.9), ipw_beta_means[,1], ui=ipw_beta_means[,1]+1.96*ipw_beta_sd[,1], li=ipw_beta_means[,1]-1.96*ipw_beta_sd[,1],col="orange",xlim=c(0.5,4.5),
       ylim=c(-1.5,1),xlab="Sample Size",ylab="",xaxt="n",pch=19,
       cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
plotCI(c(1.1, 2.1, 3.1, 4.1), mi_beta_means[,1], ui=mi_beta_means[,1]+1.96*mi_beta_sd[,1], 
       li=mi_beta_means[,1]-1.96*mi_beta_sd[,1],col="blue",add=TRUE,pch=19,
       cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
plotCI(c(1, 2, 3, 4), transformed_means[,1], ui=transformed_means[,1]+1.96*transformed_sd[,1],
       li=transformed_means[,1]-1.96*transformed_sd[,1],col="magenta",add=TRUE,pch=19,
       cex.lab=1.5, cex.axis=1.1, cex.main=1.5)

abline(h = -0.5, lty = 2, lwd = 1.5)
title(expression(beta[1]), cex.main=1.5)
legend("topright", pch = c(19, 19, 19), 
       col = c("orange", "magenta", "blue"), 
       legend = c("IPW", "Transformed MLE", "RA"))
axis(1, at = c(1,2,3,4), labels=sample_sizes, cex.axis=1.1)

dev.off()

###### beta2

pdf(file="beta2_final.pdf", width=5, height=4)

plotCI(c(0.9, 1.9, 2.9, 3.9), ipw_beta_means[,2], ui=ipw_beta_means[,2]+1.96*ipw_beta_sd[,2], li=ipw_beta_means[,2]-1.96*ipw_beta_sd[,2],col="orange",xlim=c(0.5,4.5),
       ylim=c(1,4),xlab="Sample Size",ylab="",xaxt="n",pch=19,
       cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
plotCI(c(1.1, 2.1, 3.1, 4.1), mi_beta_means[,2], ui=mi_beta_means[,2]+1.96*mi_beta_sd[,2], 
       li=mi_beta_means[,2]-1.96*mi_beta_sd[,2],col="blue",add=TRUE,pch=19,
       cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
plotCI(c(1, 2, 3, 4), transformed_means[,2], ui=transformed_means[,2]+1.96*transformed_sd[,2],
       li=transformed_means[,2]-1.96*transformed_sd[,2],col="magenta",add=TRUE,pch=19,
       cex.lab=1.5, cex.axis=1.1, cex.main=1.5)


abline(h = 2, lty = 2, lwd = 1.5)
title(expression(beta[2]), cex.main=1.5)
legend("topright", pch = c(19, 19, 19), 
       col = c("orange", "magenta", "blue"), 
       legend = c("IPW", "Transformed MLE", "RA"))
axis(1, at = c(1,2,3,4), labels=sample_sizes, cex.axis=1.1)


dev.off()

