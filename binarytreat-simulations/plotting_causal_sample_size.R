#################### plotting for causal sample size simulation
load('causal_samplesize_simulation.Rdata')

sample_sizes = c(1000, 2000, 5000, 10000)

ipw_R_matrix_means = array(NA, c(4, 3))
reg_R_matrix_means = array(NA, c(4, 3))

ipw_R_matrix_sds = array(NA, c(4, 3))
reg_R_matrix_sds = array(NA, c(4, 3))

for (kk in 1:4)
  ipw_R_matrix_means[kk,] = colMeans(ipw_R_matrix[kk,,])

for (kk in 1:4)
  ipw_R_matrix_sds[kk,] = colSds(ipw_R_matrix[kk,,])

for (kk in 1:4)
  reg_R_matrix_means[kk,] = colMeans(reg_R_matrix[kk,,])

for (kk in 1:4)
  reg_R_matrix_sds[kk,] = colSds(reg_R_matrix[kk,,])

require(plotrix)
require(matrixStats)

##### ipw

pdf(file="ipw_R_causalsim.pdf", width=5, height=4)

plotCI(c(0.9, 1.9, 2.9, 3.9), ipw_R_matrix_means[,1], ui=rowQuantiles(ipw_R_matrix[,,1], probs = rep(0.975, 4))[,1], li=rowQuantiles(ipw_R_matrix[,,1], probs = rep(0.025, 4))[,1],col="magenta",xlim=c(0.5,5),
      ylim=c(-0.2,1),xlab="Sample Size",ylab="",xaxt="n",pch=19,
      cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
plotCI(c(1, 2, 3, 4), ipw_R_matrix_means[,2], ui=rowQuantiles(ipw_R_matrix[,,2], probs = rep(0.975, 4))[,1],
      li=rowQuantiles(ipw_R_matrix[,,2], probs = rep(0.025, 4))[,1],col="darkblue",add=TRUE,pch=19,
      cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
plotCI(c(1.1, 2.1, 3.1, 4.1), ipw_R_matrix_means[,3], ui=rowQuantiles(ipw_R_matrix[,,3], probs = rep(0.975, 4))[,1],
      li=rowQuantiles(ipw_R_matrix[,,3], probs = rep(0.025, 4))[,1],col="brown",add=TRUE,pch=19,
      cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
abline(h = 0.269, lty = 2, lwd = 1.5)
title("IPW Estimators")
legend("topright", pch = c(19, 19, 19), 
       col = c("magenta", "darkblue", "brown"), 
       legend = c("IPW-R, IPW-A", "IPW-R, RA-A", "IPW-R, DR-A"),
       cex = 0.75)
axis(1, at = c(1,2,3,4), labels=sample_sizes, cex.axis=1.1)
dev.off()

###### ra


pdf(file="reg_R_causalsim.pdf", width=5, height=4)
plotCI(c(0.9, 1.9, 2.9, 3.9), reg_R_matrix_means[,1], ui=reg_R_matrix_means[,1]+1.96*reg_R_matrix_sds[,1], li=reg_R_matrix_means[,1]-1.96*reg_R_matrix_sds[,1],col="red",xlim=c(0.5,5),
       ylim=c(-0.2,1),xlab="Sample Size",ylab="",xaxt="n",pch=19,
       cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
plotCI(c(1, 2, 3, 4), reg_R_matrix_means[,2], ui=reg_R_matrix_means[,2]+1.96*reg_R_matrix_sds[,2], 
       li=reg_R_matrix_means[,2]-1.96*reg_R_matrix_sds[,2],col="blue",add=TRUE,pch=19,
       cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
plotCI(c(1.1, 2.1, 3.1, 4.1), reg_R_matrix_means[,3], ui=reg_R_matrix_means[,3]+1.96*reg_R_matrix_sds[,3],
       li=reg_R_matrix_means[,3]-1.96*reg_R_matrix_sds[,3],col="orange",add=TRUE,pch=19,
       cex.lab=1.5, cex.axis=1.1, cex.main=1.5)
abline(h = 0.269, lty = 2, lwd = 1.5)
title("Regression Adjustment Estimators")
legend("topright", pch = c(19, 19, 19), 
       col = c("red", "blue", "orange"), 
       legend = c("RA-R, IPW-A", "RA-R, RA-A", "RA-R, DR-A"),
       cex = 0.75)
axis(1, at = c(1,2,3,4), labels=sample_sizes, cex.axis=1.1)
dev.off()
