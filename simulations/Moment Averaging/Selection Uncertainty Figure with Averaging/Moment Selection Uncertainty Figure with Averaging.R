#September 27th, 2011

#IMPORTANT: The FMSC could be negative since our estimator of the squared bias could be negative. This isn't necessarily a problem for moment selection, but it could make the weights close to 0-1 for moment averaging. Here I also consider a modification in which we replace negative values of the FIC with N * b1.var, i.e. setting the bias parameter equal to zero.

#This script compares the distributions of post-selection and post-averaging estimators.

sample.size <- 500
N.reps <- 10000
endog.w <- 0.2
relevance.w <- 0.4
b.true <- 0.5
output.dir <- './Selection Uncertainty Figure with Averaging'


#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Moment Averaging"

setwd(root.dir)
setwd("./Functions")

source("reg2SLS.R")
source("compute_criteria.R")
source("linear_simulation.R")
library(mvtnorm)



#Convenience function to compute RMSE for this experiment where b = b.TRUE
rmse <- function(x){
  
	out <- sqrt((mean(x - b.TRUE))^2 + var(x))
	
	return(out)
	
}#END rmse



#Convenience function for smoothed moment averaging according to the formula from Buckland et al (1997)
weights <- function(criteria, k = 1){
  
  #For numerical stability, first subtract the minimum value of the criterion
  min.criteria <- min(criteria)
  numerators <- exp(-0.5 * k * (criteria - min.criteria))
  denominator <- sum(numerators)
  
  moment.weights <- numerators / denominator
  
  return(moment.weights)
  
}#END weights






#Initialize data.frame to hold results
results <- matrix(NA, nrow = N.reps, ncol = 20)
results <- as.data.frame(results)
names(results) <- c('b', 'b.SE', 'J', 'GMM.AIC', 'GMM.BIC', 'GMM.HQ', 'CC.AIC', 'CC.BIC', 'CC.HQ', 'FMSC', 'b1', 'b1.SE', 'J1', 'GMM.AIC.1', 'GMM.BIC.1', 'GMM.HQ.1', 'CC.AIC.1', 'CC.BIC.1', 'CC.HQ.1', 'FMSC.1')


#Set Seed for replication
set.seed(1100)

  
#Loop through replications of the simulation
for(i in  1:N.reps){
	
	sim <- linear.simulation(relevance.w, endog.w, N = sample.size)
	
	temp.results <- compute.criteria(X = sim$X, y = sim$y, Z1 = sim$Z1, Z2 = sim$Z2)
	
	results[i,] <- temp.results
	
}#END for i
	

#Construct the post-selection estimators
  	
attach(results)
  	
  
#When do the individual criteria choose to include w?
FMSC.w <- FMSC < FMSC.1
GMM.BIC.w <- GMM.BIC < GMM.BIC.1
GMM.AIC.w <- GMM.AIC < GMM.AIC.1
GMM.HQ.w <- GMM.HQ < GMM.HQ.1
CC.BIC.w <- CC.BIC < CC.BIC.1
CC.AIC.w <- CC.AIC < CC.AIC.1
CC.HQ.w <- CC.HQ < CC.HQ.1
J.90.w <- J <= qchisq(0.90, 3) #Downward J-test at 90%
J.95.w <- J <= qchisq(0.95, 3) #Downward J-test at 95%
CC.MSC.BIC.w <- CC.BIC.w & GMM.BIC.w 
CC.MSC.AIC.w <- CC.AIC.w & GMM.AIC.w
CC.MSC.HQ.w <- CC.HQ.w & GMM.HQ.w

#Function to construct post-selection estimator given a selection vector that contains 1 when the criterion chooses to include w, 0 otherwise  
post.selection <- function(selection.vector){
  
  estimator <- selection.vector * b + (1 - selection.vector) * b1
  
}#END post.selection
  

#Post-selection estimators 
b.full <- post.selection(rep(1, N.reps))  #Always choose w
b.valid <- post.selection(rep(0, N.reps)) #Never choose w
b.FMSC <- post.selection(FMSC.w)
b.BIC <- post.selection(GMM.BIC.w)
b.AIC <- post.selection(GMM.AIC.w)
b.HQ <- post.selection(GMM.HQ.w)
b.CC.BIC <- post.selection(CC.BIC.w)
b.CC.AIC <- post.selection(CC.AIC.w)
b.CC.HQ <- post.selection(CC.HQ.w)
b.J.90 <- post.selection(J.90.w)
b.J.95 <- post.selection(J.95.w)
b.CC.MSC.BIC <- post.selection(CC.MSC.BIC.w)
b.CC.MSC.AIC <- post.selection(CC.MSC.AIC.w)
b.CC.MSC.HQ <- post.selection(CC.MSC.HQ.w)


#Modified FMSC - If FMSC < 0, set it equal to N * b1.var
FMSC.star <- FMSC * (FMSC > 0) + (sample.size * (b1.SE)^2) * !(FMSC > 0)
FMSC.star.w <- FMSC.star <= FMSC.1
b.FMSC.star <- post.selection(FMSC.star.w)

#Moment AVeraging
#Weights for averaging
BIC.weights <- t(apply(cbind(GMM.BIC, GMM.BIC.1), 1, weights))
AIC.weights <- t(apply(cbind(GMM.AIC, GMM.AIC.1), 1, weights))
HQ.weights <- t(apply(cbind(GMM.HQ, GMM.HQ.1), 1, weights))
FMSC.weights <-t(apply(cbind(FMSC, FMSC.1), 1, weights, k = 0.01))
FMSC.star.weights <-t(apply(cbind(FMSC.star, FMSC.1), 1, weights))  

#Moment Average Estimators
b.BIC.avg <- apply(BIC.weights * cbind(b, b1), 1, sum)
b.AIC.avg <- apply(AIC.weights * cbind(b, b1), 1, sum)
b.HQ.avg <- apply(HQ.weights * cbind(b, b1), 1, sum)
b.FMSC.avg <- apply(FMSC.weights * cbind(b, b1), 1, sum)
b.FMSC.star.avg <- apply(FMSC.star.weights * cbind(b, b1), 1, sum)

  
detach(results)

#Function to plot the density of a given post-selection estimator alongside that of b.full and b.valid

dist.plot <- function(b.estimator, name){
  
  f.valid <- density(b.valid)
  f.full <- density(b.full)
  f.criterion <- density(b.estimator)
  
  xmin <- min(min(f.valid$x), min(f.full$x), min(f.criterion$x))
  xmax <- max(max(f.valid$x), max(f.full$x), max(f.criterion$x))
  ymin <- min(min(f.valid$y), min(f.full$y), min(f.criterion$y))
  ymax <- max(max(f.valid$y), max(f.full$y), max(f.criterion$y))
  
  plot(f.valid$x, f.valid$y, xlim = c(-0.5, xmax), ylim = c(ymin, ymax), xlab = '', ylab = 'Density', type = 'l', lty = 2)
  
  points(f.full$x, f.full$y, type = 'l', lty = 3)

  points(f.criterion$x, f.criterion$y, type = 'l', lty = 1)

  legend('topleft', c('Valid', 'Full', name), lty = c(2,3,1))
  
}#END dist.plot




dist.plot(b.FMSC, 'FMSC')
dist.plot(b.FMSC.avg, 'FMSC.avg')

dist.plot(b.BIC, 'BIC')
dist.plot(b.BIC.avg, 'BIC.avg')


dist.plot(b.BIC, 'AIC')
dist.plot(b.BIC.avg, 'AIC.avg')

dist.plot(b.FMSC.star, 'FMSC*')
dist.plot(b.FMSC.star.avg, 'FMSC* avg')

#Make plots and write to file
#setwd(root.dir)
#setwd(output.dir)

#pdf('GMM_BIC.pdf')
#dist.plot(b.BIC, 'GMM-BIC')
#dev.off()

#pdf('GMM_AIC.pdf')
#dist.plot(b.AIC, 'GMM-AIC')
#dev.off()

#pdf('GMM_HQ.pdf')
#dist.plot(b.HQ, 'GMM-HQ')
#dev.off()

#pdf('FMSC.pdf')
#dist.plot(b.FMSC, 'FMSC')
#dev.off()

#pdf('J_90.pdf')
#dist.plot(b.J.90, 'J-test (90%)')
#dev.off()

#pdf('J_95.pdf')
#dist.plot(b.J.95, 'J-test (95%)')
#dev.off()
