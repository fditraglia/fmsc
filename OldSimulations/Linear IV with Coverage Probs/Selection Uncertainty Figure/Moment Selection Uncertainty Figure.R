#September 24th, 2011

#This script makes a figure illustrating the non-normality of the distribution of estimators post-selection

sample.size <- 500
N.reps <- 10000
endog.w <- 0.2
relevance.w <- 0.4

output.dir <- './Selection Uncertainty Figure'


#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Linear IV with Coverage Probs"

setwd(root.dir)
setwd("./Functions")

source("reg2SLS.R")
source("compute_criteria.R")
source("linear_simulation.R")
library(mvtnorm)


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



#Make plots and write to file
setwd(root.dir)
setwd(output.dir)

pdf('GMM_BIC.pdf')
dist.plot(b.BIC, 'GMM-BIC')
dev.off()

pdf('GMM_AIC.pdf')
dist.plot(b.AIC, 'GMM-AIC')
dev.off()

pdf('GMM_HQ.pdf')
dist.plot(b.HQ, 'GMM-HQ')
dev.off()

pdf('FMSC.pdf')
dist.plot(b.FMSC, 'FMSC')
dev.off()

pdf('J_90.pdf')
dist.plot(b.J.90, 'J-test (90%)')
dev.off()

pdf('J_95.pdf')
dist.plot(b.J.95, 'J-test (95%)')
dev.off()
