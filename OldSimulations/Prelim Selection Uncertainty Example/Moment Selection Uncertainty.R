#Simple simulation to show how choosing moment conditions based on J-statistic (and indeed the other criteria) the standard errors of the final estimate

#To start, try the situation in which the relevance of w is 0.5 and the endogeneity is 0.2. In this situation based on RMSE we should exclude w.

sample.size <- 500
N.reps <- 10000
endog <- 0.2
rel <- 0.5


root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Selection Uncertainty Example"

setwd(root.dir)
setwd('./Functions')

library(mvtnorm)
library(sem)
source('linear_simulation.R')
source('reg2SLS.R')
source('compute_criteria.R')

#Initialize data.frame to hold results
results <- matrix(NA, nrow = N.reps, ncol = 20)
results <- as.data.frame(results)
names(results) <- c('b', 'b.var', 'J', 'GMM.AIC', 'GMM.BIC', 'GMM.HQ', 'CC.AIC', 'CC.BIC', 'CC.HQ', 'FMSC', 'b1', 'b1.var', 'J1', 'GMM.AIC.1', 'GMM.BIC.1', 'GMM.HQ.1', 'CC.AIC.1', 'CC.BIC.1', 'CC.HQ.1', 'FMSC.1')

#Makes sure that the standard errors I'm computing are right
#test.results <- matrix(NA, nrow = N.reps, ncol = 2)
#test.results <- as.data.frame(test.results)
#names(test.results) <- c('estimate', 'SE')


set.seed(1647)


#Loop through replications of the simulation
for(i in  1:N.reps){
	
  #Generate another batch of simulations from the linear model
  sim <- linear.simulation(relevance.w = rel, endog.w = endog, N = sample.size)
	
  #Calculate estimator, std error, and selection criteria
  temp.results <- compute.criteria(X = sim$X, y = sim$y, Z1 = sim$Z1, Z2 = sim$Z2)
	
  #Store the results
	results[i,] <- temp.results
  
  #Check that this is working by using the canned 2SLS package
  #y <- sim$y
  #x <- sim$X
  #z1 <- sim$Z1[,1]
  #z2 <- sim$Z1[,2]
  #z3 <- sim$Z1[,3]
  #canned.2SLS <- tsls(y ~ x - 1, ~ z1 + z2 + z3 - 1)
	
}#END for i

  
#Now we can look at the properties of various estimators
attach(results)

#My code calculates estimate of sampling variance of b rather than std. dev.
b.SE <- sqrt(b.var)
b1.SE <- sqrt(b1.var)

#When is the full model chosen by the downward J-test at 90% significance?
J.chooses.full <- J < qchisq(0.90, df = 3)
  
#When is the full model chosen by the GMM.AIC? (remember, we minimize this)
AIC.chooses.full <- GMM.AIC <= GMM.AIC.1
  
#BIC? (Note that this won't have an asymptotic distribution, but we can still look at the small sample behavior)
BIC.chooses.full <- GMM.BIC <= GMM.BIC.1  

#When is the full model chosen by the FMSC?
FMSC.chooses.full <- FMSC <= FMSC.1



#Function to construct post-selection   
post.selection <- function(selection.vector){
  
  estimator <- selection.vector * b + (1 - selection.vector) * b1
  SE <- selection.vector * b.SE + (1 - selection.vector) * b1.SE
  
  out <- data.frame(estimator = estimator, SE = SE)
  
}#END post.selection

#Post-selection estimators and standard errors ignoring moment selection
b.AIC <- post.selection(AIC.chooses.full)
b.BIC <- post.selection(BIC.chooses.full)
b.J <- post.selection(J.chooses.full)
b.FMSC <- post.selection(FMSC.chooses.full)
b.full <- post.selection(rep(1, length(b)))
b.valid <- post.selection(rep(0, length(b)))
  



#This isn't the right route to take since the coverage is impacted by the bias term. The key point here is that standard errors are TOO SMALL!
sd(b.valid$estimator) - mean(b.valid$SE) 
sd(b.full$estimator) - mean(b.full$SE) 
sd(b.BIC$estimator) - mean(b.BIC$SE) 
sd(b.AIC$estimator) - mean(b.AIC$SE) 
sd(b.FMSC$estimator) - mean(b.FMSC$SE) 
#Well, this isn't quite right either. The post-selection distributions are not only more variable but also highly non-normal. So perhaps we need to look at coverage after all...



#For coverage what we should do is find the mean of the sampling distribution of the estimator and subtract this. That is, subtract the same thing for each estimate NOT a different thing for each estimate, as this basically removes the source of variation that we're interested in trying to account for: the moment selection uncertainty.


#Function for checking coverage
coverage <- function(b.object, true.value, conf = 0.95){
  
  estimator <- b.object$estimator
  se <- b.object$SE
  
  z.quantile <- qnorm(1 - ((1 - conf)/2))
  upper <- estimator + z.quantile * se
  lower <- estimator - z.quantile * se
  
  contains.truth <- (true.value <= upper) & (true.value >= lower)
  
  out <- sum(contains.truth)/length(contains.truth)
  return(out)
  
}#END coverage


coverage(b.valid, true.value = mean(b.valid$estimator))
coverage(b.full, true.value = mean(b.full$estimator))
coverage(b.BIC, true.value = mean(b.BIC$estimator))
coverage(b.AIC, true.value = mean(b.AIC$estimator))
coverage(b.FMSC, true.value = mean(b.FMSC$estimator))
coverage(b.J, true.value = mean(b.J$estimator))




library(MASS)
truehist(b.BIC$estimator,100)
truehist(b.valid$estimator,100)
truehist(b.full$estimator,100)
truehist(b.FMSC$estimator,100)
truehist(b.AIC$estimator,100)
truehist(b.J$estimator,100)

plot(density(b.full$estimator))
plot(density(b.valid$estimator))
plot(density(b.BIC$estimator))
plot(density(b.AIC$estimator))
plot(density(b.FMSC$estimator))
plot(density(b.J$estimator))


#FMSC
f.valid <- density(b.valid$estimator)
f.full <- density(b.full$estimator)
f.fmsc <- density(b.FMSC$estimator)
xmin <- min(min(f.valid$x), min(f.full$x), min(f.fmsc$x))
xmax <- max(max(f.valid$x), max(f.full$x), max(f.fmsc$x))
ymin <- min(min(f.valid$y), min(f.full$y), min(f.fmsc$y))
ymax <- max(max(f.valid$y), max(f.full$y), max(f.fmsc$y))
plot(f.valid$x, f.valid$y, xlim = c(-0.5, xmax), ylim = c(ymin, ymax), xlab = expression(beta), ylab = 'Density', col = 'red', type = 'l')
points(f.full$x, f.full$y, type = 'l', col = 'blue')
points(f.fmsc$x, f.fmsc$y, type = 'l', col = 'black')
legend('topleft', c('Valid', 'Full', 'FMSC'), col = c('red', 'blue', 'black'), lty = 1)



#GMM-BIC
f.bic <- density(b.BIC$estimator)
xmin <- min(min(f.valid$x), min(f.full$x), min(f.bic$x))
xmax <- max(max(f.valid$x), max(f.full$x), max(f.bic$x))
ymin <- min(min(f.valid$y), min(f.full$y), min(f.bic$y))
ymax <- max(max(f.valid$y), max(f.full$y), max(f.bic$y))
plot(f.valid$x, f.valid$y, xlim = c(-0.5, xmax), ylim = c(ymin, ymax), xlab = expression(beta), ylab = 'Density', col = 'red', type = 'l')
points(f.full$x, f.full$y, type = 'l', col = 'blue')
points(f.bic$x, f.bic$y, type = 'l', col = 'black')
legend('topleft', c('Valid', 'Full', 'BIC'), col = c('red', 'blue', 'black'), lty = 1)


#GMM-AIC
f.aic <- density(b.AIC$estimator)
xmin <- min(min(f.valid$x), min(f.full$x), min(f.aic$x))
xmax <- max(max(f.valid$x), max(f.full$x), max(f.aic$x))
ymin <- min(min(f.valid$y), min(f.full$y), min(f.aic$y))
ymax <- max(max(f.valid$y), max(f.full$y), max(f.aic$y))
plot(f.valid$x, f.valid$y, xlim = c(-0.5, xmax), ylim = c(ymin, ymax), xlab = expression(beta), ylab = 'Density', col = 'red', type = 'l', lty = 2)
points(f.full$x, f.full$y, type = 'l', col = 'blue', lty = 3)
points(f.aic$x, f.aic$y, type = 'l', col = 'black', lty = 1)
legend('topleft', c('Valid', 'Full', 'AIC'), col = c('red', 'blue', 'black'), lty = c(2,3,1))





qqnorm(b.valid$estimator)
qqline(b.valid$estimator)

qqnorm(b.full$estimator)
qqline(b.full$estimator)

qqnorm(b.FMSC$estimator)
qqline(b.FMSC$estimator)

qqnorm(b.BIC$estimator)
qqline(b.BIC$estimator)

qqnorm(b.AIC$estimator)
qqline(b.AIC$estimator)


#Convenience function to compute RMSE when true value is 0.5
rmse <- function(x){return(sqrt((mean(x - 0.5))^2 + var(x)))}

rmse(b.valid$estimator)
rmse(b.full$estimator)
rmse(b.FMSC$estimator)
rmse(b.BIC$estimator)
rmse(b.AIC$estimator)
rmse(b.J$estimator)

#Same for median absolute error
#mae <- function(x){return(mean(abs(x - 0.5)))}

#mae(b.valid$estimator)
#mae(b.full$estimator)
#mae(b.FMSC$estimator)
#mae(b.BIC$estimator)
#mae(b.AIC$estimator)
#mae(b.J$estimator)



#Compare empirical CDFs?
plot(ecdf(b.AIC$estimator))
plot(ecdf(b.BIC$estimator))
plot(ecdf(b.J$estimator))

plot(ecdf(b.FMSC$estimator))
test.points <- seq(from = -3, to = 2, by = 0.01)
Fnorm.FMSC <- pnorm(test.points, mean = mean(b.FMSC$estimator), sd = mean( b.FMSC$SE))
points(test.points, Fnorm.FMSC, type = 'l', lty = 2)




