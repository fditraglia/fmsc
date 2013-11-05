#September 27th, 2011


#This function carries out the simulation study at a fixed pair of values for the endogeneity and relevance of w to evaluate the performance of smoothed GMM-AIC, BIC, and HQ moment averaging as well as FMSC moment averaging. We only compute realized RMSE since naive coverage and correct decision rates wouldn't make sense in this context.


#-----------------------------------------------------------------------#
#CONSTANT MACROS - for this simulation, the true value is b.TRUE = 0.5   #
#-----------------------------------------------------------------------#

b.TRUE <- 0.5

#-----------------------------------------------------------------------#
#                           END OF MACROS                               #
#-----------------------------------------------------------------------#




#-----------------------------------------------------------------------#
#                         CONVENIENCE FUNCTIONS                         #
#-----------------------------------------------------------------------#



#Convenience function to compute RMSE for this experiment where b = b.TRUE
rmse <- function(x){
	
	out <- sqrt((mean(x - b.TRUE))^2 + var(x))
	
	return(out)
	
}#END rmse



#Convenience function for smoothed moment averaging according to the formula from Buckland et al (1997). The parameter k allows us to adjust the weighting. When k is small, the weighting will be more uniform.
weights <- function(criteria, k = 1){
  
  #For numerical stability, first subtract the minimum value of the criterion
  min.criteria <- min(criteria)
  numerators <- exp(-0.5 * k * (criteria - min.criteria))
  denominator <- sum(numerators)
  
  moment.weights <- numerators / denominator
  
  return(moment.weights)
  
}#END weights




#Test this out applied to a matrix:
#x <- c(2,3,5)
#weights(x)
#M <- matrix(rep(c(2,3,5), each = 3), nrow = 3, ncol = 3)
#t(apply(M, 1, weights))

#-----------------------------------------------------------------------#
#                      END CONVENIENCE FUNCTIONS                        #
#-----------------------------------------------------------------------#
  
  
  
  
  

#-----------------------------------------------------------------------#
#                  FUNCTION: simulation.study.avg                       #
#-----------------------------------------------------------------------#
  

#This function carries out the simulation study for fixed values of g and rho, storing relevant quantities

simulation.study.avg <- function(N.reps, sample.size, relevance.w, endog.w){
	
	#Initialize data.frame to hold results
	results <- matrix(NA, nrow = N.reps, ncol = 20)
	results <- as.data.frame(results)
	names(results) <- c('b', 'b.SE', 'J', 'GMM.AIC', 'GMM.BIC', 'GMM.HQ', 'CC.AIC', 'CC.BIC', 'CC.HQ', 'FMSC', 'b1', 'b1.SE', 'J1', 'GMM.AIC.1', 'GMM.BIC.1', 'GMM.HQ.1', 'CC.AIC.1', 'CC.BIC.1', 'CC.HQ.1', 'FMSC.1')


  
  
	#Loop through replications of the simulation
	for(i in  1:N.reps){
	
		sim <- linear.simulation(relevance.w, endog.w, N = sample.size)
	
		temp.results <- compute.criteria(X = sim$X, y = sim$y, Z1 = sim$Z1, Z2 = sim$Z2)
	
		results[i,] <- temp.results
	
	}#END for i
	
	
  
  
	#Construct Moment Averaged Estimators
  	
  	attach(results)
  	
    #Weights for averaging
    BIC.weights <- t(apply(cbind(GMM.BIC, GMM.BIC.1), 1, weights))
    AIC.weights <- t(apply(cbind(GMM.AIC, GMM.AIC.1), 1, weights))
    HQ.weights <- t(apply(cbind(GMM.HQ, GMM.HQ.1), 1, weights))
    FMSC.weights <- t(apply(cbind(FMSC, FMSC.1), 1, weights))
    FMSC.weights.k <- t(apply(cbind(FMSC, FMSC.1), 1, weights, k = 0.01))
  
    #Moment Average Estimators
    b.BIC.avg <- apply(BIC.weights * cbind(b, b1), 1, sum)
    b.AIC.avg <- apply(AIC.weights * cbind(b, b1), 1, sum)
    b.HQ.avg <- apply(HQ.weights * cbind(b, b1), 1, sum)
    b.FMSC.avg <- apply(FMSC.weights * cbind(b, b1), 1, sum)
    b.FMSC.avg.k <- apply(FMSC.weights.k * cbind(b, b1), 1, sum)
    
    #RMSE of post-selection estimator
    RMSE.BIC.avg <- rmse(b.BIC.avg)
    RMSE.AIC.avg <- rmse(b.AIC.avg)
    RMSE.HQ.avg <- rmse(b.HQ.avg)
    RMSE.FMSC.avg <- rmse(b.FMSC.avg)
    RMSE.FMSC.avg.k <- rmse(b.FMSC.avg.k)

	  detach(results)
  
   
	out <- data.frame(relevance.w, endog.w, RMSE.BIC.avg, RMSE.AIC.avg, RMSE.HQ.avg, RMSE.FMSC.avg, RMSE.FMSC.avg.k)

	
	return(out)
	
	
}#END simulation.study.avg