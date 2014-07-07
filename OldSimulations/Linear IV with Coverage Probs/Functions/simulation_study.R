#This function carries out the simulation study at a fixed pair of values for the endogeneity and relevance of w

#Updated September 23rd, 2011 to allow for varying sample sizes and calculate coverage of 95% confidence intervals around the quasi-true values. Also include selection via a downward J-test. Also streamlined the code to make it more readable.

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

#Function for checking coverage
coverage <- function(estimator, se, true.value, conf = 0.95){
  
  z.quantile <- qnorm(1 - ((1 - conf)/2))
  upper <- estimator + z.quantile * se
  lower <- estimator - z.quantile * se
  
  contains.truth <- (true.value <= upper) & (true.value >= lower) #Elementwise
  
  out <- sum(contains.truth)/length(contains.truth)
  return(out)
  
}#END coverage




#Convenience function to compute RMSE for this experiment where b = b.TRUE
rmse <- function(x){
	
	out <- sqrt((mean(x - b.TRUE))^2 + var(x))
	
	return(out)
	
}#END rmse

#-----------------------------------------------------------------------#
#                      END CONVENIENCE FUNCTIONS                        #
#-----------------------------------------------------------------------#
  
  
  
  
  

#-----------------------------------------------------------------------#
#                     FUNCTION: simulation.study                        #
#-----------------------------------------------------------------------#
  

#This function carries out the simulation study for fixed values of g and rho, storing relevant quantities

simulation.study <- function(N.reps, sample.size, relevance.w, endog.w){
	
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
	
	
  
  
	#Collect Quantities for Comparison of MSCs: true MSE, realized MSEs of each procedure, and correct decision rates.
  	
  	attach(results)
  	
  	#True RMSEs of each model
  	RMSE <- rmse(b) #Full Instrument Set
  	RMSE1 <- rmse(b1) #Valid Instrument Set
  	
    #Coverage of full and valid models, appropriately recentered
    COVERAGE <- coverage(estimator = b, se = b.SE, true.value = mean(b))
    COVERAGE1 <- coverage(estimator = b1, se = b1.SE, true.value = mean(b1))
  	
  	#True Minimum-MSE estimator - should we use w?
  	use.w <- (RMSE < RMSE1)
  
  
  
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
  	
  
  	#When does CCIC combined with MSC choose to include w?
  	
  	#CCIC chooses largest set of non-redundant instruments while MSC chooses largest set of valid instruments. In this example, since there are only two moment sets under consideration, the order in which the criteria are applied doesn't matter. The combined criterion chooses to include w precisely when both of the individual criteria choose to include w. Recall that a single ampersand means elementwise AND.
  	CC.MSC.BIC.w <- CC.BIC.w & GMM.BIC.w
  	CC.MSC.AIC.w <- CC.AIC.w & GMM.AIC.w
  	CC.MSC.HQ.w <- CC.HQ.w & GMM.HQ.w
  
  
  
  
  
  #Function to calculate correct decision rates of the various criteria.
  #If we really should be including w based on the finite-sample RMSE given the current values of relevance.w and endog.w, then the correct decision rate is simply the frequency with which w is included by a given criterion. Otherwise it is the frequency with which w is excluded.
	if(use.w){
    
    correct.rate <- function(criterion.includes.w){
      
      out <- sum(criterion.includes.w) / N.reps
      
    }#END correct.rate
    
	}else{
      
    correct.rate <- function(criterion.includes.w){
      
      out <- 1 - sum(criterion.includes.w) / N.reps
      
    }#END correct.rate
    
	}#END ifelse 
    
    
    
    #Now, set up a function that collects all the desired information from a given criterion's selections. This includes bias, RMSE, correct decision rate, and the actual coverage achieved by a naive procedure that ignores moment selection uncertainty when constructing a 95% interval. For the last item, we need to recenter to remove any bias that is present. This is so that we can isolate coverage problems arising from moment selection (and the non-normality that results) from the bias (and hence incorrect centering of the confidence interval) that results from including invalid moment conditions.
    summarize.criterion <- function(criterion.includes.w){
      
      #Post-selection estimator
      b.criterion <- c(b[criterion.includes.w], b1[!criterion.includes.w])
      
      #Bias of post-selection estimator
      BIAS <- mean(b.criterion - b.TRUE)
      
      #RMSE of post-selection estimator
      RMSE <- rmse(b.criterion)
      
      #Post-selection naive standard error
      SE.criterion <- c(b.SE[criterion.includes.w], b1.SE[!criterion.includes.w])
      
      #Coverage for naive confidence interval where we take the truth to be (biased) mean of the finite sample distribution of the post-selection estimator
      COVERAGE <- coverage(estimator = b.criterion, se = SE.criterion, true.value = mean(b.criterion))
      
      #Correct Decision Rate
      CORRECT <- correct.rate(criterion.includes.w)
      
      return(data.frame(BIAS, RMSE, COVERAGE, CORRECT))
      
    }#END summarize.criterion
    
    
    
    #Apply the function summarize.criterion to each of the selection procedures:
    FMSC.summary <- summarize.criterion(FMSC.w)
		GMM.BIC.summary <- summarize.criterion(GMM.BIC.w)
		GMM.AIC.summary <- summarize.criterion(GMM.AIC.w)
		GMM.HQ.summary <- summarize.criterion(GMM.HQ.w)
		CC.BIC.summary <- summarize.criterion(CC.BIC.w)
		CC.AIC.summary <- summarize.criterion(CC.AIC.w)
		CC.HQ.summary <- summarize.criterion(CC.HQ.w)
		CC.MSC.BIC.summary <- summarize.criterion(CC.MSC.BIC.w)
		CC.MSC.AIC.summary <- summarize.criterion(CC.MSC.AIC.w)
		CC.MSC.HQ.summary <- summarize.criterion(CC.MSC.HQ.w)
    J.90.summary <- summarize.criterion(J.90.w)
    J.95.summary <- summarize.criterion(J.95.w)
  

	  detach(results)
  
    
  #Names of criteria
  criterion.names <- c('FMSC', 'GMM.BIC', 'GMM.AIC', 'GMM.HQ', 'CC.BIC', 'CC.AIC', 'CC.HQ', 'CC.MSC.BIC', 'CC.MSC.AIC', 'CC.MSC.HQ', 'J.90', 'J.95')
    
  #Names of summary statistics
  summary.names <- c('BIAS', 'RMSE', 'COVERAGE', 'CORRECT')
    
  #Names for output 
  out.names <- paste(rep(criterion.names, each = length(summary.names)), rep(summary.names, length(criterion.names)), sep = '.')
   
	out <- data.frame(	relevance.w, 
						endog.w, 
						RMSE, 
						RMSE1,
            COVERAGE,
            COVERAGE1,
						FMSC.summary, 
						GMM.BIC.summary, 
						GMM.AIC.summary, 
						GMM.HQ.summary, 
						CC.BIC.summary, 
						CC.AIC.summary, 
						CC.HQ.summary, 
						CC.MSC.BIC.summary, 
						CC.MSC.AIC.summary, 
						CC.MSC.HQ.summary,
            J.90.summary,
            J.95.summary)

  names(out) <- c('relevance.w', 'endog.w', 'RMSE', 'RMSE1', 'COVERAGE', 'COVERAGE1', out.names)
	
	return(out)
	
	
}#END simulation.study