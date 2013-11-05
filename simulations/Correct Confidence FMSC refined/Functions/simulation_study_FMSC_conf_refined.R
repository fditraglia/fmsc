#October 7th, 2011

#This function examines the coverage of a refined confidence interval correction procedure based for FMSC moment selection based on a preliminary confidence interval for tau.

#-----------------------------------------------------------------------#
#CONSTANT MACROS - for this simulation, the true value is b.TRUE = 0.5   #
#-----------------------------------------------------------------------#

b.TRUE <- 0.5

#-----------------------------------------------------------------------#
#                           END OF MACROS                               #
#-----------------------------------------------------------------------#

#Debugging code
#N.reps <- 100
#sample.size <- 500
#relevance.w <- 0.4
#endog.w <- 0.2
  

#-----------------------------------------------------------------------#
#          FUNCTION: simulation.study.FMSC.conf.refined                #
#-----------------------------------------------------------------------#
  

#This function carries out the simulation study for fixed values of g and rho, storing relevant quantities

simulation.study.FMSC.conf.refined <- function(N.reps, sample.size, relevance.w, endog.w){
	
	#Initialize data.frame to hold results
	results <- matrix(NA, nrow = N.reps, ncol = 8)
	results <- as.data.frame(results)
	names(results) <- c('b', 'b.SE', 'FMSC', 'b1', 'b1.SE', 'FMSC.1', 'lower.refined', 'upper.refined')


  
  
	#Loop through replications of the simulation
	for(i in  1:N.reps){
		
    sim <- linear.simulation(relevance.w, endog.w, N = sample.size)
	
		temp.results <- compute.criteria.refined(X = sim$X, y = sim$y, Z1 = sim$Z1, Z2 = sim$Z2)
	
		results[i,] <- temp.results
	
	}#END for i
 
 attach(results)
 
 #Post-FMSC selection estimator
 FMSC.w <- FMSC < FMSC.1
 b.FMSC <- c(b[FMSC.w], b1[!FMSC.w])
 FMSC.mean <- mean(b.FMSC)

  #Coverge of corrected interval for mu.true that uses tau.hat for tau
  contains.truth <- (lower.refined <= b.TRUE) & (b.TRUE <= upper.refined) #elementwise AND

 refined.correct.cover <- sum(contains.truth)/length(contains.truth)
 
 detach(results)
 
 out <- data.frame(relevance.w, endog.w, refined.correct.cover)
 names(out) <- c('relevance.w', 'endog.w', 'refined.FMSC.correct.cover')
 
 return(out)
	
	
}#END simulation.study.FMSC.conf.refined