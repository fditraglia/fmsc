#October 5th, 2011

#This function examines the coverage of a naive confidence interval correction procedure based for FMSC moment selection that treats tau.hat as the true value of tau

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
#            FUNCTION: simulation.study.FMSC.conf.naive                 #
#-----------------------------------------------------------------------#
  

#This function carries out the simulation study for fixed values of g and rho, storing relevant quantities

simulation.study.FMSC.conf.naive <- function(N.reps, sample.size, relevance.w, endog.w){
	
	#Initialize data.frame to hold results
	results <- matrix(NA, nrow = N.reps, ncol = 8)
	results <- as.data.frame(results)
	names(results) <- c('b', 'b.SE', 'FMSC', 'b1', 'b1.SE', 'FMSC.1', 'lower.naive', 'upper.naive')


  
  
	#Loop through replications of the simulation
	for(i in  1:N.reps){
	
		sim <- linear.simulation(relevance.w, endog.w, N = sample.size)
	
		temp.results <- compute.criteria.naive(X = sim$X, y = sim$y, Z1 = sim$Z1, Z2 = sim$Z2)
	
		results[i,] <- temp.results
	
	}#END for i
 
 attach(results)
 
 #Post-FMSC selection estimator
 FMSC.w <- FMSC < FMSC.1
 b.FMSC <- c(b[FMSC.w], b1[!FMSC.w])
 FMSC.mean <- mean(b.FMSC)

  #Coverge of corrected interval for mu.true that uses tau.hat for tau
  contains.truth <- (lower.naive <= b.TRUE) & (b.TRUE <= upper.naive) #elementwise AND

 naive.correct.cover <- sum(contains.truth)/length(contains.truth)
 
 detach(results)
 
 out <- data.frame(relevance.w, endog.w, naive.correct.cover)
 names(out) <- c('relevance.w', 'endog.w', 'naive.FMSC.correct.cover')
 
 return(out)
	
	
}#END simulation.study.FMSC.conf.naive