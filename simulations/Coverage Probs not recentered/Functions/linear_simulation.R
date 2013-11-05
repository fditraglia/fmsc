#This function generates simulations 

linear.simulation <- function(relevance.w, endog.w, N = 500, b0 = 0.5, p1 = 0.1, p2 = 0.1, p3 = 0.1){
	
	#Abbreviated names for parameters passed to this function	
	g <- relevance.w
	rho <- endog.w
	
	#Set the covariance between two error terms to maintain a constant covariance between x and the error in its regression equation of 0.5
	s <- 0.5 - g * rho
	
	#Covariance matrix of errors and "bad" instrument
	S <- matrix(c(1,   s, rho,
					s,   1,  0,
					rho, 0,  1), nrow = 3, ncol = 3)
	
	#Draw realizations for errors and "bad" instrument
	errors <- rmvnorm(N, mean = c(0,0,0), sigma = S, method = "chol")
	u <- errors[,1]
	e <- errors[,2]
	w <- errors[,3]
	
	#Draw realizations for "good" instruments
	z1 <- rnorm(N)
	z2 <- rnorm(N)
	z3 <- rnorm(N)
	
	#First stage
	x <- (p1 * z1) + (p2 * z2) + (p3 * z3) + (g * w) + e
	
	#Second stage
	y <- (b0 * x) + u
	
	out <- list(X = x, y = y, Z1 = cbind(z1, z2, z3), Z2 = w) 
	
	return(out)
	
	
}#END linear.simulation