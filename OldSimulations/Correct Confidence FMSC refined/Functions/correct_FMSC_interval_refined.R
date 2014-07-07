#October 7th, 2011

#This function corrects confidence intervals for FMSC moment selection by first constructing a confidence interval for tau.

compute.criteria.refined <- function(X, y, Z1, Z2, B = 1000, alpha = 0.05){
  
	
	#-----------------#
	#--PRELIMINARIES--#
	#-----------------#

  Z <- cbind(Z1, Z2)	#Full set of instruments
	
	#Capitalized version of ncol treats an array as a one column matrix
	q1 <- NCOL(Z1)	#Number of instruments in Z1
	q2 <- NCOL(Z2)	#Number of instruments in Z2
	
	q <- q1 + q2		#Number of instruments in Z
	p <- NCOL(X)		#Number of regression parameters
	p1 <- p			#Same number of parameters in valid model
	
	N <- length(y)	#Sample size
	
	
	#-------------------------------------------------------#
	#--TWO STAGE LEAST SQUARES CALCULATIONS FOR FULL MODEL--#
	#-------------------------------------------------------#
	
	#Fit 2SLS using full instrument set
	full.model <- reg2SLS(X, y, Z, first.stage = FALSE)
	b <- full.model$b
	C <- full.model$C
		
	#Estimate covariance matrix S of moment conditions for full model using centered, heteroscedasticity-robust estimator
	u <- as.vector(y - X %*% b)	#Residuals from full model
	D <- diag(u^2) #diag is picky about being passed a vector versus 1-col matrix!
	uu <- outer(u, u)
	S <- t(Z) %*% ( (D / N) - (uu / N^2) ) %*% Z
	
	#J-test Statistic for full model
	#Z.u <- t(Z) %*% u
	#J <- (1/N) * t(Z.u) %*% solve(S, Z.u)
	
	#Estimated Variance for b (full model estimate) based on S
	b.var <- (1/N) * C %*% S %*% t(C)
	b.SE <- sqrt(b.var)
	
	
	#--------------------------------------------------------#
	#--TWO STAGE LEAST SQUARES CALCULATIONS FOR VALID MODEL--#
	#--------------------------------------------------------#
	
	#Fit 2SLS using valid instrument set
	valid.model <- reg2SLS(X, y, Z1, first.stage = TRUE)
	b1 <- valid.model$b
	C1 <- valid.model$C
	X.star1 <- valid.model$X.star #Projection of X into colspace(Z1)
	
	#Estimate covariance matrix S1 of moment conditions for valid model without centering but allowing for heteroscedasticity
	u1 <- as.vector(y - X %*% b1)
	D1 <- diag(u1^2)
	S1 <- t(Z1) %*% ( D1 / N ) %*% Z1
	
	#J-test Statistic for valid model
	#Z1.u <- t(Z1) %*% u1
	#J1 <- (1/N) * t(Z1.u) %*% solve(S1, Z1.u)
	
	#Estimated Variance and std dev for b1 (valid model estimate) based on S1
	b1.var <- (1/N) * C1 %*% S1 %*% t(C1)
	b1.SE <- sqrt(b1.var)
	
	

	
	#---------------------------------------------#
	#--FOCUSED MOMENT SELECTION CRITERION (FMSC)--#
	#---------------------------------------------#

	#Construct Asymptotically Unbiased Estimate of outer product of bias vector
	tau <- t(Z2) %*% u1 / sqrt(N) 
	A <- (-1/N) * t(Z2) %*% X %*% C1
	Psi <- cbind(A, diag(nrow = nrow(A))) #Diag defaults to identity matrix
	Sigma <- Psi %*% S %*% t(Psi)
	T.hat <- as.matrix(outer(tau, tau)) - Sigma
	
	#In our calculation of the Asymptotic Bias, T.hat needs to be embedded as the bottom right block in a matrix of the same size as S that has zeros in all other locations. 
	Bias.matrix <- matrix(0, nrow = q, ncol = q)	
	Bias.matrix[(q1+1):q, (q1+1):q] <- T.hat
	
	#Full Model -- sum of estimated squared Asymptotic Bias and Asymptotic Variance
	AVAR <- N * b.var
	ABIAS.sq <- C %*% Bias.matrix %*% t(C)
	FMSC <- ABIAS.sq + AVAR
	
	#Valid Model -- by definition this is assumed to have no Asymptotic Bias so the FMSC is simply an estimate of the AVAR of sqrt(N)*(b1 - b.true)
	FMSC.1 <- N * b1.var
	

  #---------------------------------------------#
	#----REFINED CONFIDENCE INTERVAL FOR FMSC-----#
	#---------------------------------------------#
  
  #Construct 95% Confidence Interval for tau
  critical.value <- qchisq(p = 0.95, df = q2)
  half.width <- sqrt(critical.value * Sigma)
  tau.lower <- tau - half.width
  tau.upper <- tau + half.width
  
  #Now we create a function that returns the upper or lower end of a confidence interval for Lambda at a given value of tau. Eventually, we'll use this function to create a conservative confidence interval for Lambda by combining it with the 95% confidence interval for tau
  lambda.conf <- function(tau, tail, alpha = 0.05, B = 1000){
    #Simulate from distribution of M
    M.var <- S
    M.mean <- c(rep(0, q1), tau)
    
    #Generate Samples from M
    M.sims <- t(rmvnorm(B, M.mean, M.var)) #Transpose so that columns are realizations of the random vector
    
    #The simulations corresponding to the valid instruments only
    M1.sims <- M.sims[1:3,]
    
    #Calculate the limit version of T.hat for each sample
    Psi.M <- Psi %*% M.sims
    Psi.M.M.Psi <- Psi.M^2 #This next step only works if tau is a scalar so that Psi.M is a vector! (i.e. Psi %*% M_j is a scalar)
    T.hat.limit <- as.vector(Psi.M.M.Psi) - as.vector(Sigma)
    
    #Each realization of T.hat.limit yields a different value of FMSC.full.limit so we'll define a function to apply to the rows of T.hat.limit, treating it as a column vector
    
    get.FMSC.full.limit <- function(limit.bias){
      
      #Need to embed limit.bias as the bottom right block in a matrix of the same size as S that has zeros in all other locations. 
      bias.matrix <- matrix(0, nrow = q, ncol = q)
      bias.matrix[(q1+1):q, (q1+1):q] <- limit.bias
      
      bias.sq <- C %*% bias.matrix %*% t(C)
      out <- bias.sq + AVAR
      
      return(out)
      
      }#END get.FMSC.full.limit
    
    FMSC.full.limit <- sapply(T.hat.limit, get.FMSC.full.limit)
    choose.full <- as.vector(FMSC.full.limit) <= as.vector(FMSC.1)
    
    #Simulations from the limiting distribution of the target parameter, evaulated at the bias value tau
    Lambda.sim <- c(C %*% M.sims[,choose.full], C1 %*% M1.sims[,!choose.full])   
  
  #This is a confidence interval for sqrt(n)(mu.hat - mu.true)
    lower.lambda <- quantile(Lambda.sim, alpha/2)
    upper.lambda <- quantile(Lambda.sim, 1 - alpha/2)
    
    if(tail == 'lower'){return(lower.lambda)}
    if(tail == 'upper'){return(upper.lambda)}
    
  }#END lambda.conf

  
  lower.lambda <- optimize(f = lambda.conf, interval = c(tau.lower, tau.upper), tail = 'lower', maximum = FALSE)$objective
  
  upper.lambda <- optimize(f = lambda.conf, interval = c(tau.lower, tau.upper), tail = 'upper', maximum = TRUE)$objective
  


  #We want a confidence interval for mu.true
  FMSC.chooses.full <- FMSC <= FMSC.1
  mu.hat <- FMSC.chooses.full * b + (!FMSC.chooses.full) * b1 #post-FMSC estimator
  lower.refined<- mu.hat - (upper.lambda)/sqrt(N)
  upper.refined <- mu.hat - (lower.lambda)/sqrt(N)
  
  
	#---------------------------------------------#
	#--COLLECT AND RETURN CRITERIA AND ESTIMATES--#
	#---------------------------------------------#
	
  #Need to return the FMSC so that we can calculate the bias of the procedure which we'll need in the step where we see if the coverage is correct.
	out <- data.frame(b, b.SE, FMSC, b1, b1.SE, FMSC.1, lower.refined, upper.refined, row.names = NULL)
	
	return(out)
	
	
}#END compute.criteria