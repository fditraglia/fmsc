#October 12th, 2011

#This script implements FMSC moment selection and confidence interval correction for Linear IV models 


#------------------------------------------------------#
#-------------------FUNCTION: reg2SLS------------------#
#------------------------------------------------------#
#This function returns the 2SLS regression coefficient b for a regression of y on X with instruments Z, as well as the product needed for estimating the variance matrix. If first.stage is TRUE, it also returns the ``transformed'' X variables that result from the first-stage regression.
reg2SLS <- function(X, y, Z, first.stage = FALSE){
  
	#Calculate QR Decomposition of Z
	QR.z <- qr(Z)
	Qz <- qr.Q(QR.z)
	Rz <- qr.R(QR.z)
	
	#Calculate (Z'Z)^{-1}
	Rz.inv <- backsolve(Rz, diag(1, nrow = nrow(Rz), ncol = ncol(Rz)))
	ZZ.inv <- Rz.inv %*% t(Rz.inv)
	
	#Transform X and y
	X.star <- t(Qz) %*% X
	y.star <- t(Qz) %*% y
	
	#Calculate QR decomposition of transformed X
	QR<- qr(X.star)
	Q <- qr.Q(QR)
	R <- qr.R(QR)
	
	#Calculate (X'X)^{-1} for the transformed X
	R.inv <- backsolve(R, diag(1, nrow = nrow(R), ncol = ncol(R)))
	XX.star.inv <- R.inv %*% t(R.inv)
	
	#Calculate C, the matrix used to estimate the variance via CSC'
	N <- length(y)
	C <- N * XX.star.inv %*% t(X) %*% Z %*% ZZ.inv
	
	#Calculate 2SLS coefficient
	c <- t(Q) %*% y.star
	b <- backsolve(R, c)
	
	if(first.stage == TRUE){
		
		out <- list(b = b, C = C, X.star = X.star)
		
	}else{
		
		out <- list(b = b, C = C)
		
	}#END if else (first.stage == TRUE)
	
	return(out)
	
}#END reg2SLS
 
 
#------------------------------------------------------#
#----------------END FUNCTION: reg2SLS-----------------#
#------------------------------------------------------#


FMSC.2SLS <- function(X, y, Z.valid, Z.additional, instrument.blocks, nabla.mu){
  
  Z <- cbind(Z.valid, Z.additional)
  Z1 <- Z.valid
  Z2 <- Z.additional
  
  #Capitalized version of ncol treats an array as a one column matrix
	q1 <- NCOL(Z1)	#Number of instruments in Z1
	q2 <- NCOL(Z2)	#Number of instruments in Z2
	p <- NCOL(X)		#Number of regression parameters
	N <- length(y)	#Sample size
	
	#Fit 2SLS using full instrument set
	full.model <- reg2SLS(X, y, Z)
	b <- full.model$b
	C <- full.model$C
		
	#Estimate covariance matrix S of moment conditions for full model using centered, heteroscedasticity-robust estimator
	u <- as.vector(y - X %*% b)	#Residuals from full model
	D <- diag(u^2) #diag is picky about being passed a vector versus 1-col matrix!
	uu <- outer(u, u)
	Omega <- t(Z) %*% ( (D / N) - (uu / N^2) ) %*% Z
	
	#Estimated Var and SE for b (full model estimate) based on Omega
	b.var <- (1/N) * C %*% Omega %*% t(C)
	
	#Fit 2SLS using valid instrument set
	valid.model <- reg2SLS(X, y, Z1)
	b1 <- valid.model$b
	C1 <- valid.model$C
	
	#Estimate covariance matrix S1 of moment conditions for valid model without centering but allowing for heteroscedasticity
	u1 <- as.vector(y - X %*% b1)
	D1 <- diag(u1^2)
	Omega1 <- t(Z1) %*% ( D1 / N ) %*% Z1
	
	#Estimated Var and SE for b1 (valid model estimate) based on S1
	b1.var <- (1/N) * C1 %*% Omega1 %*% t(C1)
	
	#Construct Asymptotically Unbiased Estimate of outer product of bias vector
	tau <- t(Z2) %*% u1 / sqrt(N) 
	A <- (-1/N) * t(Z2) %*% X %*% C1
	Psi <- cbind(A, diag(nrow = nrow(A))) #Diag defaults to identity matrix
	Sigma <- Psi %*% Omega %*% t(Psi) #Variance matrix of tau.hat
	T.hat <- tau %*% t(tau) - Sigma
	
	#In our calculation of the Asymptotic Bias, T.hat needs to be embedded as the bottom right block in a matrix of the same size as S that has zeros in all other locations. 
	Bias.matrix <- matrix(0, nrow = (q1+q2), ncol = (q1+q2))	
	Bias.matrix[(q1+1):(q1+q2), (q1+1):(q1+q2)] <- T.hat
	
  #Evaluate the derivative of mu at the valid model estimator
  nabla.mu.valid <- nabla.mu(b1)
  #nabla.mu.valid <- c(0,0,1) #Coeff on malfal is the target parameter
  
  #The inner matrix used to calculate the FMSC for all of the potentially mis-specified models
  M.inner <- Bias.matrix + Omega
  
  #Calculate the FMSC for the various submodels
  
  #Initialize vector to store the submodel FMSC values and list to store corresponding K.S matrices for confidence interval correction.
  if(is.matrix(instrument.blocks)){N.blocks <- nrow(instrument.blocks)}
  if(is.vector(instrument.blocks)){N.blocks <- 1}
  
  FMSC.submodel <- rep(NA, N.blocks)
  K.S.list <- vector(mode = "list", length = N.blocks)
  
  for(i in 1:N.blocks){
    
    #Include only the columns of Z2 that have a corresponding '1' in the ith column of instrument.blocks
    if(N.blocks == 1){
      
      S <- which(instrument.blocks == TRUE)
      
      }else{
        
        S <- which(instrument.blocks[i,] == TRUE)
        
      }#END elseif
    
    Z.S <- cbind(Z1, Z2[,S])
    
    #Set up the matrix Xi.S that extracts elements corresponding to Z.S
    Xi.S <- matrix(0, nrow = (q1 + length(S)), ncol = (q1 + q2))
    Xi.S[1:q1,1:q1] <- diag(q1) #diag defaults to identity matrix
    pi.S <- diag(q2)[S,]
    Xi.S[(q1+1):(q1+length(S)),(q1+1):(q1+q2)] <- pi.S
    
    #Calculate K.S
    C.S <- reg2SLS(X, y, Z.S)$C
    K.S <- -1 * C.S %*% Xi.S
    
    outer.S <- t(nabla.mu.valid) %*% K.S
    FMSC.S <-  outer.S %*% M.inner %*% t(outer.S)
    
    FMSC.submodel[i] <- FMSC.S
    K.S.list[[i]] <- K.S
    
  }#END for
  
  
	#Full Model -- sum of estimated squared Asymptotic Bias and Asymptotic Variance
	AVAR <- N * b.var
	ABIAS.sq <- C %*% Bias.matrix %*% t(C)
	FMSC.full <- t(nabla.mu.valid) %*% (ABIAS.sq + AVAR) %*% nabla.mu.valid
	
	#Valid Model -- by definition this is assumed to have no Asymptotic Bias so the FMSC is simply an estimate of the AVAR of sqrt(N)*(mu.hat - mu.true)
	FMSC.valid <- t(nabla.mu.valid) %*% (N * b1.var) %*% nabla.mu.valid
  
  #Collect and return values of FMSC
  FMSC.submodel <- as.matrix(FMSC.submodel)
  row.names(FMSC.submodel) <- row.names(instrument.blocks)
  FMSC.valid <- as.matrix(FMSC.valid)
  row.names(FMSC.valid) <- c('Valid')
  FMSC.full <- as.matrix(FMSC.full)
  row.names(FMSC.full) <- c('Full')
  FMSC <- rbind(FMSC.valid, FMSC.submodel, FMSC.full)
  colnames(FMSC) <- c('FMSC')

  return(FMSC)
  
}#END FMSC.2SLS
 



library(mvtnorm)

FMSC.2SLS.conf.naive <- function(X, y, Z.valid, Z.additional, instrument.blocks, mu, nabla.mu, Lambda.conf = 0.95, B = 10000){
  
  
    
  Z <- cbind(Z.valid, Z.additional)
  Z1 <- Z.valid
  Z2 <- Z.additional
  
  #Capitalized version of ncol treats an array as a one column matrix
  q1 <- NCOL(Z1)  #Number of valid instruments
	q2 <- NCOL(Z2)	#Number of potentially invalid instruments
	p <- NCOL(X)		#Number of regression parameters
	N <- length(y)	#Sample size
	
	#Fit 2SLS using full instrument set
	full.model <- reg2SLS(X, y, Z)
	b <- full.model$b
	C <- full.model$C
		
	#Estimate covariance matrix S of moment conditions for full model using centered, heteroscedasticity-robust estimator
	u <- as.vector(y - X %*% b)	#Residuals from full model
	D <- diag(u^2) #diag is picky about being passed a vector versus 1-col matrix!
	uu <- outer(u, u)
	Omega <- t(Z) %*% ( (D / N) - (uu / N^2) ) %*% Z
	
	#Estimated Var and SE for b (full model estimate) based on Omega
	b.var <- (1/N) * C %*% Omega %*% t(C)
	
	#Fit 2SLS using valid instrument set
	valid.model <- reg2SLS(X, y, Z1)
	b1 <- valid.model$b
	C1 <- valid.model$C
	
	#Estimate covariance matrix S1 of moment conditions for valid model without centering but allowing for heteroscedasticity
	u1 <- as.vector(y - X %*% b1)
	D1 <- diag(u1^2)
	Omega1 <- t(Z1) %*% ( D1 / N ) %*% Z1
	
	#Estimated Var and SE for b1 (valid model estimate) based on S1
	b1.var <- (1/N) * C1 %*% Omega1 %*% t(C1)
	
	#Construct Asymptotically Unbiased Estimate of outer product of bias vector
	tau <- t(Z2) %*% u1 / sqrt(N) 
	A <- (-1/N) * t(Z2) %*% X %*% C1
	Psi <- cbind(A, diag(nrow = nrow(A))) #Diag defaults to identity matrix
	Sigma <- Psi %*% Omega %*% t(Psi) #Variance matrix of tau.hat
	T.hat <- tau %*% t(tau) - Sigma
	
	#In our calculation of the Asymptotic Bias, T.hat needs to be embedded as the bottom right block in a matrix of the same size as S that has zeros in all other locations. 
	Bias.matrix <- matrix(0, nrow = (q1+q2), ncol = (q1+q2))	
	Bias.matrix[(q1+1):(q1+q2), (q1+1):(q1+q2)] <- T.hat
	
  #Evaluate the derivative of mu at the valid model estimator
  nabla.mu.valid <- nabla.mu(b1)
  #nabla.mu.valid <- c(0,0,1) #Coeff on malfal is the target parameter
  
  #The inner matrix used to calculate the FMSC for all of the potentially mis-specified models
  M.inner <- Bias.matrix + Omega
  
  #Calculate the FMSC and target parameter estimates for the various submodels
  
  #Initialize vector to store the submodel FMSC values and list to store corresponding K.S matrices for confidence interval correction.
  if(is.matrix(instrument.blocks)){N.blocks <- nrow(instrument.blocks)}
  if(is.vector(instrument.blocks)){N.blocks <- 1}
  
  FMSC.submodel <- rep(NA, N.blocks)
  mu.submodel <- rep(NA, N.blocks)
  K.S.list <- vector(mode = "list", length = N.blocks)
  
  for(i in 1:N.blocks){
    
    #Include only the columns of Z2 that have a corresponding '1' in the ith column of instrument.blocks
    if(N.blocks == 1){
      
      S <- which(instrument.blocks == TRUE)
      
      }else{
        
        S <- which(instrument.blocks[i,] == TRUE)
        
      }#END elseif
    
    Z.S <- cbind(Z1, Z2[,S])
    
    #Set up the matrix Xi.S that extracts elements corresponding to Z.S
    Xi.S <- matrix(0, nrow = (q1 + length(S)), ncol = (q1 + q2))
    Xi.S[1:q1,1:q1] <- diag(q1) #diag defaults to identity matrix
    pi.S <- diag(q2)[S,]
    Xi.S[(q1+1):(q1+length(S)),(q1+1):(q1+q2)] <- pi.S
    
    #Carry out 2SLS for the submodel to calculate b.S and K.S
    temp.reg <- reg2SLS(X, y, Z.S)
    b.S <- temp.reg$b
    C.S <- temp.reg$C
    K.S <- -1 * C.S %*% Xi.S
    
    outer.S <- t(nabla.mu.valid) %*% K.S
    FMSC.S <-  outer.S %*% M.inner %*% t(outer.S)
    
    #Estimate of the target parameter
    mu.S <- mu(b.S)
    mu.submodel[i] <- mu.S
    
    FMSC.submodel[i] <- FMSC.S
    K.S.list[[i]] <- K.S
    
  }#END for
  
  
	#Full Model -- sum of estimated squared Asymptotic Bias and Asymptotic Variance
	AVAR <- N * b.var
	ABIAS.sq <- C %*% Bias.matrix %*% t(C)
	FMSC.full <- t(nabla.mu.valid) %*% (ABIAS.sq + AVAR) %*% nabla.mu.valid
	
	#Valid Model -- by definition this is assumed to have no Asymptotic Bias so the FMSC is simply an estimate of the AVAR of sqrt(N)*(mu.hat - mu.true)
	FMSC.valid <- t(nabla.mu.valid) %*% (N * b1.var) %*% nabla.mu.valid
  
  #Collect FMSC values for the various instrument sets
  FMSC.submodel <- as.matrix(FMSC.submodel)
  row.names(FMSC.submodel) <- row.names(instrument.blocks)
  FMSC.valid <- as.matrix(FMSC.valid)
  row.names(FMSC.valid) <- c('Valid')
  FMSC.full <- as.matrix(FMSC.full)
  row.names(FMSC.full) <- c('Full')
  FMSC <- rbind(FMSC.valid, FMSC.submodel, FMSC.full)
  colnames(FMSC) <- c('FMSC')

  #Calculate target paramter estimates for valid and full models
  mu.full <- mu(b)
  mu.valid <- mu(b1)
  mu.submodel <- as.matrix(mu.submodel)
  row.names(mu.submodel) <- row.names(instrument.blocks)
  mu.valid <- as.matrix(mu.valid)
  row.names(mu.valid) <- c('Valid')
  mu.full <- as.matrix(mu.full)
  row.names(mu.full) <- c('Full')
  estimates <- rbind(mu.valid, mu.submodel, mu.full)
  colnames(estimates) <- c('Estimate')
  
  #Create augmented list of K.S matrices including K.valid and K.full
  Xi.full <- diag(q1 + q2)
  K.full <- -1 * C %*% Xi.full
  Xi.valid <- cbind(diag(q1), matrix(0, nrow = q1, ncol = q2))
  K.valid <- -1 * C1 %*% Xi.valid
  K.additional <- list(K.valid, K.full)
  K.list <- c(K.additional, K.S.list)
  
 
  
  #Generate Samples from M evaulated at tau.hat
  M.var <- Omega
  M.mean <- c(rep(0, q1), tau)
  M.sims <- t(rmvnorm(B, M.mean, M.var))
    
  #Initialize vector to store simulated values of Lambda
  Lambda <- rep(NA, B)
    
  #Loop through the simulations from M.sims
  for(j in 1:B){
    
    M.j <- M.sims[,j]
      
    #Lower right block of limit version of bias.matrix
    lower.block <- Psi %*% (M.j %*% t(M.j) - Omega) %*% t(Psi)
      
    #Construct limit version of bias.matrix
    limit.bias.matrix <- matrix(0, nrow = (q1+q2), ncol = (q1+q2))
    limit.bias.matrix[(q1+1):(q1+q2), (q1+1):(q1+q2)] <- lower.block
      
    #Construct limit.inner.matrix with limit.bias.matrix
    limit.inner.matrix <- limit.bias.matrix + Omega
      
    #Now we can get the limit version of the FMSC for each of the submodels by creating a function to apply over K.S.list
    FMSC.limit <- function(K.S){
        
      outer.term <- t(nabla.mu.valid) %*% K.S 
      result <- outer.term %*% limit.inner.matrix %*% t(outer.term)
        
      return(result)
        
    }#END FMSC.limit
      

    FMSC.limit.values <- lapply(K.list, FMSC.limit)
      
    #Choose the model that minimizes the limit FMSC
    K.selected <- K.list[[which.min(FMSC.limit.values)]]

    #Simulation from the limiting distribution of the target parameter, evaulated at the bias value tau
    Lambda[j] <- -1 * nabla.mu.valid %*% K.selected %*% M.j
      
      
    }#end for
    
    #This is a confidence interval for sqrt(n)(mu.hat - mu.true)
    lower.lambda <- quantile(Lambda, (1 - Lambda.conf)/2)
    upper.lambda <- quantile(Lambda, 1 - (1 - Lambda.conf)/2)
    
    #What we want is a confidence interval for mu.hat, so we need to center and re-scale 
    
    #Which estimator minimizes the FMSC? This is the point around which we construct a confidence interval  
    FMSC.results <- data.frame(cbind(FMSC, estimates))
    mu.hat <- FMSC.results$Estimate[which.min(FMSC.results$FMSC)]
  
    #Now we center and scale appropriately to get the desired interval
    lower.mu.hat <- as.numeric(mu.hat - (upper.lambda)/sqrt(N)) 
    upper.mu.hat <- as.numeric(mu.hat - (lower.lambda)/sqrt(N))
    
    post.select <- rbind(lower.mu.hat, mu.hat, upper.mu.hat)
    row.names(post.select) <- c('lower', 'estimate', 'upper')
    
    #Sort FMSC.results by FMSC in ascending order
    FMSC.results <- FMSC.results[order(FMSC, decreasing = TRUE),]
  
    out <- list(FMSC = FMSC.results, CI = post.select)   
  
  
}#END FMSC.2SLS.conf.naive



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


#Now a version of the previous function that produces a valid post-selection confidence interval

#Load packages for optimization subject to nonlinear constraints
library(alabama)
library(DEoptim)
#Load package for simulating from a multivariate normal distribution


FMSC.2SLS.conf <- function(X, y, Z.valid, Z.additional, instrument.blocks, mu, nabla.mu, tau.conf = 0.95, Lambda.conf = 0.95, B = 1000){
  
  Z <- cbind(Z.valid, Z.additional)
  Z1 <- Z.valid
  Z2 <- Z.additional
  
  #Capitalized version of ncol treats an array as a one column matrix
  q1 <- NCOL(Z1)	#Number of valid instruments
	q2 <- NCOL(Z2)	#Number of potentially invalid instruments
	p <- NCOL(X)		#Number of regression parameters
	N <- length(y)	#Sample size
	
	#Fit 2SLS using full instrument set
	full.model <- reg2SLS(X, y, Z)
	b <- full.model$b
	C <- full.model$C
		
	#Estimate covariance matrix S of moment conditions for full model using centered, heteroscedasticity-robust estimator
	u <- as.vector(y - X %*% b)	#Residuals from full model
	D <- diag(u^2) #diag is picky about being passed a vector versus 1-col matrix!
	uu <- outer(u, u)
	Omega <- t(Z) %*% ( (D / N) - (uu / N^2) ) %*% Z
	
	#Estimated Var and SE for b (full model estimate) based on Omega
	b.var <- (1/N) * C %*% Omega %*% t(C)
	
	#Fit 2SLS using valid instrument set
	valid.model <- reg2SLS(X, y, Z1)
	b1 <- valid.model$b
	C1 <- valid.model$C
	
	#Estimate covariance matrix S1 of moment conditions for valid model without centering but allowing for heteroscedasticity
	u1 <- as.vector(y - X %*% b1)
	D1 <- diag(u1^2)
	Omega1 <- t(Z1) %*% ( D1 / N ) %*% Z1
	
	#Estimated Var and SE for b1 (valid model estimate) based on S1
	b1.var <- (1/N) * C1 %*% Omega1 %*% t(C1)
	
	#Construct Asymptotically Unbiased Estimate of outer product of bias vector
	tau <- t(Z2) %*% u1 / sqrt(N) 
	A <- (-1/N) * t(Z2) %*% X %*% C1
	Psi <- cbind(A, diag(nrow = nrow(A))) #Diag defaults to identity matrix
	Sigma <- Psi %*% Omega %*% t(Psi) #Variance matrix of tau.hat
	T.hat <- tau %*% t(tau) - Sigma
	
	#In our calculation of the Asymptotic Bias, T.hat needs to be embedded as the bottom right block in a matrix of the same size as S that has zeros in all other locations. 
	Bias.matrix <- matrix(0, nrow = (q1+q2), ncol = (q1+q2))	
	Bias.matrix[(q1+1):(q1+q2), (q1+1):(q1+q2)] <- T.hat
	
  #Evaluate the derivative of mu at the valid model estimator
  nabla.mu.valid <- nabla.mu(b1)
  #nabla.mu.valid <- c(0,0,1) #Coeff on malfal is the target parameter
  
  #The inner matrix used to calculate the FMSC for all of the potentially mis-specified models
  M.inner <- Bias.matrix + Omega
  
  #Calculate the FMSC and target parameter estimates for the various submodels
  
  #Initialize vector to store the submodel FMSC values and list to store corresponding K.S matrices for confidence interval correction.
  if(is.matrix(instrument.blocks)){N.blocks <- nrow(instrument.blocks)}
  if(is.vector(instrument.blocks)){N.blocks <- 1}
  
  FMSC.submodel <- rep(NA, N.blocks)
  mu.submodel <- rep(NA, N.blocks)
  K.S.list <- vector(mode = "list", length = N.blocks)
  
  for(i in 1:N.blocks){
    
    #Include only the columns of Z2 that have a corresponding '1' in the ith column of instrument.blocks
    if(N.blocks == 1){
      
      S <- which(instrument.blocks == TRUE)
      
      }else{
        
        S <- which(instrument.blocks[i,] == TRUE)
        
      }#END elseif
    
    Z.S <- cbind(Z1, Z2[,S])
    
    #Set up the matrix Xi.S that extracts elements corresponding to Z.S
    Xi.S <- matrix(0, nrow = (q1 + length(S)), ncol = (q1 + q2))
    Xi.S[1:q1,1:q1] <- diag(q1) #diag defaults to identity matrix
    pi.S <- diag(q2)[S,]
    Xi.S[(q1+1):(q1+length(S)),(q1+1):(q1+q2)] <- pi.S
    
    #Carry out 2SLS for the submodel to calculate b.S and K.S
    temp.reg <- reg2SLS(X, y, Z.S)
    b.S <- temp.reg$b
    C.S <- temp.reg$C
    K.S <- -1 * C.S %*% Xi.S
    
    outer.S <- t(nabla.mu.valid) %*% K.S
    FMSC.S <-  outer.S %*% M.inner %*% t(outer.S)
    
    #Estimate of the target parameter
    mu.S <- mu(b.S)
    mu.submodel[i] <- mu.S
    
    FMSC.submodel[i] <- FMSC.S
    K.S.list[[i]] <- K.S
    
  }#END for
  
  
	#Full Model -- sum of estimated squared Asymptotic Bias and Asymptotic Variance
	AVAR <- N * b.var
	ABIAS.sq <- C %*% Bias.matrix %*% t(C)
	FMSC.full <- t(nabla.mu.valid) %*% (ABIAS.sq + AVAR) %*% nabla.mu.valid
	
	#Valid Model -- by definition this is assumed to have no Asymptotic Bias so the FMSC is simply an estimate of the AVAR of sqrt(N)*(mu.hat - mu.true)
	FMSC.valid <- t(nabla.mu.valid) %*% (N * b1.var) %*% nabla.mu.valid
  
  #Collect FMSC values for the various instrument sets
  FMSC.submodel <- as.matrix(FMSC.submodel)
  row.names(FMSC.submodel) <- row.names(instrument.blocks)
  FMSC.valid <- as.matrix(FMSC.valid)
  row.names(FMSC.valid) <- c('Valid')
  FMSC.full <- as.matrix(FMSC.full)
  row.names(FMSC.full) <- c('Full')
  FMSC <- rbind(FMSC.valid, FMSC.submodel, FMSC.full)
  colnames(FMSC) <- c('FMSC')

  #Calculate target paramter estimates for valid and full models
  mu.full <- mu(b)
  mu.valid <- mu(b1)
  mu.submodel <- as.matrix(mu.submodel)
  row.names(mu.submodel) <- row.names(instrument.blocks)
  mu.valid <- as.matrix(mu.valid)
  row.names(mu.valid) <- c('Valid')
  mu.full <- as.matrix(mu.full)
  row.names(mu.full) <- c('Full')
  estimates <- rbind(mu.valid, mu.submodel, mu.full)
  colnames(estimates) <- c('Estimate')

  #Bind together and sort by FMSC?
  
  
  #Constraint Function for optimization problems
  critical <- qchisq(p = tau.conf, df = q2)
  
  #The constraint is the confidence region: 
  #t(tau.test - tau)%*%solve(Sigma)%*%(tau.test - tau) <= critical
  #Which we can rewrite as
  #critical - t(tau.test - tau)%*%solve(Sigma)%*%(tau.test - tau) > 0
  #This puts the constraint in the form we'll need to pass it to the minimization algorithm
  
  
  #To compute this efficiently, we use a Cholesky Decomposition rather than inverting Sigma directly. To be even more efficient, since the constraint function will be called repeatedly, we do the cholesky decomposition OUTSIDE of the function.
  R <- chol(Sigma) #Upper triangular Cholesky factor of Sigma
  
  constraint <- function(tau.test){
    
    tau.diff <- tau.test - tau
    
    #Solve Sigma %*% c = (tau.test - tau) efficiently in two steps
    c <- backsolve(R, forwardsolve(t(R), tau.diff))
    
    #Recall that: c = solve(Sigma) %*% (tau.test - tau)
    #What we want is: t(tau.test - tau) %*% solve(Sigma) %*% (tau.test - tau)
    out <-  critical - t(tau.diff) %*% c
      
    return(out)
    
  }#END constraint
  

  constraint.jacobian <- function(tau.test){
    
  #The derivative of the constraint function is: 
  #-2 * solve(Sigma) %*% (tau.test - tau)
  #But this is simply -2 * c, the intermediate step in the constraint function!
    tau.diff <- tau.test - tau
    c <- backsolve(R, forwardsolve(t(R), tau.diff))
    
    constraint.derivative <- -2 * c
    
    #Since we have a single constraint function, its Jacobian is simply the transpose of its column derivative vector.
    out <- t(constraint.derivative)
    
    return(out)
    
  }#END constraint.jacobian 
    
  #Now we need the objective function. This returns the upper or lower end of a confidence interval for Lambda, evaluated at a given value of tau.test
    
  #Create augmented list of K.S matrices including K.valid and K.full
  Xi.full <- diag(q1 + q2)
  K.full <- -1 * C %*% Xi.full
  Xi.valid <- cbind(diag(q1), matrix(0, nrow = q1, ncol = q2))
  K.valid <- -1 * C1 %*% Xi.valid
  K.additional <- list(K.valid, K.full)
  K.list <- c(K.additional, K.S.list)
  
  
    
  lambda.conf <- function(tau.test, which.tail){

    M.var <- Omega
    M.mean <- c(rep(0, q1), tau.test)
    
    #Generate Samples from M
    M.sims <- t(rmvnorm(B, M.mean, M.var))
    
    #Initialize vector to store simulated values of Lambda
    Lambda <- rep(NA, B)
    
    #Loop through the simulations from M.sims
    for(j in 1:B){
      
      M.j <- M.sims[,j]
      
      #Lower right block of limit version of bias.matrix
      lower.block <- Psi %*% (M.j %*% t(M.j) - S) %*% t(Psi)
      
      #Construct limit version of bias.matrix
      limit.bias.matrix <- matrix(0, nrow = (q1+q2), ncol = (q1+q2))
      limit.bias.matrix[(q1+1):(q1+q2), (q1+1):(q1+q2)] <- lower.block
      
      #Construct limit.inner.matrix with limit.bias.matrix
      limit.inner.matrix <- limit.bias.matrix + Omega
      
      #Now we can get the limit version of the FMSC for each of the submodels by creating a function to apply over K.S.list
      FMSC.limit <- function(K.S){
        
        outer.term <- t(nabla.mu.valid) %*% K.S 
        result <- outer.term %*% limit.inner.matrix %*% t(outer.term)
        
        return(result)
        
      }#END FMSC.limit
      

      FMSC.limit.values <- lapply(K.list, FMSC.limit)
      
      #Choose the model that minimizes the limit FMSC
      K.selected <- K.list[[which.min(FMSC.limit.values)]]

      #Simulation from the limiting distribution of the target parameter, evaulated at the bias value tau
      Lambda[j] <- -1 * nabla.mu.valid %*% K.selected %*% M.j
      
      
    }#end for
    
    #This is a confidence interval for sqrt(n)(mu.hat - mu.true)
    lower.lambda <- quantile(Lambda, (1 - Lambda.conf)/2)
    upper.lambda <- quantile(Lambda, 1 - (1 - Lambda.conf)/2)
    
    if(which.tail == 'lower'){return(lower.lambda)}
    if(which.tail == 'upper'){return(upper.lambda)}
    
  }#END lambda.conf
              

    
  lambda.lower.penalty <- function(tau.test){
    
    tau.diff <- tau.test - tau
    #Solve Sigma %*% c = (tau.test - tau) efficiently in two steps
    c <- backsolve(R, forwardsolve(t(R), tau.diff))
    #Recall that: c = solve(Sigma) %*% (tau.test - tau)
    #What we want is: t(tau.test - tau) %*% solve(Sigma) %*% (tau.test - tau)
    
    penalty <- 0
    
    if((t(tau.diff) %*% c) > critical){penalty <- Inf}
    
    objective <- lambda.conf(tau.test, which.tail = 'lower') + penalty
    
  }#END lambda.lower.penalty
  
  
  #Construct a bounding box for the constraint set by constructing componentwise confidence intervals for tau. Recall that tau.hat has a normal distribution with mean tau and variance matrix Sigma
  z.crit <- qnorm(1 - (1 - tau.conf)/2)
  lower.tau <- tau - z.crit * diag(Sigma)
  upper.tau <- tau + z.crit * diag(Sigma)
  
  DEoptim(fn = lambda.lower.penalty, lower = lower.tau, upper = upper.tau)

    
  
    
  #Now we construct a confidence region for Lambda, which is the limit to which sqrt(N) * (mu.hat - mu.true) converges in distribution. Each value of tau.test gives a different confidence interval. We examine all of the values of tau.test that lie in the confidence region for tau constructed above. 
    
  #To get a conservative interval, we take the min over the lower bounds of all the intervals, and the max of all the upper bounds.
  #starting.value <- rep(0, length(tau)) #auglag doesn't require a feasible starting value
    
    
    #Feasible starting value by construction of confidence region
  
  #To make things easier, construct two new objective functions from lambda.conf: one that gives the lower bound of a confidence interval for Lambda, and one the upper
    #lambda.lower <- function(tau.test){
      
      #return(lambda.conf(tau.test, tail = 'lower'))
      
    #}#END lambda.lower
    
    #lambda.upper <- function(tau.test){
      
      #return(lambda.conf(tau.test, tail = 'upper'))
      
    #}#END lambda.upper
    
    #Carry out the constrained optimization problems
    #lower.lambda <- auglag(par = starting.value, fn = lambda.lower, hin = constraint, hin.jac = constraint.jacobian)$value
    
    #upper.lambda <- auglag(par = starting.value, fn = lambda.upper, hin = constraint, hin.jac = constraint.jacobian)$value
  
  #What we actually want is a confidence interval for mu.hat, the post-selection estimator of the target parameter, not for Lambda. First we need to identify the post-selection estimator, namely the estimator of the target parameter under the instrument set that minimizes FMSC
  FMSC.results <- data.frame(cbind(FMSC, estimates))
  mu.hat <- FMSC.results$Estimate[which.min(FMSC.results$FMSC)]
  
  #Now we center and scale appropriately to get the desired interval
  #lower.mu.hat <- mu.hat - (upper.lambda)/sqrt(N)
  #upper.mu.hat <- mu.hat - (lower.lambda)/sqrt(N)
    
  #Perhaps return results sorted by FMSC? 
  #Should return what the implied minimum coverage is as a function of the user-specified coverages for Lambda and tau
  
}#END FMSC.2SLS.conf

 