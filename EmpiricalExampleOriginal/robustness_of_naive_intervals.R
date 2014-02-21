#------------------------------------------------------#
#-------------------FUNCTION: reg2SLS------------------#
#------------------------------------------------------#
#This function returns the 2SLS regression coefficient b for a regression of y on X with instruments Z, as well as the product needed for estimating the variance matrix. If first.stage is TRUE, it also returns the ``transformed'' X variables that result from the first-stage regression. I have updated this to be even more efficient by using the functions crossprod and tcrossprod. The old version is given in the comments.
reg2SLS <- function(X, y, Z, first.stage = FALSE){
  
  #Calculate QR Decomposition of Z
	QR.z <- qr(Z)
	Qz <- qr.Q(QR.z)
	Rz <- qr.R(QR.z)
	
	#Calculate (Z'Z)^{-1}
	Rz.inv <- backsolve(Rz, diag(1, nrow = nrow(Rz), ncol = ncol(Rz)))
	ZZ.inv <- tcrossprod(Rz.inv)
  #ZZ.inv <- Rz.inv %*% t(Rz.inv)
	
	#Transform X and y
	X.star <- crossprod(Qz, X)
  y.star <- crossprod(Qz, y)
  #X.star <- t(Qz) %*% X
	#y.star <- t(Qz) %*% y
	
	#Calculate QR decomposition of transformed X
	QR<- qr(X.star)
	Q <- qr.Q(QR)
	R <- qr.R(QR)
	
	#Calculate (X'X)^{-1} for the transformed X
	R.inv <- backsolve(R, diag(1, nrow = nrow(R), ncol = ncol(R)))
	XX.star.inv <- tcrossprod(R.inv)
  #XX.star.inv <- R.inv %*% t(R.inv)
	
	#Calculate C, the matrix used to estimate the variance via CSC'
	N <- length(y)
	C <- N * XX.star.inv %*% crossprod(X, Z) %*% ZZ.inv
  #C <- N * XX.star.inv %*% t(X) %*% Z %*% ZZ.inv
	
	#Calculate 2SLS coefficient
	c <- crossprod(Q, y.star)
  #c <- t(Q) %*% y.star
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




#------------------------------------------------------#
#-----------FUNCTION: FMSC.2SLS.return.Ks--------------#
#------------------------------------------------------#

FMSC.2SLS.return.Ks <- function(X, y, Z.valid, Z.additional, instrument.blocks, mu, nabla.mu){
    
  Z <- cbind(Z.valid, Z.additional)
  Z1 <- Z.valid
  Z2 <- Z.additional
  
  #Capitalized version of ncol treats an array as a one column matrix
  q1 <- NCOL(Z1)  #Number of valid instruments
  q2 <- NCOL(Z2)  #Number of potentially invalid instruments
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
    
  #Which estimator minimizes the FMSC? This is the point around which we construct a confidence interval  
  FMSC.results <- data.frame(cbind(FMSC, estimates))
  mu.hat <- FMSC.results$Estimate[which.min(FMSC.results$FMSC)]
  
  out <- list(FMSC.results, mu.hat, K.list, Omega, Psi, tau, q1, q2, nabla.mu.valid, N)
  return(out)
  
}#END FMSC.2SLS.return.Ks


#------------------------------------------------------#
#---------END FUNCTION: FMSC.2SLS.return.Ks------------#
#------------------------------------------------------#







conf.robustness <- function(X, y, Z.valid, Z.additional, instrument.blocks, mu, nabla.mu, B = 10000){
  
  first.step <- FMSC.2SLS.return.Ks(X, y, Z.valid, Z.additional, instrument.blocks, mu, nabla.mu)
  
  FMSC.results <- first.step[[1]]
  mu.hat <- first.step[[2]]
  K.list <- first.step[[3]]
  Omega <- first.step[[4]]
  Psi <- first.step[[5]]
  tau.hat <- first.step[[6]]
  q1 <- first.step[[7]]
  q2 <- first.step[[8]]
  nabla.mu.valid <- first.step[[9]]
  N <- first.step[[10]]
  
  Sigma <- Psi %*% Omega %*% t(Psi)
  
  #Confidence level for Lambda confidence interval
  Lambda.conf <- 0.95

  mu.conf <- function(tau.test){

    #Generate Samples from M evaulated at tau.test
    M.var <- Omega
    M.mean <- c(rep(0, q1), tau.test)
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
  
      #convert to a confidence interval for mu
      lower.mu <- mu.hat - (upper.lambda)/sqrt(N)
      upper.mu <- mu.hat - (lower.lambda)/sqrt(N)
  
      return(data.frame(lower = lower.mu, upper = upper.mu))
  
  }#END mu.conf
  
  #Endpoints of element-wise 95% confidence intervals for tau
  z.crit <- qnorm(1 - 0.5/2)
  lower.tau <- tau.hat - z.crit * diag(Sigma)
  upper.tau <- tau.hat + z.crit * diag(Sigma)

  
  
  tau.zero <- mu.conf(rep(0, length(tau.hat)))
  tau.est <- mu.conf(tau.hat)
  tau.max <- mu.conf(upper.tau)
  tau.min <-  mu.conf(lower.tau)
  
  out <- rbind(tau.zero, tau.est, tau.max, tau.min)
  row.names(out) <- c('tau.0', 'tau.hat', 'tau.upper', 'tau.lower')
  return(out)
  
}#END conf.robustness



#-------------------------------------------------------#
#----------------------PRELIMINARIES--------------------#
#-------------------------------------------------------#

#Set Working Directory
#office computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Empirical Examples/Carstensen and Gundlach (Development)"

#Home computer
#root.dir <- "/Users/Frank/Dropbox/PhD Dissertation/Job Market Paper/Empirical Examples/Carstensen and Gundlach (Development)"

#Load package for two-stage least squares
library(sem)

#Load my functions for FMSC moment selection and confidence interval correction
#setwd(root.dir)
#source("FMSC_linear_2SLS.R")


#-------------------------------------------------------#
#-----------------LOAD AND CLEAN DATASET----------------#
#-------------------------------------------------------#

#Load data: NAs are coded as -999.99
input.data <- read.csv("Carstensen_Gundlach.csv", header = TRUE, na.strings = "-999.999")

#Change country from factor data to character data
#is.character(data$country) #FALSE
input.data$country <- as.character(input.data$country)

#Some of the variable names in the dataset don't match those in the paper. Change them so they do.
attach(input.data)
rule <- kaufman
malfal <- mfalrisk
exprop <- exprop2
lngdpc <- lngdpc95
trade <- frarom
latitude <- lat
coast <- landsea
clean.data <- data.frame(country, lngdpc, rule, malfal, malrisk, gadp, exprop, lnmort, maleco, frost, humid, latitude, eurfrac, engfrac, coast, trade)
detach(input.data)
clean.data$country <- as.character(clean.data$country)




#-------------------------------------------------------#
#-----------------FMSC MOMENT SELECTION-----------------#
#-------------------------------------------------------#

#Instrument blocks considered in Table 2 of the paper
attach(clean.data)
baseline <- cbind(lnmort, maleco)
openness <- cbind(trade, coast)
climate <- cbind(frost, humid, latitude)
europe <- cbind(eurfrac, engfrac)
detach(clean.data)


#My functions use matrix notation rather than R's formula objects to specify the 2SLS model, so I need to create a matrix of all complete cases of the instruments and regressors I'll use.
Full.data <- na.omit(data.frame(lngdpc = clean.data$lngdpc, rule = clean.data$rule, malfal = clean.data$malfal, baseline, climate, europe, openness))

#The outcome variable: log real gdp per capita at PPP in 1995 in international dollars
y <- Full.data[,1]

#Include a constant term by appending a column of ones
constant <- rep(1, NROW(Full.data))

#The regressors
X <- data.frame(constant, Full.data[,2:3])
#The columns of X are: constant, rule, malfal
X <- as.matrix(X) #My functions require matrix input

#All exogenous regressors go in the instrument set: here only the constant is exogenous (we instrument both rule and malfal)

#The baseline instruments are assumed to be valid
Z.valid <- data.frame(constant, Full.data[,4:5]) 
#The columns of Z.valid are: constant, lnmort, maleco
Z.valid <- as.matrix(Z.valid) #My functions require matrix input

#The additional instruments may be invalid
Z.additional <- data.frame(Full.data[,6:12])
#The columns of Z.additional are: 
#frost, humid, latitude, eurfrac, engfrac, trade, coast
Z.additional <- as.matrix(Z.additional) #My functions require matrix input

#Encode the instrument blocks used in Table 2: each row of this matrix corresponds to the columns of Z.additional (the potentially invalid instruments that we're considering including). A zero in a given column means that the instrument in that column of Z.additional is not included; a one means that it is included.

#Recall that the columns of Z.additional are:
#frost, humid, latitude, eurfrac, engfrac, trade, coast
#Thus we have:
#climate block: columns 1-3 of Z.additional
#Europe block: columns 4-5 of Z.additional
#openness block: columns 6-7 of Z.additional

climate.instruments <- c(1, 1, 1, 0, 0, 0, 0) #include climate block
europe.instruments <- c(0, 0, 0, 1, 1, 0, 0) #include Europe block
openness.instruments <- c(0, 0, 0, 0, 0, 1, 1) #include openness block

additional.instruments <- as.matrix(rbind(climate.instruments,
                                                europe.instruments,
                                                openness.instruments)) 
#We don't need to specify the baseline model or the model including all instruments, as these are included automatically by my function.

#In the paper, the stated target is the disease-income link, that is the coefficient on malfal. Since the columns of X are constant, rule, malfal this means we're interested in the third coefficient
mu.malfal <- function(x){return(x[3])}
nabla.malfal <- function(x){return(c(0, 0, 1))}


#What if the target were rule?
mu.rule <- function(x){return(x[2])}
nabla.rule <- function(x){return(c(0, 1, 0))}


#In each of these cases, the FMSC tells us to use the full model.



#Variation #1: what if we consider including some of the instrument blocks together?
climate.only <- c(1, 1, 1, 0, 0, 0, 0) 
europe.only <- c(0, 0, 0, 1, 1, 0, 0) 
openness.only <- c(0, 0, 0, 0, 0, 1, 1)
climate.europe <- c(1, 1, 1, 1, 1, 0, 0)
climate.openness <- c(1, 1, 1, 0, 0, 1, 1)
europe.openness <- c(0, 0, 0, 1, 1, 1, 1)

additional.instruments.var1 <- as.matrix(rbind(climate.only,
                                               europe.only,
                                               openness.only,
                                               climate.europe,
                                               climate.openness,
                                               europe.openness))




#Variation #2: What about at about including squares of the endogenous regressors? 
#Recall that the columns of X are: constant, rule, malfal
rule.sq <- X[,2]^2
malfal.sq <- X[,3]^2

Z.additional.var2 <- cbind(rule.sq, malfal.sq)
#The columns of Z.additional.var2 are: rule.sq, malfal.sq
Z.additional.var2 <- as.matrix(Z.additional.var2)

rule.squared <- c(1, 0)
malfal.squared <- c(0, 1)

additional.instruments.var2 <- as.matrix(rbind(rule.squared,
                                               malfal.squared))








#Variation #3: What if we also consider the instrument blocks from Table 2?

Z.additional.var3 <- cbind(Z.additional, Z.additional.var2)
#The columns of Z.additional are: 
#frost, humid, latitude, eurfrac, engfrac, trade, coast
#The columns of Z.additional.var2 are: rule.sq, malfal.sq
#Hence, the columns of Z.additional.var3 are:
#frost, humid, latitude, eurfrac, engfrac, trade, coast, rule.sq, malfal.sq

climate.block <-    c(1, 1, 1, 0, 0, 0, 0, 0, 0) 
europe.block <-     c(0, 0, 0, 1, 1, 0, 0, 0, 0) 
openness.block <-   c(0, 0, 0, 0, 0, 1, 1, 0, 0)
rule.sq <-          c(0, 0, 0, 0, 0, 0, 0, 1, 0)
malfal.sq <-        c(0, 0, 0, 0, 0, 0, 0, 0, 1)
endog.sq <-         c(0, 0, 0, 0, 0, 0, 0, 1, 1)

additional.instruments.var3 <- as.matrix(rbind(climate.block,
                                               europe.block,
                                               openness.block,
                                               rule.sq,
                                               malfal.sq,
                                               endog.sq))


                                               
                                               
                                               
       
#FMSC results and naive confidence intervals
FMSC.malfal <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional, additional.instruments, mu.malfal, nabla.malfal)

FMSC.rule <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional, additional.instruments, mu.rule, nabla.rule)                                            
                                               
FMSC.malfal.var1 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional, additional.instruments.var1, mu.malfal, nabla.malfal)

FMSC.rule.var1 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional, additional.instruments.var1, mu.rule, nabla.rule)                                           
                                               
                                               
FMSC.malfal.var2 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional.var2, additional.instruments.var2, mu.malfal, nabla.malfal)

FMSC.rule.var2 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional.var2, additional.instruments.var2, mu.rule, nabla.rule) 
                                               

FMSC.malfal.var3 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional.var3, additional.instruments.var3, mu.malfal, nabla.malfal)

FMSC.rule.var3 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional.var3, additional.instruments.var3, mu.rule, nabla.rule)

    
#Now, how robust are these confidence intervals to differences in tau?
intervals.malfal <- conf.robustness(X, y, Z.valid, Z.additional, additional.instruments, mu.malfal, nabla.malfal)                                   

intervals.rule <- conf.robustness(X, y, Z.valid, Z.additional, additional.instruments, mu.rule, nabla.rule)
                                   
intervals.malfal.var1 <- conf.robustness(X, y, Z.valid, Z.additional, additional.instruments.var1, mu.malfal, nabla.malfal)

intervals.rule.var1 <- conf.robustness(X, y, Z.valid, Z.additional, additional.instruments.var1, mu.rule, nabla.rule)                                           
                                                                                      
intervals.malfal.var2 <- conf.robustness(X, y, Z.valid, Z.additional.var2, additional.instruments.var2, mu.malfal, nabla.malfal)

intervals.rule.var2 <- conf.robustness(X, y, Z.valid, Z.additional.var2, additional.instruments.var2, mu.rule, nabla.rule) 
                                               
intervals.malfal.var3 <- conf.robustness(X, y, Z.valid, Z.additional.var3, additional.instruments.var3, mu.malfal, nabla.malfal)

intervals.rule.var3 <- conf.robustness(X, y, Z.valid, Z.additional.var3, additional.instruments.var3, mu.rule, nabla.rule)
                                               
                                              
                                            
