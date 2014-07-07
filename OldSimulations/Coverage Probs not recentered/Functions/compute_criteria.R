#September 18th, 2011

#Z1 contains the instruments assumed to be valid whereas Z2 contains the instruments that may be invalid. This code assumes in several places that x is a scalar so that X is really a vector rather an a matrix

#All quantities associated with the valid model, that is the model using only those instruments in Z1, are indicated by a "1" while those for the full model, using all the instruments in Z1 and Z2, are not. For example, b is the 2SLS estimate from the full model, while b1 is the 2SLS estimate from the valid model.

compute.criteria <- function(X, y, Z1, Z2, center.coverage = FALSE){
	
	
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
	full.model <- reg2SLS(X, y, Z, first.stage = TRUE)
	b <- full.model$b
	C <- full.model$C
	X.star <- full.model$X.star #Projection of X into colspace(Z) 
		
	#Estimate covariance matrix S of moment conditions for full model using centered, heteroscedasticity-robust estimator
	u <- as.vector(y - X %*% b)	#Residuals from full model
	D <- diag(u^2) #diag is picky about being passed a vector versus 1-col matrix!
	uu <- outer(u, u)
	S <- t(Z) %*% ( (D / N) - (uu / N^2) ) %*% Z
	
	#J-test Statistic for full model
	Z.u <- t(Z) %*% u
	J <- (1/N) * t(Z.u) %*% solve(S, Z.u)
	
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
	Z1.u <- t(Z1) %*% u1
	J1 <- (1/N) * t(Z1.u) %*% solve(S1, Z1.u)
	
	#Estimated Variance and std dev for b1 (valid model estimate) based on S1
	b1.var <- (1/N) * C1 %*% S1 %*% t(C1)
	b1.SE <- sqrt(b1.var)
	
	
	#-------------------------------------------------------#
	#--ANDREWS (1999) GMM Moment SELECTION CRITERION (MSC)--#
	#-------------------------------------------------------#
	
	#MSC Full Model
	AIC.penalty <- (q - p) * 2 
	BIC.penalty <- (q - p) * log(N)
	HQ.penalty <- (q - p) * 2.01 * log(log(N))
	GMM.AIC <- J - AIC.penalty
	GMM.BIC <- J - BIC.penalty
	GMM.HQ <- J - HQ.penalty
	
	#Valid Model
	AIC.penalty1 <- (q1 - p1) * 2
	BIC.penalty1 <- (q1 - p1) * log(N)
	HQ.penalty1 <- (q1 - p1) * 2.01 * log(log(N))
	GMM.AIC.1 <- J1 - AIC.penalty1
	GMM.BIC.1 <- J1 - BIC.penalty1
	GMM.HQ.1 <- J1 - HQ.penalty1
	
	
	
	#-------------------------------------------------------#
	#--CANONICAL CORRELATIONS INFORMATION CRITERION (CCIC)--#
	#-------------------------------------------------------#
	
	#Note that the CCIC uses the same "penalty terms" as Andrews' GMM-MSC given above but ADDS rather than subtracts these quantities. We choose an instrument by minimizing each criterion. This means that the GMM-MSC "penalties" are really bonus terms since they make the overall criterion smaller. In contrast, these quantities are ADDED to the CCIC, so in this case they are true penalties: they make the criterion larger.
	
	#Calculate first-stage R2 for each model. We use uncentered R2 since the model doesn't contain a constant
	xx <- X %*% X
	Rsq <- (t(X.star) %*% X.star) / xx
	Rsq1 <- (t(X.star1) %*% X.star1) / xx
	
	#CCIC Full Model
	Xi <-  N * log(1 - Rsq) #First term is common to all variants of CCIC 
	CC.AIC <- Xi + AIC.penalty
	CC.BIC <- Xi + BIC.penalty
	CC.HQ <- Xi + HQ.penalty
	
	#CCIC Valid Model
	Xi.1 <- N * log(1- Rsq1) #First term is common to all variants of CCIC 
	CC.AIC.1 <- Xi.1 + AIC.penalty1
	CC.BIC.1 <- Xi.1 + BIC.penalty1
	CC.HQ.1 <- Xi.1 + HQ.penalty1
	
	
	
	#---------------------------------------------#
	#--FOCUSED MOMENT SELECTION CRITERION (FMSC)--#
	#---------------------------------------------#

	#Construct Asymptotically Unbiased Estimate of outer product of bias vector
	tau <- t(Z2) %*% u1 / sqrt(N) 
	A <- (-1/N) * t(Z2) %*% X %*% C1
	M <- cbind(A, diag(nrow = nrow(A)))
	Sigma <- M %*% S %*% t(M)
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
	#--COLLECT AND RETURN CRITERIA AND ESTIMATES--#
	#---------------------------------------------#
	
	out <- data.frame(b, b.SE, J, GMM.AIC, GMM.BIC, GMM.HQ, CC.AIC, CC.BIC, CC.HQ, FMSC, b1, b1.SE, J1, GMM.AIC.1, GMM.BIC.1, GMM.HQ.1, CC.AIC.1, CC.BIC.1, CC.HQ.1, FMSC.1)
	
	return(out)
	
	
}#END compute.criteria