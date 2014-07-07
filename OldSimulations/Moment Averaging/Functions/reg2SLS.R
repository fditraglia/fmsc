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

