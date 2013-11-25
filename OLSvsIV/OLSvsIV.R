#Filename:        OLSvsIV.R
#Author:          Frank DiTraglia
#First Version:   2013-21-11
#This Version:    2013-22-11

#This script is a preliminary implementation of the OLS versus IV example that I am adding to the fmsc paper.


dgp <- function(b, PI, V.e, V.z, n){
  
  #------------------------------------------------------------------
  #Arguments:
  # b       coefficient on x in the second stage
  # PI      vector of coefficients on z in the second stage
  # V.e     variance-covariance matrix of the errors epsilon and v
  # V.z     variance-covariance matrix of the instruments z
  # n       sample size
  #
  #Returns:
  # data    list of matrices x, y and z containing simulated dataset
  #------------------------------------------------------------------ 
  
  n.z <- length(PI)

  errors <- chol(V.e) %*% matrix(rnorm(2 * n), 2, n)  
  z <- chol(V.z) %*% matrix(rnorm(n.z * n), n.z, n)
  
  x <- t(z) %*% PI + t(errors)[,2]
  y <- b * x + t(errors)[,1]
  
  data <- list(x = x, y = y, z = t(z))
  return(data)
  
}




#Faster C++ version of dgp: requires RcppArmadillo and setwd("~/fmsc/OLSvsIV")
sourceCpp("dgp.cpp") 




fmsc.ols.iv <- function(x, y, z, DHW.levels = NULL){
  
  #------------------------------------------------------------------
  #NOTE: This function assumes that x is a column of observations for 
  #      a single endogenous regressor. In other words, it assumes 
  #      that all exogenous regressors, including the constant, have 
  #      been "projected out" of the system.
  #------------------------------------------------------------------
  #Arguments:
  # x               matrix of observations of single endog. regressor
  # y               matrix of observations for outcome variable
  # z               matrix of observations for the instruments
  # DHW.levels      optional vector of significance levels for a 
  #                   Durbin-Hausman-Wu tests. If specified, the 
  #                   function returns, in addition to the OLS, IV 
  #                   and FMSC-selected estimators, the corresponding 
  #                   DHW pretest estimators.
  #------------------------------------------------------------------
  
  n <- length(y)
  
  xx <- crossprod(x)
  b.ols <- crossprod(x,y) / xx
  #lm(y ~ x - 1)
  
  zx <- crossprod(z, x)
  zz.inv <- chol2inv(qr.R(qr(z))) 
  x.hat <- z %*% (zz.inv %*% zx)
  b.tsls <- crossprod(x.hat, y) / crossprod(x.hat)
  #tsls.fit <- tsls(y ~ x - 1, ~ z - 1)
  
  tsls.residuals <- y - x %*% b.tsls
  s.e.squared <- crossprod(tsls.residuals) / n
  #crossprod(tsls.fit$residuals) / n
  
  tau <- crossprod(x, tsls.residuals) / sqrt(n)
  #tau <- crossprod(x, y - x %*% b.tsls) / sqrt(n)
  
  s.x.squared <- xx / n
  g.squared <- t(zx) %*% zz.inv %*% zx / n
  s.v.squared <- s.x.squared - g.squared
  
  Tfmsc <- tau^2 * g.squared / (s.v.squared * s.e.squared * s.x.squared)
  #This is numerically equivalent to a Durbin-Hausman-Wu Test Statistic
  #n * (b.ols - b.tsls)^2 / (s.e.squared * (1 / g.squared - 1 / s.x.squared))
  
  #Estimator chosen by FMSC
  b.fmsc <- ifelse(Tfmsc < 2, b.ols, b.tsls)
  
  #If DHW-test levels are specified, calculate the pre-test estimators
  if(!is.null(DHW.levels)){
    DHW.crit <- qchisq(1 - DHW.levels, 1)
    DHW <- (as.vector(Tfmsc) <= DHW.crit) * b.ols + (as.vector(Tfmsc) > DHW.crit) * b.tsls
    names(DHW) <- paste('b.DHW.', 100 *DHW.levels, sep ="")
  }else{
    DHW <- NULL
  }

  
  #Estimated Optimal Weight for OLS - plug in asymptotically unbiased estimator of tau-squared
  tau.squared.est <- tau^2 - s.e.squared * s.x.squared * s.v.squared / g.squared
  bias.est <- max(0, tau.squared.est / s.x.squared^2)
  var.diff <- s.e.squared * (1/g.squared - 1/s.x.squared)
  omega.star <- 1 / (1 + (bias.est / var.diff))
  if(omega.star > 1){omega.star <- 1}
  if(omega.star < 0){omega.start <- 0}
  
  b.star <- omega.star * b.ols + (1 - omega.star) * b.tsls
  
  out <- c(b.ols, b.tsls, b.fmsc, b.star) 
  names(out) <- c('b.ols', 'b.tsls', 'b.fmsc', 'b.star')
  
  #Append DHW pretest estimators (NULL by default) 
  out <- c(out, DHW)
 
  return(out)
  
}




#Potentitally faster C++ version
sourceCpp("fmsc_ols_iv.cpp")





simple.sim <- function(p, r, n){
  
  sim.data <- dgp_cpp(1, p * rep(1, 3), matrix(c(1, r, r, 1), 2, 2), diag(rep(1, 3)), n)
  b <- fmsc.ols.iv(sim.data$x, sim.data$y, sim.data$z)
  return(b)
}


mse <- function(x, truth){mean((x - truth)^2)}



mse.compare <- function(p, r, n, n.reps = 10000){
  
  sim.results <- replicate(n.reps, simple.sim(p, r, n))
  out <- t(apply(sim.results, 1, mse, truth = 1))
  return(out)
}



r.seq <- seq(0, 0.2, 0.01)

set.seed(4938)
fooCpp <- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse.compare(p = 0.3, r, n = 250)}, mc.set.seed = FALSE))#mclapply outputs a list of vectors. Combine them into a matrix.

matplot(r.seq, apply(fooCpp, 2, sqrt))
        #, col = c('black', 'red', 'blue'), xlab = 'Cor(e,v)', ylab = 'RMSE', type =  'l', lty = 1)
# legend("topleft", c("FMSC", "OLS", "IV"), fill = c("black", "red", "blue"))




