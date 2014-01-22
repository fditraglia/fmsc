#--------------------------------------------------------------
# Filename:        simulation_functions_OLSvsIV.R
# Author:          Frank DiTraglia
# First Version:   2014-01-22
# This Version:    2013-11-22
#--------------------------------------------------------------
# R code to carry out the OLS vs. IV simulation experiments.
# All functions here correspond to C++ versions contained in 
# the file simulation_functions_OLSvsIV.cpp.
#--------------------------------------------------------------


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
  #tsls.fit <- tsls(y ~ x - 1, ~ z - 1) #Debugging code
  
  tsls.residuals <- y - x %*% b.tsls
  s.e.squared <- crossprod(tsls.residuals) / n
  #crossprod(tsls.fit$residuals) / n
  
  tau <- crossprod(x, tsls.residuals) / sqrt(n)
  #tau <- crossprod(x, y - x %*% b.tsls) / sqrt(n) #Debugging code
  
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





simple.sim <- function(p, r, n){
  
  sim.data <- dgp(1, p * rep(1, 3), matrix(c(1, r, r, 1), 2, 2), diag(rep(1, 3)), n)
  b <- fmsc.ols.iv(sim.data$x, sim.data$y, sim.data$z)
  return(b)
}



mse <- function(x, truth){mean((x - truth)^2)}



mse.compare <- function(p, r, n, n.reps = 10000){
  
  sim.results <- replicate(n.reps, simple.sim(p, r, n))
  out <- t(apply(sim.results, 1, mse, truth = 1))
  return(out)
}



#Faster version that calls a C++ version
mse.compare.cpp <- function(p, r, n, n.reps = 10000){
  
  sim.results <- replicate(n.reps, simple_sim_cpp(p, r, n))
  out <- t(apply(sim.results, 1, mse, truth = 1))
  return(out)
}
