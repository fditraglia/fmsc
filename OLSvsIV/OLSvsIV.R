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



fmsc.ols.iv <- function(x, y, z){
  
  #------------------------------------------------------------------
  #NOTE: This function assumes that x is a column of observations for 
  #      a single endogenous regressor. In other words, it assumes 
  #      that all exogenous regressors, including the constant, have 
  #      been "projected out" of the system.
  #------------------------------------------------------------------
  #Arguments:
  # x       matrix of observations of single endogenous regressor
  # y       matrix of observations for outcome variable
  # z       matrix of observations for the instruments
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
  
  #Estimated chosen by FMSC
  b.fmsc <- ifelse(Tfmsc < 2, b.ols, b.tsls)
  
  #Estimated Optimal Weight for OLS
  #tau.squared.est <- tau^2 - s.e.squared * s.x.squared * s.v.squared / g.squared
  #abias.squared.est.ols <- tau.squared.est / s.x.squared^2
  #numerator <- max(0, abias.squared.est.ols)
  #avar.est.tsls <- s.e.squared / g.squared
  #avar.est.ols <- s.e.squared / s.x.squared
  #denominator <- avar.est.ols - avar.est.tsls
  #omega.star <- 1 / (1 - (numerator / denominator))
  #b.star <- omega.star * b.ols + (1 - omega.star) * b.tsls
  
  out <- c(b.fmsc, b.ols, b.tsls) 
  names(out) <- c('b.fmsc', 'b.ols', 'b.tsls')
  return(out)
  
}


simple.sim <- function(p, r, n){

  PI <- p * rep(1, 3)
  V.z <- diag(rep(1, 3))
  V.e <- diag(rep(1, 2)) + matrix(c(0, r, r, 0), 2, 2)
  
  sim.data <- dgp(1, PI, V.e, V.z, n)
  x <- sim.data$x
  y <- sim.data$y
  z <- sim.data$z
  b <- fmsc.ols.iv(x, y, z)
  return(b)
  
}


mse <- function(x, truth){mean((x - truth)^2)}


#This function runs in parallel on Linux/Mac machines with more than one core using mclapply. It assumes that the library multicore has been loaded.
mse.compare <- function(p, r, n, n.reps = 10000){
  
  sim.results <- lapply(X = 1:n.reps, FUN = function(.){simple.sim(p, r, n)}) 
  sim.results <- do.call(rbind, sim.results) #mclapply outputs a list of vectors. Combine them into a matrix.
  out <- apply(sim.results, 2, mse, truth = 1)  
  return(out)
}


#(Potentially) Faster C++ version of dgp. Load Rcpp and RcppArmadillo first. Also set the appropriate directory: "~/fmsc/OLSvsIV"
sourceCpp("dgp.cpp")

#Corresponding version of simple.sim
simple.sim.cpp <- function(p, r, n){
  
  PI <- p * rep(1, 3)
  V.z <- diag(rep(1, 3))
  V.e <- diag(rep(1, 2)) + matrix(c(0, r, r, 0), 2, 2)
  
  sim.data <- dgp_cpp(1, PI, V.e, V.z, n)
  x <- sim.data$x
  y <- sim.data$y
  z <- sim.data$z
  b <- fmsc.ols.iv(x, y, z)
  return(b)
  
}

#Corresponding version of mse.compare
mse.compare.cpp <- function(p, r, n, n.reps = 10000){
  
  sim.results <- lapply(X = 1:n.reps, FUN = function(.){simple.sim.cpp(p, r, n)}) 
  sim.results <- do.call(rbind, sim.results) #mclapply outputs a list of vectors. Combine them into a matrix.
  out <- apply(sim.results, 2, mse, truth = 1)  
  return(out)
}


set.seed(3029)
mse.compare.cpp(0.1, 0.2, 100, 10000)
set.seed(3029)
mse.compare.cpp(0.1, 0.2, 1000, 10000)

#Example of the kind of plot I'll use in the paper
#r.seq <- seq(0, 0.2, 0.01)
#set.seed(3728)
#mse.values <- t(mapply(mse.compare, p = 0.1, r = r.seq, n = 250))
#set.seed(3728)
#mse.values.cpp <- t(mapply(mse.compare.cpp, p = 0.1, r = r.seq, n = 250))
#set.seed(3728)
#mse.values.cpp.alt <- t(mapply(mse.compare.cpp.alt, p = 0.1, r = r.seq, n = 250))


# 
# matplot(r.seq, apply(mse.values, 2, sqrt), col = c('black', 'red', 'blue'), xlab = 'Cor(e,v)', ylab = 'RMSE', type =  'l', lty = 1)
# legend("topleft", c("FMSC", "OLS", "IV"), fill = c("black", "red", "blue"))
# 
# matplot(r.seq, apply(mse.values.cpp, 2, sqrt), col = c('black', 'red', 'blue'), xlab = 'Cor(e,v)', ylab = 'RMSE', type =  'l', lty = 1)
# legend("topleft", c("FMSC", "OLS", "IV"), fill = c("black", "red", "blue"))
# 
# matplot(r.seq, apply(mse.values.cpp.alt, 2, sqrt), col = c('black', 'red', 'blue'), xlab = 'Cor(e,v)', ylab = 'RMSE', type =  'l', lty = 1)
# legend("topleft", c("FMSC", "OLS", "IV"), fill = c("black", "red", "blue"))




