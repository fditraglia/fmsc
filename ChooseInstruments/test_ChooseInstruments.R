setwd("~/fmsc/ChooseInstruments")

library(Rcpp)
library(RcppArmadillo)

sourceCpp("simulation_functions_ChooseInstruments.cpp")

library(sem)


n <- 100
a0 <- 0.5
a1 <- 1
a2 <- 1
b0 <- 0.75
b1 <- 1.5

my_sim <- function(){
  z0 <- rep(1, n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  
  x1 <- a0 + a1 * z1 + a2 * z2 + rnorm(n)
  x0 <- rep(1, n)
  y <- b0 * x0 + b1 * x1 + rnorm(n) 
  
  X <- cbind(x0, x1)
  Z <- cbind(z0, z1, z2)  
  return(tsls_cpp(X, y, Z)) 
}


sem_sim <- function(){  
  z0 <- rep(1, n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  
  x1 <- a0 + a1 * z1 + a2 * z2 + rnorm(n)
  x0 <- rep(1, n)
  y <- b0 * x0 + b1 * x1 + rnorm(n) 
  
  X <- cbind(x0, x1)
  Z <- cbind(z0, z1, z2)  
  return(as.matrix(tsls(y ~ x1, ~ z1 + z2)$coef)) 
}


set.seed(7436)
system.time(foo <- replicate(1000, sem_sim()))

set.seed(7436)
system.time(bar <- replicate(1000, sem_sim()))

all.equal(foo, bar)
