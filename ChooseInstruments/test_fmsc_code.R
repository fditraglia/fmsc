setwd("~/fmsc/ChooseInstruments")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("simulation_functions_ChooseInstruments.cpp")

simple_dgp <- function(){
  n <- 100
  a0 <- 0.5
  a1 <- 1
  a2 <- 1
  b0 <- 0.75
  b1 <- 1.5 
  z0 <- rep(1, n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  x1 <- a0 + a1 * z1 + a2 * z2 + rnorm(n)
  x0 <- rep(1, n)
  y <- b0 * x0 + b1 * x1 + rnorm(n)
  out <- list(x0 = x0, x1 = x1, z0 = z0, z1 = z1, z2 = z2, y = y)
  return(out)
}

test.data <- simple_dgp()
X <- with(test.data, cbind(x0, x1))
y <- test.data$y
Z1 <- with(test.data, cbind(z0, z1))
Z2 <- as.matrix(test.data$z2)
rm(test.data)

tsls_est_cpp(X, y, Z1)
tsls_est_cpp(X, y, cbind(Z1, Z2))

testy <- fmsc_test(X, y, Z1, Z2, as.matrix(0))
