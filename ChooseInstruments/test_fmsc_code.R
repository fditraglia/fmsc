setwd("~/fmsc/ChooseInstruments")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("simulation_functions_ChooseInstruments.cpp")

test.dgp <- function(n = 100){
  a0 <- 0.5
  a1 <- 1
  a2 <- 1
  a3 <- 1
  b0 <- 0.75
  b1 <- 1.5 
  z0 <- rep(1, n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  z3 <- rnorm(n)
  x1 <- a0 + a1 * z1 + a2 * z2 + a3 * z3 + rnorm(n)
  x0 <- rep(1, n)
  y <- b0 * x0 + b1 * x1 + rnorm(n)
  out <- list(x0 = x0, x1 = x1, 
              z0 = z0, z1 = z1, z2 = z2, z3 = z3, 
              y = y)
  return(out)
}

set.seed(283)
test.data <- test.dgp()
X <- with(test.data, cbind(x0, x1))
y <- test.data$y
Z1 <- with(test.data, cbind(z0, z1))
Z2 <- with(test.data, cbind(z2, z3))

tsls_est_cpp(X, y, Z1) #Valid model estimates
tsls_est_cpp(X, y, cbind(Z1, Z2)) #Full model estimates
#Additional candidates
tsls_est_cpp(X, y, cbind(Z1, test.data$z2))
tsls_est_cpp(X, y, cbind(Z1, test.data$z3))

#Specify additional candidates using indicators for the columns of Z2 to be used in estimation
add.cand <- cbind(c(1, 0), c(0, 1))
foo <- fmsc_test(X, y, Z1, Z2, as.matrix(0))
bar <- fmsc_test(X, y, Z1, Z2, add.cand)
