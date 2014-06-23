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

#Specify additional candidates using indicators for the columns of Z2 to be used in estimation
add.cand <- cbind(c(1, 0), c(0, 1))
foo <- fmsc_test(X, y, Z1, Z2, as.matrix(0))
bar <- fmsc_test(X, y, Z1, Z2, add.cand)

#Check that we're getting the right estimates
#stored in the correct order
all.equal(foo$Estimates, bar$Estimates[,c(1,4)])

est.valid <- tsls_est_cpp(X, y, Z1) #Valid model estimates
est1 <- tsls_est_cpp(X, y, cbind(Z1, test.data$z2)) #Add. candidate #1
est2 <- tsls_est_cpp(X, y, cbind(Z1, test.data$z3)) #Add. candidate #2
est.full <- tsls_est_cpp(X, y, cbind(Z1, Z2)) #Full model estimates
est <- cbind(est.valid, est1, est2, est.full)

all.equal(est, bar$Estimates)

#tau, Psi, tau.outer and Bias.mat
#should be the same regardless of
#whether we include add.cand
all.equal(foo$tau, bar$tau)
all.equal(foo$Psi, bar$Psi)
all.equal(foo$tau.outer, bar$tau.outer)
all.equal(foo$Bias.mat, bar$Bias.mat)

#Valid and Full model avar.inner
#and sq.bias.inner should be the same
#regardless of whether we include
#add.cand as should K and Omega
all.equal(foo$avar.inner, bar$avar.inner[,,c(1,4)])
all.equal(foo$sq.bias.inner, bar$sq.bias.inner[,,c(1,4)])
all.equal(foo$K[[1]], bar$K[[1]])
all.equal(foo$K[[2]], bar$K[[4]])


#Now that I knows the results are consistent between foo and bar, I need to check that the calculations are in fact correct. I'll do this by directly implementing the FMSC calculations inline

#Calculation of tau
n <- length(y)
tau <- t(Z2) %*% (y - X %*% est.valid) / sqrt(n)
all.equal(tau, bar$tau, check.attributes = FALSE)

#Function to calculate K - Neither efficient nor particularly stable but fine for testing this simple example and identical to the formula from the paper
K.calc <- function(Z){
  ZZ.inv <- solve(t(Z) %*% Z)
  P.Z <- Z %*% ZZ.inv %*% t(Z)
  out <- n * solve(t(X) %*% P.Z %*% X) %*% t(X) %*% Z %*% ZZ.inv
  colnames(out) <- NULL
  rownames(out) <- NULL
  return(out)
}

#Test K matrices
K.valid <- K.calc(Z1)
K.1 <- K.calc(cbind(Z1, test.data$z2))
K.2 <- K.calc(cbind(Z1, test.data$z3))
K.full <- K.calc(cbind(Z1, Z2))

all.equal(K.valid, bar)

#Calculation of Psi
Psi <- cbind(-1 * t(Z2) %*% BLAH , diag(1, ncol(Z2)))

