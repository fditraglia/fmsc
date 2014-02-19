setwd("~/fmsc/ChooseInstruments")

library(Rcpp)
library(RcppArmadillo)

sourceCpp("simulation_functions_ChooseInstruments.cpp")

library(sem)

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

cpp_est_sim <- function(){
  sim_data <- simple_dgp()
  X <- with(sim_data, cbind(x0, x1))
  Z <- with(sim_data, cbind(z0, z1, z2))
  y <- sim_data$y
  out <- as.vector(tsls_est_cpp(X, y, Z))
  return(out)
}

r_est_sim <- function(){  
  sim_data <- simple_dgp()
  out <- tsls(y ~ x1, ~ z1 + z2, data = sim_data)$coef
  names(out) <- NULL
  return(out)
}


set.seed(7436)
system.time(foo <- replicate(1000, r_est_sim()))

set.seed(7436)
system.time(bar <- replicate(1000, cpp_est_sim()))

all.equal(foo, bar)




cpp_se_sim <- function(){
  sim_data <- simple_dgp()
  X <- with(sim_data, cbind(x0, x1))
  Z <- with(sim_data, cbind(z0, z1, z2))
  y <- sim_data$y
  out <- as.vector(tsls_SE_cpp(X, y, Z))
  return(out)  
}


r_se_sim <- function(){
  sim_data <- simple_dgp()
  out <- tsls(y ~ x1, ~ z1 + z2, data = sim_data)$V
  out <- sqrt(diag(out))
  names(out) <- NULL
  return(out)
}


set.seed(7436)
system.time(foo <- replicate(1000, r_se_sim()))

set.seed(7436)
system.time(bar <- replicate(1000, cpp_se_sim()))

all.equal(foo, bar)

