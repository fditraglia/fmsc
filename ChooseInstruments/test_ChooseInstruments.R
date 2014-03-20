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


#Test C++ estimates against tsls
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



#Test C++ textbook standard errors against tsls
cpp_se_sim <- function(){
  sim_data <- simple_dgp()
  X <- with(sim_data, cbind(x0, x1))
  Z <- with(sim_data, cbind(z0, z1, z2))
  y <- sim_data$y
  out <- as.vector(tsls_SE_textbook_cpp(X, y, Z))
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

#tsls doesn't seem to have an option for robust or centered standard errors so we can't test the C++ against it. However, for this DGP the robust and centered standard errors should be very close to the textbook ones.
set.seed(821)
sim_data <- simple_dgp()
X <- with(sim_data, cbind(x0, x1))
Z <- with(sim_data, cbind(z0, z1, z2))
y <- sim_data$y

tsls_SE_textbook_cpp(X, y, Z)
tsls_SE_robust_cpp(X, y, Z)
tsls_SE_center_cpp(X, y, Z)


#Test the dgp function
set.seed(352)
test_dgp(0.2, 0.1, 100)


#Test CCIC class
set.seed(389)
baz <- CCIC_test(1,0.1)

set.seed(389)
testy <- dgp_cpp(1,0.1)
cc <- cancor(testy$x, cbind(testy$z1, testy$z2))
n <- length(testy$x)
r <- cc$cor
bar <- lm(testy$x ~ testy$z1 + testy$z2 - 1)
all.equal(r^2, summary(bar)$r.squared)
first.term <- n * log(1 - r^2)
overid <- ncol(cbind(testy$z1, testy$z2)) - ncol(testy$x)
CC.BIC <- first.term + overid * log(n)
CC.AIC <- first.term + overid * 2
CC.HQ <- first.term + overid * 2.01 * log(log(n))
foo <- matrix(c(CC.BIC, CC.AIC, CC.HQ), 3, 1)
foo
baz
all.equal(foo, baz)

