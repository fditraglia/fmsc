library(Rcpp)
library(RcppArmadillo)
library(parallel)
#library(tikzDevice)

setwd("~/fmsc/OLSvsIV")
sourceCpp("simulation_functions_OLSvsIV.cpp")

# set.seed(1928)
# system.time(foo1 <- mse_compare_default_cpp(0.3, 0.2, 250, 10000))
# foo1
# 
# Vz <- diag(1/3, 3, 3)
# p <- 0.3
# pvec <- p * c(1, 1, 1)
# r <- 0.2
# Ve <- matrix(c(1, r * sqrt(1 - p^2),
#                r * sqrt(1 - p^2), 1 - p^2), 2, 2) 
# set.seed(1928)
# system.time(foo2 <- mse_compare_cpp(0.5, pvec, Ve, Vz, 250, 10000))
# foo2
# 
# all.equal(foo1, foo2)
# 
# r.seq <- seq(0, 0.6, 0.01)
# 
# set.seed(4938)
# system.time(bar <- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse_compare_default_cpp(0.3, r, 250, 20000)}, mc.set.seed = FALSE)))
# 
# 
# 
# matplot(r.seq, apply(bar, 2, sqrt), xlab = 'Cor(e,v)', ylab = 'RMSE',
#         type = 'l', lty = 1, lwd = 2)
# legend("topleft", c("OLS", "2SLS", "FMSC", "AVG", "DHW90", "DHW95"), fill = 1:6)
# 
# #Plot expressed in percentage points relative to the oracle
# oracle <- pmin(bar[,1], bar[,2])
# rel.rmse <- (bar - oracle)/oracle * 100
# matplot(r.seq, rel.rmse[,-c(1:2)], xlab = 'Cor(e,v)', ylab = 'RMSE Relative to Oracle (%)', type = 'l', lty = 1:4, lwd = 4, col = c("black", "blue", "red", "orange"))
# legend("topright", c("FMSC", "AVG", "DHW90", "DHW95"), fill = c("black", "blue", "red", "orange"), lty = 1:4, lwd = 3)

set.seed(2819)

n.reps <- 1000
rho.fine <- seq(0, 0.6, 0.1)
pi.fine <- seq(0.1, 0.6, 0.1)
rho.coarse <- seq(0.1, 0.2, 0.1)
pi.coarse <- seq(0.1, 0.2, 0.1)

sample.size.grid <- c(50, 100, 250)

#Run the simulation over all values in rho.fine for a fixed value of pi and sample size
fix.pi <- function(pi.value, sample.size){
  out <- lapply(rho.fine, function(r) mse_compare_default_cpp(pi.value, r, sample.size, n.reps))
  out <- do.call(rbind, out) #Convert list to matrix
  nRows <- length(rho.fine)
  cbind(n = rep(sample.size, nRows), p = rep(pi.value, nRows), r = rho.fine, out)
}

#Run the simulation over all values in pi.fine for a fixed value of rho and sample size
fix.rho <- function(rho.value, sample.size){
  out <- lapply(pi.fine, function(p) mse_compare_default_cpp(p, rho.value, sample.size, n.reps))
  out <- do.call(rbind, out) #Convert list to matrix
  nRows <- length(pi.fine)
  cbind(n = rep(sample.size, nRows), p = pi.fine, r = rep(rho.value, nRows), out)
}

#Run fix.pi over all values in pi.coarse for a fixed sample size
panels.coarse.pi <- function(sample.size){
  lapply(pi.coarse, function(p) fix.pi(p, sample.size))
} 

#Run fix.rho over all values in rho.coarse for a fixed sample size
panels.coarse.rho <- function(sample.size){
  lapply(rho.coarse, function(r) fix.rho(r, sample.size))
}

results.coarse.rho <- lapply(sample.size.grid, panels.coarse.rho)
results.coarse.pi <- lapply(sample.size.grid, panels.coarse.pi)