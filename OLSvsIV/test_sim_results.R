library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(tikzDevice)

setwd("~/fmsc/OLSvsIV")
sourceCpp("simulation_functions_OLSvsIV.cpp")

set.seed(1928)
system.time(foo1 <- mse_compare_default_cpp(0.3, 0.2, 250, 10000))
foo1

Vz <- diag(1/3, 3, 3)
p <- 0.3
pvec <- p * c(1, 1, 1)
r <- 0.2
Ve <- matrix(c(1, r * sqrt(1 - p^2),
               r * sqrt(1 - p^2), 1 - p^2), 2, 2) 
set.seed(1928)
system.time(foo2 <- mse_compare_cpp(0.5, pvec, Ve, Vz, 250, 10000))
foo2

all.equal(foo1, foo2)

r.seq <- seq(0, 0.6, 0.01)

set.seed(4938)
system.time(bar <- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse_compare_default_cpp(0.3, r, 250, 20000)}, mc.set.seed = FALSE)))



matplot(r.seq, apply(bar, 2, sqrt), xlab = 'Cor(e,v)', ylab = 'RMSE',
        type = 'l', lty = 1, lwd = 2)
legend("topleft", c("OLS", "2SLS", "FMSC", "AVG", "DHW90", "DHW95"), fill = 1:6)

#Plot expressed in percentage points relative to the oracle
oracle <- pmin(bar[,1], bar[,2])
rel.rmse <- (bar - oracle)/oracle * 100
matplot(r.seq, rel.rmse[,-c(1:2)], xlab = 'Cor(e,v)', ylab = 'RMSE Relative to Oracle (%)', type = 'l', lty = 1:4, lwd = 4, col = c("black", "blue", "red", "orange"))
legend("topright", c("FMSC", "AVG", "DHW90", "DHW95"), fill = c("black", "blue", "red", "orange"), lty = 1:4, lwd = 3)

