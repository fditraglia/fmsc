library(Rcpp)
library(RcppArmadillo)
library(parallel)

setwd("~/fmsc/OLSvsIV")
sourceCpp("simulation_functions_OLSvsIV.cpp")

set.seed(9281)
system.time(foo <- mse_compare_default_cpp(0.3, 0.2, 250, 10000))

set.seed(438)
system.time(bar <- mse_compare_default_cpp(0.3, 0.2, 250, 10000))

sum(abs(foo - bar))


system.time(testy <- test_CIs_cpp(0.1, 0.1, 250, 10000))
sum(testy[,1] < 1 & testy[,2] > 1)/nrow(testy)
median(testy[,2] - testy[,1])


r.seq <- seq(0, 0.6, 0.005)

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

