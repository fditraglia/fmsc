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


system.time(testy <- test_CIs_cpp(0.3, 0.15, 250, 10000))
sum(testy[,1] < 1 & testy[,2] > 1)/nrow(testy)
median(testy[,2] - testy[,1])


r.seq <- seq(0, 0.2, 0.01)

set.seed(4938)
system.time(bar <- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse_compare_default_cpp(0.3, r, 250, 10000)}, mc.set.seed = FALSE)))


matplot(r.seq, apply(bar, 2, sqrt), xlab = 'Cor(e,v)', ylab = 'RMSE',
        type = 'l', lty = 1, lwd = 2)
legend("bottomright", c("OLS", "2SLS", "FMSC", "AVG", "DHW90", "DHW95"), fill = 1:6)
