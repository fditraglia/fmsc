library(Rcpp)
library(RcppArmadillo)
library(parallel)


setwd("~/fmsc/OLSvsIV")
sourceCpp("simulation_functions_OLSvsIV.cpp")
source("simulation_functions_OLSvsIV.R")


r.seq <- seq(0, 0.2, 0.01)

set.seed(4938)
fooR <- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse.compare(0.3, r, 250, 10000)}, mc.set.seed = FALSE))#mclapply outputs a list of vectors. Combine them into a matrix.

set.seed(4938)
fooCpp <- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse_compare_default_cpp(0.3, r, 250, 10000)}, mc.set.seed = FALSE))

all.equal(fooR, fooCpp)


matplot(r.seq, apply(fooCpp, 2, sqrt), xlab = 'Cor(e,v)', ylab = 'RMSE',
        type = 'l', lty = 1, lwd = 2)
legend("bottomright", c("OLS", "2SLS", "FMSC", "AVG", "DHW90", "DHW95"), fill = 1:6)

