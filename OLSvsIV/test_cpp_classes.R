library(Rcpp)
library(RcppArmadillo)

setwd("~/fmsc/OLSvsIV")
sourceCpp("simulation_functions_OLSvsIV_classes.cpp")
sourceCpp("simulation_functions_OLSvsIV.cpp")

set.seed(438)
system.time(foo <- mse_compare_default_cpp(0.3, 0.2, 250, 10000))


set.seed(438)
system.time(bar <- mse_compare_default_classes(0.3, 0.2, 250, 10000))

all.equal(foo, bar)


r.seq <- seq(0, 0.2, 0.01)

set.seed(4938)
foo<- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse_compare_default_cpp(0.3, r, 250, 10000)}, mc.set.seed = FALSE))

set.seed(4938)
bar <- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse_compare_default_classes(0.3, r, 250, 10000)}, mc.set.seed = FALSE))

all.equal(foo, bar)

matplot(r.seq, apply(bar, 2, sqrt), xlab = 'Cor(e,v)', ylab = 'RMSE',
        type = 'l', lty = 1, lwd = 2)
legend("bottomright", c("OLS", "2SLS", "FMSC", "AVG", "DHW90", "DHW95"), fill = 1:6)
