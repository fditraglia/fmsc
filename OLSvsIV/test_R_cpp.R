library(Rcpp)
library(RcppArmadillo)
library(parallel)


setwd("~/fmsc/OLSvsIV")
sourceCpp("simulation_functions_OLSvsIV.cpp")
source("simulation_functions_OLSvsIV.R")


set.seed(2)
foo <- simple.sim(0.3, 0.2, 250)
set.seed(2)
bar <- simple_sim_cpp(0.3, 0.2, 250)
all.equal(foo, bar)



library(microbenchmark)
#microbenchmark(mse.compare.cpp(0.3, 0.2, 250), mse.compare(0.3, 0.2, 250))

set.seed(3728)
barR <- mse.compare(0.3, 0.2, 250)
set.seed(3728)
barCpp <- mse_compare_cpp(1, 0.3 * rep(1, 3), matrix(c(1, 0.2, 0.2, 1), 2, 2), diag(rep(1, 3)), 250, 10000)
  
sum(abs(barR - barCpp))


r.seq <- seq(0, 0.2, 0.01)

system.time(do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse.compare(p = 0.3, r, n = 250)}, mc.set.seed = FALSE)))

system.time(do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse.compare.cpp(p = 0.3, r, n = 250)}, mc.set.seed = FALSE)))


set.seed(4938)
fooR <- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse.compare(p = 0.3, r, n = 250)}, mc.set.seed = FALSE))#mclapply outputs a list of vectors. Combine them into a matrix.

set.seed(4938)
fooCpp <- do.call(rbind, mclapply(X = r.seq, FUN = function(r){mse.compare.cpp(p = 0.3, r, n = 250)}, mc.set.seed = FALSE))

all.equal(fooR, fooCpp)


matplot(r.seq, apply(fooCpp, 2, sqrt))
#, col = c('black', 'red', 'blue'), xlab = 'Cor(e,v)', ylab = 'RMSE', type =  'l', lty = 1)
# legend("topleft", c("FMSC", "OLS", "IV"), fill = c("black", "red", "blue"))

