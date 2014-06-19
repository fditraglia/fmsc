library(Rcpp)
library(RcppArmadillo)

sourceCpp("functions.cpp")


set.seed(8273)
testy <- CIs_compare_default_cpp(0.2, 0.2, 500, 1000)
