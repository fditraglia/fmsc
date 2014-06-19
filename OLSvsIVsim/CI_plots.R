library(Rcpp)
library(RcppArmadillo)

sourceCpp("functions.cpp")


set.seed(8273)
testy <- CIs_compare_default_cpp(0.2, 0.2, 1000, 1000)
