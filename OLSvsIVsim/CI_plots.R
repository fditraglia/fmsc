library(Rcpp)
library(RcppArmadillo)

sourceCpp("functions.cpp")


set.seed(8273)


testy <- test_CIs_cpp(0.1, 0.1, 250, 1000)
