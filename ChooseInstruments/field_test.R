setwd("~/fmsc/ChooseInstruments")

library(Rcpp)
library(RcppArmadillo)

sourceCpp("field_test.cpp")

mat0 <- diag(1, 3)
mat1 <- matrix(rep(3, 4), 2, 2)
mat2 <- matrix(rnorm(16), 4, 4)

field_test(mat0, mat1, mat2, 0)
field_test(mat0, mat1, mat2, 1)
field_test(mat0, mat1, mat2, 2)

