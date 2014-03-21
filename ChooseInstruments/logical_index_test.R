#This script tests the use of logical indices to subset
#the columns of a matrix in armadillo. The documentation
#on this feature is a little confusing so I wanted to
#construct a simple example.

setwd("~/fmsc/ChooseInstruments")

library(Rcpp)
library(RcppArmadillo)

sourceCpp("logical_index_test.cpp")

M <- matrix(rep(1:5, each = 5), 5, 5)

col_indices <- cbind(c(TRUE, FALSE, FALSE, FALSE, TRUE),
                     c(FALSE, TRUE, TRUE, TRUE, FALSE),
                     c(TRUE, TRUE, TRUE, TRUE, TRUE),
                     c(TRUE, FALSE, FALSE, FALSE, FALSE))

M[,col_indices[,1]]
M[,col_indices[,2]]
M[,col_indices[,3]]
M[,col_indices[,4]]

#Aha! 
logical_subset(M, c(TRUE, FALSE, TRUE, FALSE, TRUE))
