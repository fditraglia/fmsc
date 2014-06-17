library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(tikzDevice)

setwd("~/fmsc/OLSvsIV")
sourceCpp("simulation_functions_OLSvsIV.cpp")

set.seed(1928)