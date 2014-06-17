library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(tikzDevice)

sourceCpp("functions.cpp")

set.seed(1928)

#Calculate RMSE results for the simulation
#and store them in rmse_results.Rdata
source("rmse_calculations.R")

#Open rmse_results.Rdata
#plot the results 
#and export in tikz format
#source("rmse_plots.R")

#source("CI_calculations.R")

#source("CI_plots.R")