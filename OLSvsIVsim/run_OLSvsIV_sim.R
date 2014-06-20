library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(tikzDevice)

sourceCpp("functions.cpp")

set.seed(1928)

#Calculate RMSE results for the simulation
#and store them in rmse_results.Rdata
#source("mse_calculations.R")

#Open rmse_results.Rdata
#plot the results 
#and export in tikz format
#source("rmse_plots.R")

#Calculate Confidence Interval results for simulation
#and store them in CI_results.Rdata
#source("CI_calculations.R")

#Load CI_results.Rdata 
#plot results and export
#in tikz format
source("CI_plots.R")