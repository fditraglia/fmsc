library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(tikzDevice)

sourceCpp("functions.cpp")

set.seed(1928)

#Calculate RMSE results for the simulation
#and store them in ./Results/rmse_results.Rdata
source("mse_calculations.R")

#Open ./Results/rmse_results.Rdata
#plot the results 
#and export in tikz format to ./Results/
source("rmse_plots.R")

#Calculate Confidence Interval results for simulation
#and store them in ./Results/CI_results.Rdata
source("CI_calculations.R")

#Load CI_results.Rdata and 
#format into two lists of tables
#One for coverage probability
#and another for median CI width relative
#to TSLS
source("CI_Rtables.R")

#Convert the lists of CI tables
#to LaTeX and output to ./Results/
source("CI_TeXtables.R")
