#October 5th, 2011

#This script carries out the simulation study with a sample size of 500 to see how well the naive confidence interval correction for FMSC moment selection performs

#----------------------------------------------------------#
#                    CONSTANT MACROS                       #
#----------------------------------------------------------#

N <- 500 #sample size of each simulation
replications <- 10000 #Number of replications
outfilename <- "sim_results_FMSC_naive_correct_N_500.csv"
output.dir <- './Test Run'

#----------------------------------------------------------#
#                   END CONSTANT MACROS                    #
#----------------------------------------------------------#

#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Correct Confidence FMSC Naive"

output.dir <- './N = 500'
outfile.name <- "FMSC_naive_conf_N_500"

setwd(root.dir)
setwd("./Functions")

source("reg2SLS.R")
source("correct_FMSC_interval_naive.R")
source("linear_simulation.R")
source("simulation_study_FMSC_conf_naive.R")
library(mvtnorm)


#Different values for endogeneity of w
rho <- seq(from = 0, to = 0.4, by = 0.05)


#Values for relevance of w
g <- seq(from = 0, to = 1.3, by = 0.1)


#Create dataframe listing all pairs of simulation parameters
rho.column <- rep(rho, each = length(g))
g.column <- rep(g, times = length(rho))
parameters <- data.frame(rho = rho.column, g = g.column)


#Initialize Dataframe to store results
results <- matrix(NA, nrow = nrow(parameters), ncol = 3)
results <- as.data.frame(results)


names(results) <- c('relevance.w', 'endog.w', 'naive.FMSC.correct.cover')


#Starting time
start.time <- proc.time()

#Set seed for replication
set.seed(1100)

for(i in 1:nrow(parameters)){
  
	results[i,] <- simulation.study.FMSC.conf.naive(N.reps = replications, sample.size = N, relevance.w = parameters$g[i], endog.w = parameters$rho[i])
	
	
}#END for i

#How long did all of that take?
proc.time() - start.time

#Store results
setwd(root.dir)
setwd(output.dir)

write.csv(results, outfilename, row.names = FALSE)


#Clean up
#rm(list = ls())

