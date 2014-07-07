#September 27th, 2011

#This script carries out the simulation study with a sample size of 500, incorporating my latest changes, adding moment smooth FMSC moment averaging and smoothed GMM-AIC moment averaging

#----------------------------------------------------------#
#                    CONSTANT MACROS                       #
#----------------------------------------------------------#

N <- 500 #sample size of each simulation
replications <- 10000 #Number of replications
outfilename <- "sim_results_avg_500.csv"
output.dir <- './N = 500'

#----------------------------------------------------------#
#                   END CONSTANT MACROS                    #
#----------------------------------------------------------#




#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Moment Averaging"

setwd(root.dir)
setwd("./Functions")

source("reg2SLS.R")
source("compute_criteria.R")
source("linear_simulation.R")
source("simulation_study_averaging_only.R")
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
results <- matrix(NA, nrow = nrow(parameters), ncol = 7)
results <- as.data.frame(results)


names(results) <- c('relevance.w', 'endog.w', 'RMSE.BIC.avg', 'RMSE.AIC.avg', 'RMSE.HQ.avg', 'RMSE.FMSC.avg', 'RMSE.FMSC.avg.k')


#Starting time
start.time <- proc.time()

#Set seed for replication
set.seed(1100)

for(i in 1:nrow(parameters)){
	
	results[i,] <- simulation.study.avg(N.reps = replications, sample.size = N, relevance.w = parameters$g[i], endog.w = parameters$rho[i])
	
	
}#END for i

#How long did all of that take?
proc.time() - start.time

#Store results
setwd(root.dir)
setwd(output.dir)

write.csv(results, outfilename, row.names = FALSE)


#Clean up
#rm(list = ls())