#September 23rd, 2011

#This script carries out the simulation study with a sample size of 100, incorporating my latest changes.

#----------------------------------------------------------#
#                    CONSTANT MACROS                       #
#----------------------------------------------------------#

N <- 50 #sample size of each simulation
replications <- 10000 #Number of replications
outfilename <- "sim_results_50.csv"
output.dir <- './N = 50'

#----------------------------------------------------------#
#                   END CONSTANT MACROS                    #
#----------------------------------------------------------#




#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Weak Instruments"

setwd(root.dir)
setwd("./Functions")

source("reg2SLS.R")
source("compute_criteria.R")
source("linear_simulation.R")
source("simulation_study.R")
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
results <- matrix(NA, nrow = nrow(parameters), ncol = 54)
results <- as.data.frame(results)

criterion.names <- c('FMSC', 'GMM.BIC', 'GMM.AIC', 'GMM.HQ', 'CC.BIC', 'CC.AIC', 'CC.HQ', 'CC.MSC.BIC', 'CC.MSC.AIC', 'CC.MSC.HQ', 'J.90', 'J.95')

summary.names <- c('BIAS', 'RMSE', 'COVERAGE', 'CORRECT')
names2 <- paste(rep(criterion.names, each = length(summary.names)), rep(summary.names, length(criterion.names)), sep = '.')

names(results) <- c('relevance.w', 'endog.w', 'RMSE', 'RMSE1', 'COVERAGE', 'COVERAGE1', names2)


#Starting time
start.time <- proc.time()

#Set seed for replication
set.seed(1100)

for(i in 1:nrow(parameters)){
	
	results[i,] <- simulation.study(N.reps = replications, sample.size = N, relevance.w = parameters$g[i], endog.w = parameters$rho[i])
	
	
}#END for i

#How long did all of that take?
proc.time() - start.time

#Store results
setwd(root.dir)
setwd(output.dir)

write.csv(results, outfilename, row.names = FALSE)

#true.rmse <- results[,1:4]
#advantage.w <- round((true.rmse$RMSE - true.rmse$RMSE1), 2) 
#true.rmse <- cbind(true.rmse, advantage.w)

#Clean up
#rm(list = ls())