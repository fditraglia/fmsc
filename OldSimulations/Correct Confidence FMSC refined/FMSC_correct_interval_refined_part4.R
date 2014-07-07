#October 11th, 2011

#This script carries out the simulation study with a sample size of 500 to see how well the refined confidence interval correction for FMSC moment selection performs

#Since it would otherwise take too long, we break the parameter space up into different pieces and give each to a single core. Each point on the simulation grid takes a little under 5 hours on a single core.


#----------------------------------------------------------#
#                 WHICH PART IS THIS?                      #
#----------------------------------------------------------#

part <- 4

#----------------------------------------------------------#
#                 WHICH PART IS THIS?                      #
#----------------------------------------------------------#





#----------------------------------------------------------#
#                    CONSTANT MACROS                       #
#----------------------------------------------------------#

N <- 500 #sample size of each simulation
replications <- 10000 #Number of replications


#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Correct Confidence FMSC Refined"

#For the F drive
#root.dir <- "F:/Correct Confidence FMSC Refined"

output.dir <- root.dir
outfile.name <- paste("FMSC_refined_conf_test_N_500", "_part", part, ".csv", sep = '')

#----------------------------------------------------------#
#                   END CONSTANT MACROS                    #
#----------------------------------------------------------#



setwd(root.dir)
setwd("./Functions")

source("reg2SLS.R")
source("correct_FMSC_interval_refined.R")
source("linear_simulation.R")
source("simulation_study_FMSC_conf_refined.R")
library(mvtnorm)


#Different values for endogeneity of w
rho <- seq(from = 0, to = 0.4, by = 0.1)
#rho <- seq(from = 0, to = 0.4, by = 0.05)
#rho <- 0.2


#Values for relevance of w
g <- seq(from = 0, to = 1.3, by = 0.2)
#g <- seq(from = 0, to = 1.3, by = 0.1)
#g <- 0.4 


#Create dataframe listing all pairs of simulation parameters
rho.column <- rep(rho, each = length(g))
g.column <- rep(g, times = length(rho))
full.parameters <- data.frame(rho = rho.column, g = g.column)

#Break the parameter space into five blocks and create a list from the five blocks
part1 <- full.parameters[1:7,]
part2 <- full.parameters[8:14,]
part3 <- full.parameters[15:21,]
part4 <- full.parameters[22:28,]
part5 <- full.parameters[29:35,]
parts <- list(part1, part2, part3, part4, part5)

#Which part to use for this simulation?
parameters <- parts[[part]]


#Initialize Dataframe to store results
results <- matrix(NA, nrow = nrow(parameters), ncol = 3)
results <- as.data.frame(results)


names(results) <- c('relevance.w', 'endog.w', 'refined.FMSC.correct.cover')


#Starting time
start.time <- proc.time()

#Set seed for replication
set.seed(1100)

for(i in 1:nrow(parameters)){
  
	results[i,] <- simulation.study.FMSC.conf.refined(N.reps = replications, sample.size = N, relevance.w = parameters$g[i], endog.w = parameters$rho[i])
	
	
}#END for i

#How long did all of that take?
proc.time() - start.time

#Store results
setwd(root.dir)
setwd(output.dir)

write.csv(results, outfile.name, row.names = FALSE)


#Clean up
#rm(list = ls())

