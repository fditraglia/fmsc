#October 16th, 2011

#This script creates tables summarizing the simulation results

#---------------------------------------------#
#                   MACROS                    #
#---------------------------------------------#

sim.results.dir <- './Tables for Paper/Simulation Summary Tables'
infile.prefix <- 'sim_results_'
outfilename <- 'summary_tables_raw.tex'

#---------------------------------------------#
#                 END MACROS                  #
#---------------------------------------------#


#Package to output tables to LaTeX
library(Hmisc)

#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/"

setwd(root.dir)
setwd(sim.results.dir)

#Read simulation results
results.500 <- read.csv(paste(infile.prefix, '500.csv', sep = ''), header = TRUE)
results.100 <- read.csv(paste(infile.prefix, '100.csv', sep = ''), header = TRUE)
results.50 <- read.csv(paste(infile.prefix, '50.csv', sep = ''), header = TRUE)



r <- results_500$endog.w
g <- results_500$relevance.w


#Function that converts from long format dataframe (each combination of r and g gets its own row) to a wide format table (r varies over columns, g over rows).
make.table <- function(input.data){

  out <- data.frame(r, g, input.data)
  names(out) <- c('r', 'g', 'input.data')
  out <- reshape(out, v.names='input.data', idvar = 'g', timevar='r', direction="wide")
  names(out) <- c('g', paste('r=', unique(r), sep = '')) 
  return(out)
}#END make.table


#Function to collect desired quantities from a given results dataset
summary.results <- function(input.results){
  
  attach(input.results)
  
    avg.RMSE.valid <- mean(RMSE1)
    avg.RMSE.full <- mean(RMSE)
    avg.RMSE.FMSC <- mean(FMSC.RMSE)
    avg.RMSE.BIC <- mean(GMM.BIC.RMSE)
    avg.RMSE.HQ <- mean(GMM.HQ.RMSE)
    avg.RMSE.AIC <- mean(GMM.AIC.RMSE)
    avg.RMSE.J90 <- mean(J.90.RMSE)
    avg.RMSE.J95 <- mean(J.95.RMSE)
    avg.RMSE.CC.MSC.BIC <- mean(CC.MSC.BIC.RMSE)
    avg.RMSE.CC.MSC.HQ <- mean(CC.MSC.HQ.RMSE)
    avg.RMSE.CC.MSC.AIC <- mean(CC.MSC.AIC.RMSE)
    
    worst.RMSE.valid <- max(RMSE1)
    worst.RMSE.full <- max(RMSE)
    worst.RMSE.FMSC <- max(FMSC.RMSE)
    worst.RMSE.BIC <- max(GMM.BIC.RMSE)
    worst.RMSE.HQ <- max(GMM.HQ.RMSE)
    worst.RMSE.AIC <- max(GMM.AIC.RMSE)
    worst.RMSE.J90 <- max(J.90.RMSE)
    worst.RMSE.J95 <- max(J.95.RMSE)
    worst.RMSE.CC.MSC.BIC <- max(CC.MSC.BIC.RMSE)
    worst.RMSE.CC.MSC.HQ <- max(CC.MSC.HQ.RMSE)
    worst.RMSE.CC.MSC.AIC <- max(CC.MSC.AIC.RMSE)
    
    avg.CORRECT.FMSC <- mean(FMSC.CORRECT)
    avg.CORRECT.BIC <- mean(GMM.BIC.CORRECT)
    avg.CORRECT.HQ <- mean(GMM.HQ.CORRECT)
    avg.CORRECT.AIC <- mean(GMM.AIC.CORRECT)
    avg.CORRECT.J90 <- mean(J.90.CORRECT)
    avg.CORRECT.J95 <- mean(J.95.CORRECT)
    avg.CORRECT.CC.MSC.BIC <- mean(CC.MSC.BIC.CORRECT)
    avg.CORRECT.CC.MSC.HQ <- mean(CC.MSC.HQ.CORRECT)
    avg.CORRECT.CC.MSC.AIC <- mean(CC.MSC.AIC.CORRECT)
    
  
  detach(input.results)
  
  out <- rbind(avg.RMSE.valid, 
    avg.RMSE.full, 
    avg.RMSE.FMSC,
    avg.RMSE.BIC,
    avg.RMSE.HQ,
    avg.RMSE.AIC, 
    avg.RMSE.J90,
    avg.RMSE.J95, 
    avg.RMSE.CC.MSC.BIC, 
    avg.RMSE.CC.MSC.HQ, 
    avg.RMSE.CC.MSC.AIC, 
    worst.RMSE.valid, 
    worst.RMSE.full, 
    worst.RMSE.FMSC, 
    worst.RMSE.BIC,
    worst.RMSE.HQ, 
    worst.RMSE.AIC, 
    worst.RMSE.J90, 
    worst.RMSE.J95, 
    worst.RMSE.CC.MSC.BIC, 
    worst.RMSE.CC.MSC.HQ, 
    worst.RMSE.CC.MSC.AIC, 
    avg.CORRECT.FMSC, 
    avg.CORRECT.BIC, 
    avg.CORRECT.HQ, 
    avg.CORRECT.AIC, 
    avg.CORRECT.J90, 
    avg.CORRECT.J95, 
    avg.CORRECT.CC.MSC.BIC, 
    avg.CORRECT.CC.MSC.HQ, 
    avg.CORRECT.CC.MSC.AIC)
    
    return(round(out,2))  
}#END summary.results


N.500 <- data.frame(N500 = summary.results(results.500))
N.100 <- data.frame(N100 = summary.results(results.100))
N.50 <- data.frame(N50 = summary.results(results.50))

out.table <- cbind(N.50, N.100, N.500)



#Assign to a bogus variable to prevent the function from trying to load LaTex itself and compile the tables it creates
table <- latex(out.table, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = row.names(out.table), append = FALSE)


#Clean up
rm(list = ls())
