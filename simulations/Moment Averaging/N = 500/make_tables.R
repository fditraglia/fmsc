#September 30th, 2011

#This script creates tables summarizing the simulation results

#---------------------------------------------#
#                   MACROS                    #
#---------------------------------------------#

#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Moment Averaging"

sim.results.dir <- './N = 500'
infilename <- 'sim_results_avg_500.csv'
outfilename <- 'sim_tables_avg_raw_500.tex'

#---------------------------------------------#
#                 END MACROS                  #
#---------------------------------------------#


#Package to output tables to LaTeX
library(Hmisc)


setwd(root.dir)
setwd(sim.results.dir)

#Read simulation results
results <- read.csv(infilename, header = TRUE)

#Results from moment selection simulation, for comparison with moment averaging
select.results <- read.csv('sim_results_500.csv', header = TRUE)

attach(results)

r <- endog.w
g <- relevance.w


#Function that converts from long format dataframe (each combination of r and g gets its own row) to a wide format table (r varies over columns, g over rows).
make.table <- function(input.data){

  out <- data.frame(r, g, input.data)
  names(out) <- c('r', 'g', 'input.data')
  out <- reshape(out, v.names='input.data', idvar = 'g', timevar='r', direction="wide")
  names(out) <- c('g', paste('r=', unique(r), sep = '')) 
  return(out)
}#END make.table


#RMSE of different procedures
rmse.FMSC.avg <- make.table(round(RMSE.FMSC.avg, 2))
rmse.FMSC.avg.k <- make.table(round(RMSE.FMSC.avg.k, 2))
rmse.BIC.avg <- make.table(round(RMSE.BIC.avg, 2))
rmse.AIC.avg <- make.table(round(RMSE.AIC.avg, 2))
rmse.HQ.avg <- make.table(round(RMSE.HQ.avg, 2))



#Comparisons of RMSE against FMSC
#How much lower is the realized MSE of the FMSC than that of other procedures?
rmse.vs.gmm.bic.avg <- make.table(round(RMSE.FMSC.avg - RMSE.BIC.avg, 2))
rmse.vs.gmm.aic.avg <- make.table(round(RMSE.FMSC.avg - RMSE.AIC.avg, 2))
rmse.vs.gmm.hq.avg <- make.table(round(RMSE.FMSC.avg - RMSE.HQ.avg, 2))

#How much lower is the realized MSE of the FMSC than that of other procedures? Now with k = 1/100
rmse.k.vs.gmm.bic.avg <- make.table(round(RMSE.FMSC.avg.k - RMSE.BIC.avg, 2))
rmse.k.vs.gmm.aic.avg <- make.table(round(RMSE.FMSC.avg.k - RMSE.AIC.avg, 2))
rmse.k.vs.gmm.hq.avg <- make.table(round(RMSE.FMSC.avg.k - RMSE.HQ.avg, 2))


detach(results)


#How does RMSE of selection procedure compare to equivalent averaging procedure?
decrease.RMSE.FMSC.avg <- make.table(round(results$RMSE.FMSC.avg.k - select.results$FMSC.RMSE, 2))

decrease.RMSE.BIC.avg <- make.table(round(results$RMSE.BIC.avg - select.results$GMM.BIC.RMSE, 2))

decrease.RMSE.HQ.avg <- make.table(round(results$RMSE.HQ.avg - select.results$GMM.HQ.RMSE, 2))

decrease.RMSE.AIC.avg <- make.table(round(results$RMSE.AIC.avg - select.results$GMM.AIC.RMSE, 2))

avg.RMSE.FMSC <- c(mean(results$RMSE.FMSC.avg.k), mean(select.results$FMSC.RMSE))
avg.RMSE.BIC <- c(mean(results$RMSE.BIC.avg), mean(select.results$GMM.BIC.RMSE))
avg.RMSE.HQ <- c(mean(results$RMSE.HQ.avg), mean(select.results$GMM.HQ.RMSE))
avg.RMSE.AIC <- c(mean(results$RMSE.AIC.avg), mean(select.results$GMM.AIC.RMSE))

avg.RMSE <- rbind(avg.RMSE.FMSC, avg.RMSE.BIC, avg.RMSE.HQ, avg.RMSE.AIC)
colnames(avg.RMSE) <- c('Average', 'Select')

worst.RMSE.FMSC <- c(max(results$RMSE.FMSC.avg.k), max(select.results$FMSC.RMSE))
worst.RMSE.BIC <- c(max(results$RMSE.BIC.avg), max(select.results$GMM.BIC.RMSE))
worst.RMSE.HQ <- c(max(results$RMSE.HQ.avg), max(select.results$GMM.HQ.RMSE))
worst.RMSE.AIC <- c(max(results$RMSE.AIC.avg), max(select.results$GMM.AIC.RMSE))

worst.RMSE <- rbind(worst.RMSE.FMSC, worst.RMSE.BIC, worst.RMSE.HQ, worst.RMSE.AIC)
colnames(worst.RMSE) <- c('Average', 'Select')

RMSE.average.vs.select <- round(rbind(avg.RMSE, worst.RMSE),2)



#Output all of the tables in LaTeX tabular format to a table.

#Assign to a bogus variable to prevent the function from trying to load LaTex itself and compile the tables it creates
table <- latex(rmse.FMSC.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = FALSE)

#Append the remaining tables
table <- latex(rmse.FMSC.avg.k, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(rmse.BIC.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.AIC.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.HQ.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.gmm.bic.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.gmm.aic.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.gmm.hq.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.k.vs.gmm.bic.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.k.vs.gmm.aic.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.k.vs.gmm.hq.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)


table <- latex(RMSE.average.vs.select, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

table <- latex(decrease.RMSE.FMSC.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(decrease.RMSE.BIC.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(decrease.RMSE.HQ.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(decrease.RMSE.AIC.avg, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)


#Clean up
rm(list = ls())
