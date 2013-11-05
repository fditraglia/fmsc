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
results_500 <- read.csv(paste(infile.prefix, '500.csv', sep = ''), header = TRUE)
results_100 <- read.csv(paste(infile.prefix, '100.csv', sep = ''), header = TRUE)
results_50 <- read.csv(paste(infile.prefix, '50.csv', sep = ''), header = TRUE)



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


#RMSE of different procedures
rmse.valid <- make.table(round(RMSE1, 2))
rmse.full <- make.table(round(RMSE, 2))

rmse.FMSC <- make.table(round(FMSC.RMSE, 2))
rmse.BIC <- make.table(round(GMM.BIC.RMSE, 2))
rmse.AIC <- make.table(round(GMM.AIC.RMSE, 2))
rmse.HQ <- make.table(round(GMM.HQ.RMSE, 2))
rmse.CC.BIC <- make.table(round(CC.BIC.RMSE, 2))
rmse.CC.AIC <- make.table(round(CC.AIC.RMSE, 2))
rmse.CC.HQ <- make.table(round(CC.HQ.RMSE, 2))
rmse.J.90 <- make.table(round(J.90.RMSE, 2))
rmse.J.95 <- make.table(round(J.95.RMSE, 2))
rmse.CC.MSC.BIC <- make.table(round(CC.MSC.BIC.RMSE, 2))
rmse.CC.MSC.AIC <- make.table(round(CC.MSC.AIC.RMSE, 2))
rmse.CC.MSC.HQ <- make.table(round(CC.MSC.HQ.RMSE, 2))

#Comparisons of RMSE
#How much lower is true RMSE when w is included in the instrument set?
w.advantage <- make.table(round(RMSE - RMSE1, 2))

#Comparisons of RMSE against FMSC
#How much lower is the realized MSE of the FMSC than that of other procedures?
rmse.vs.valid <- make.table(round(FMSC.RMSE - RMSE1, 2))
rmse.vs.full <- make.table(round(FMSC.RMSE - RMSE, 2))
rmse.vs.J.90 <- make.table(round(FMSC.RMSE - J.90.RMSE, 2))
rmse.vs.J.95 <- make.table(round(FMSC.RMSE - J.95.RMSE, 2))
rmse.vs.gmm.bic <- make.table(round(FMSC.RMSE - GMM.BIC.RMSE, 2))
rmse.vs.gmm.aic <- make.table(round(FMSC.RMSE - GMM.AIC.RMSE, 2))
rmse.vs.gmm.hq <- make.table(round(FMSC.RMSE - GMM.HQ.RMSE, 2))
rmse.vs.cc.bic <- make.table(round(FMSC.RMSE - CC.BIC.RMSE, 2))
rmse.vs.cc.aic <- make.table(round(FMSC.RMSE - CC.AIC.RMSE, 2))
rmse.vs.cc.hq <- make.table(round(FMSC.RMSE - CC.HQ.RMSE, 2))
rmse.vs.cc.msc.bic <- make.table(round(FMSC.RMSE - CC.MSC.BIC.RMSE, 2))
rmse.vs.cc.msc.aic <- make.table(round(FMSC.RMSE - CC.MSC.AIC.RMSE, 2))
rmse.vs.cc.msc.hq <- make.table(round(FMSC.RMSE - CC.MSC.HQ.RMSE, 2))



#Correct Decision Rates of different procedures
correct.FMSC <- make.table(round(FMSC.CORRECT, 2)*100)
correct.BIC <- make.table(round(GMM.BIC.CORRECT, 2)*100)
correct.AIC <- make.table(round(GMM.AIC.CORRECT, 2)*100)
correct.HQ <- make.table(round(GMM.HQ.CORRECT, 2)*100)
correct.CC.BIC <- make.table(round(CC.BIC.CORRECT, 2)*100)
correct.CC.AIC <- make.table(round(CC.AIC.CORRECT, 2)*100)
correct.CC.HQ <- make.table(round(CC.HQ.CORRECT, 2)*100)
correct.J.90 <- make.table(round(J.90.CORRECT, 2)*100)
correct.J.95 <- make.table(round(J.95.CORRECT, 2)*100)
correct.CC.MSC.BIC <- make.table(round(CC.MSC.BIC.CORRECT, 2)*100)
correct.CC.MSC.AIC <- make.table(round(CC.MSC.AIC.CORRECT, 2)*100)
correct.CC.MSC.HQ <- make.table(round(CC.MSC.HQ.CORRECT, 2)*100)

#Comparisons of Correct Decision Rates to FMSC
correct.vs.J.90 <- make.table(round(FMSC.CORRECT - J.90.CORRECT, 2)*100)
correct.vs.J.95 <- make.table(round(FMSC.CORRECT - J.95.CORRECT, 2)*100)
correct.vs.gmm.bic <- make.table(round(FMSC.CORRECT - GMM.BIC.CORRECT, 2)*100)
correct.vs.gmm.aic <- make.table(round(FMSC.CORRECT - GMM.AIC.CORRECT, 2)*100)
correct.vs.gmm.hq <- make.table(round(FMSC.CORRECT - GMM.HQ.CORRECT, 2)*100)
correct.vs.cc.bic <- make.table(round(FMSC.CORRECT - CC.BIC.CORRECT, 2)*100)
correct.vs.cc.aic <- make.table(round(FMSC.CORRECT - CC.AIC.CORRECT, 2)*100)
correct.vs.cc.hq <- make.table(round(FMSC.CORRECT - CC.HQ.CORRECT, 2)*100)
correct.vs.cc.msc.bic <- make.table(round(FMSC.CORRECT - CC.MSC.BIC.CORRECT, 2)*100)
correct.vs.cc.msc.aic <- make.table(round(FMSC.CORRECT - CC.MSC.AIC.CORRECT, 2)*100)
correct.vs.cc.msc.hq<- make.table(round(FMSC.CORRECT - CC.MSC.HQ.CORRECT, 2)*100)


#Coverage of naive 95% confidence intervals around (biased) quasi-true value
cover.valid <- make.table(round(COVERAGE1, 2))
cover.full <- make.table(round(COVERAGE, 2))
cover.FMSC <- make.table(round(FMSC.COVERAGE, 2))
cover.BIC <- make.table(round(GMM.BIC.COVERAGE, 2))
cover.AIC <- make.table(round(GMM.AIC.COVERAGE, 2))
cover.HQ <- make.table(round(GMM.HQ.COVERAGE, 2))
cover.CC.BIC <- make.table(round(CC.BIC.COVERAGE, 2))
cover.CC.AIC <- make.table(round(CC.AIC.COVERAGE, 2))
cover.CC.HQ <- make.table(round(CC.HQ.COVERAGE, 2))
cover.J.90 <- make.table(round(J.90.COVERAGE, 2))
cover.J.95 <- make.table(round(J.95.COVERAGE, 2))
cover.CC.MSC.BIC <- make.table(round(CC.MSC.BIC.COVERAGE, 2))
cover.CC.MSC.AIC <- make.table(round(CC.MSC.AIC.COVERAGE, 2))
cover.CC.MSC.HQ <- make.table(round(CC.MSC.HQ.COVERAGE, 2))



detach(results)


#Output all of the tables in LaTeX tabular format to a table.

#Assign to a bogus variable to prevent the function from trying to load LaTex itself and compile the tables it creates
table <- latex(w.advantage, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = FALSE)

#Append the remaining tables
table <- latex(rmse.valid, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.full, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.FMSC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.BIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.AIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.HQ, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.CC.BIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.CC.AIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.CC.HQ, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.J.90, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.J.95, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.CC.MSC.BIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.CC.MSC.AIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.CC.MSC.HQ, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 



table <- latex(rmse.vs.valid, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.full, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.J.90, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.J.95, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.gmm.bic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.gmm.aic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.gmm.hq, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.cc.bic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.cc.aic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.cc.hq, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.cc.msc.bic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(rmse.vs.cc.msc.aic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(rmse.vs.cc.msc.hq, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 



table <- latex(correct.FMSC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.BIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.AIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.HQ, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.CC.BIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.CC.AIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.CC.HQ, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.J.90, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.J.95, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.CC.MSC.BIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.CC.MSC.AIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(correct.CC.MSC.HQ, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)



table <- latex(correct.vs.J.90, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.J.95, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.gmm.bic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.gmm.aic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.gmm.hq, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.cc.bic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.cc.aic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.cc.hq, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.cc.msc.bic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.cc.msc.aic, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(correct.vs.cc.msc.hq, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 



table <- latex(cover.valid, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.full, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.FMSC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.BIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.AIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.HQ, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.CC.BIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.CC.AIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.CC.HQ, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.J.90, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(cover.J.95, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.CC.MSC.BIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE)

table <- latex(cover.CC.MSC.AIC, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 

table <- latex(cover.CC.MSC.HQ, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = TRUE) 


#Clean up
rm(list = ls())
