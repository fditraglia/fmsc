#September 24th, 2011

#This script creates tables summarizing the simulation results

#---------------------------------------------#
#                   MACROS                    #
#---------------------------------------------#

sim.results.dir <- './N = 500'
infilename <- 'sim_results_500_not_centered.csv'
outfilename <- 'sim_tables_raw_500_not_centered.tex'

#---------------------------------------------#
#                 END MACROS                  #
#---------------------------------------------#


#Package to output tables to LaTeX
library(Hmisc)

#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Coverage Probs not recentered"

setwd(root.dir)
setwd(sim.results.dir)

#Read simulation results
results <- read.csv(infilename, header = TRUE)

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
table <- latex(cover.valid, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = FALSE)

#Append the remaining tables
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
