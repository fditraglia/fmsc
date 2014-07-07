#December 6th, 2011



#---------------------------------------------#
#                   MACROS                    #
#---------------------------------------------#

infile.dir <- './N = 500'
outfilename <- 'sim_tables_bayeslike_correct_conf_FMSC_500.tex'

#---------------------------------------------#
#                 END MACROS                  #
#---------------------------------------------#


#Package to output tables to LaTeX
library(Hmisc)

#For my work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Simulations/Bayes-like Correct Confidence"

setwd(root.dir)
setwd(infile.dir)

#Read simulation results
results <- read.csv('sim_results_FMSC_bayeslike_correct_N_500.csv', header = TRUE)

  


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



conf.bayeslike <- make.table(round(bayeslike.FMSC.correct.cover, 2))



detach(results)


#Output all of the tables in LaTeX tabular format to a table.

#Assign to a bogus variable to prevent the function from trying to load LaTex itself and compile the tables it creates
table <- latex(conf.bayeslike, file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, rowname = NULL, append = FALSE)


#Clean up
rm(list = ls())
