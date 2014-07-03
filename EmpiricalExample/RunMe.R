#Frank DiTraglia
#Last Updated: July 2nd, 2014
#This script runs the empirical example.

#My empirical example is based on Carstensen & Gundlach (2006) which I abbreviate C&G in the comments below. Specifically, I examine the instrument selection exercise from Table 2 of their paper. This table uses lngdpc as the dependent variable and rule and malfal as the independent variables throughout. The panels of the table look at how adding additional instruments changes the results.

#Note that when estimating with the baseline instruments we could use 45 observations (as CG do in their paper), since lnmort and maleco are observed for Vietnam. We choose not to, however, so that everything is fixed across different specifications except the instruments used. The difference is miniscule in any case.

library(sem) #contains tsls routine
library(Rcpp)
library(RcppArmadillo)

sourceCpp("functions.cpp")

source("load_and_clean_data.R")

source("calculate_2SLS_results.R")

source("make_2SLS_table.R")

source("calculate_fmsc_ingredients.R")

source("calculate_fmsc_values.R")

source("make_fmsc_table.R")

#source("calculate_fmsc_CIs.R")

#source("make_CI_table.R")