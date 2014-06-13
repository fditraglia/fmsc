#Frank DiTraglia
#Last Updated: June 13th, 2014
#This script runs the empirical example.


#My empirical example is based on Carstensen & Gundlach (2006) which I abbreviate C&G in the comments below. Specifically, I examine the instrument selection exercise from Table 2 of their paper. This table uses lngdpc as the dependent variable and rule and malfal as the independent variables throughout. The panels of the table look at how adding additional instruments changes the results.

#The instrument blocks considered are:
#     Baseline: lnmort, maleco
#     Climate:  frost, humid, latitude
#     Europe:   eurfrac, engfrac
#     Openness: coast, trade

#Note that there are some missing values in the dataset, so Table 2 is based on a subset including only 44 observations. When estimating with the baseline instruments we could use 45 observations, since lnmort and maleco are observed for Vietnam.


setwd("~/fmsc/EmpiricalExample/")

#Load packages
library(sem) #contains tsls routine

source("load_and_clean_data.R")

source("table_2SLS_results.R")