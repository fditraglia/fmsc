#Frank DiTraglia
#Last Updated: October 25th, 2011

#This script applies FMSC moment selection and confidence interval correction to the Carstensen and Gundlach (2006) dataset. 



#-------------------------------------------------------#
#----------------------PRELIMINARIES--------------------#
#-------------------------------------------------------#

#Set Working Directory
#office computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Empirical Examples/Carstensen and Gundlach (Development)"

#Home computer
#root.dir <- "/Users/Frank/Dropbox/PhD Dissertation/Job Market Paper/Empirical Examples/Carstensen and Gundlach (Development)"

#Load package for two-stage least squares
library(sem)

#Load my functions for FMSC moment selection and confidence interval correction
setwd(root.dir)
source("FMSC_linear_2SLS.R")


#-------------------------------------------------------#
#-----------------LOAD AND CLEAN DATASET----------------#
#-------------------------------------------------------#

#Load data: NAs are coded as -999.99
input.data <- read.csv("Carstensen_Gundlach.csv", header = TRUE, na.strings = "-999.999")

#Change country from factor data to character data
#is.character(data$country) #FALSE
input.data$country <- as.character(input.data$country)

#Some of the variable names in the dataset don't match those in the paper. Change them so they do.
attach(input.data)
rule <- kaufman
malfal <- mfalrisk
exprop <- exprop2
lngdpc <- lngdpc95
trade <- frarom
latitude <- lat
coast <- landsea
clean.data <- data.frame(country, lngdpc, rule, malfal, malrisk, gadp, exprop, lnmort, maleco, frost, humid, latitude, eurfrac, engfrac, coast, trade)
detach(input.data)
clean.data$country <- as.character(clean.data$country)

#-------------------------------------------------------#
#-------------------FUNCTION: reg.table-----------------#
#-------------------------------------------------------#

#Convenience function to summarize the output of the two-stage least squares function tsls so that it matches the format of the tables in the paper.
reg.table <- function(tsls.object, intercept = FALSE, conf = 0.95, rounded = 2){
  
  #Estimated coefficients
  coeffs <- tsls.object$coefficients
  
  #Estimated covariance matrix
  V <- tsls.object$V
  
  #Standard errors (1-a)*100%
  SE <- sqrt(diag(V))
  
  #(1-a) * 100% Confidence Interval
  a <- 1 - conf
  N <- tsls.object$n
  df <- N - length(coeffs)
  z <- qt(1 - (a/2), df)
  lower <- coeffs - z * SE
  upper <- coeffs + z * SE
  
  results <- rbind(coeffs, SE, lower, upper) 
  results <- round(results, rounded)
  
  #If the option has been invoked, exclude the intercept from the results
  if(intercept == FALSE){results <- results[,-1]}
  
  return(results)
  
}#END reg.table



#-------------------------------------------------------#
#-----------REPLICATE TABLE 2 FROM THE PAPER------------#
#-------------------------------------------------------#

#Table 2 from the paper uses the regressors rule and malfal and looks at how adding additional instruments to the baseline instrument set (lnmort and maleco) changes the results.

#The additional instruments are given in the following blocks:
#     Climate:  frost, humid, latitude
#     Europe:   eurfrac, engfrac
#     Openness: coast, trade

#Are any values of these missing? 
attach(clean.data)
country[is.na(engfrac)] #Vietnam 
country[is.na(eurfrac)] #Vietnam
country[is.na(frost)] #Czechoslovakia, East Germany, USSR, Yugoslavia 
country[is.na(humid)] #Czechoslovakia, East Germanny, USSR, Yugoslavia
country[is.na(latitude)] #Vietnam 
country[is.na(coast)] #Czechoslovakia, East Germany, Puerto Rico, USSR
country[is.na(trade)] #East Germany, Vietnam

#At least one of these missing:
something.missing <- is.na(engfrac) | is.na(eurfrac) | is.na(frost) | is.na(humid) | is.na(latitude) | is.na(coast) | is.na(trade)
country[something.missing] #Czechoslovakia, East Germany, Puerto Rico, USSR, Vietnam, Yugoslavia

#Are there any countries missing values for these variables that nonetheless have observations for lnmort, maleco, rule, and maleco?
missing.baseline <- is.na(lnmort) | is.na(maleco) | is.na(rule) | is.na(maleco)
country[!missing.baseline & something.missing] #Vietnam
detach(clean.data)

#So for this Table we're working with the same countries as in the First Baseline Model from Table 1. We have lnmort for 46 countries, namely: 
# [1,] "Angola"    
# [2,] "Argentina" 
# [3,] "Australia" 
# [4,] "Bangladesh"
# [5,] "Bolivia"   
# [6,] "Brazil"    
# [7,] "Burkina Fa"
# [8,] "Cameroon"  
# [9,] "Canada"    
#[10,] "Chile"     
#[11,] "Colombia"  
#[12,] "Costa Rica"
#[13,] "Dom. Rep." 
#[14,] "Ecuador"   
#[15,] "El Salvado"
#[16,] "Guatemala" 
#[17,] "Guinea"    
#[18,] "Honduras"  
#[19,] "Hong Kong" 
#[20,] "India"     
#[21,] "Indonesia" 
#[22,] "Cote d'Ivo"
#[23,] "Jamaica"   
#[24,] "Kenya"     
#[25,] "Madagascar"
#[26,] "Malaysia"  
#[27,] "Mauritania"
#[28,] "Mexico"    
#[29,] "Morocco"   
#[30,] "New Zealan"
#[31,] "Nigeria"   
#[32,] "Pakistan"  
#[33,] "Panama"    
#[34,] "Paraguay"  
#[35,] "Peru"      
#[36,] "Senegal"   
#[37,] "Singapore" 
#[38,] "South Afri"
#[39,] "Sri Lanka" 
#[40,] "Tanzania"  
#[41,] "Trinidad a"
#[42,] "Tunisia"   
#[43,] "Uruguay"   
#[44,] "USA"       
#[45,] "Venezuela" 
#[46,] "Vietnam
#From these we lack an observation of rule for Mauritania, various of the additional instruments for Vietnam. 

#Baseline instruments only
baseline <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco, clean.data)

#Baseline instruments plus climate block: frost, humid, latitude
instruments1 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + frost + humid + latitude, clean.data)

#Baseline instruments plus Europe block: eurfrac, engfrac
instruments2 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + eurfrac + engfrac, clean.data)

#Baseline instruments plus openness block: coast, trade
instruments3 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + coast + trade, clean.data)

#All instruments
instruments4 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + frost + humid + latitude + eurfrac + engfrac + coast + trade, clean.data)


#Summarize the results and replicate Table 2 (with baseline model from Table 1 added in)
baseline.results <- reg.table(baseline)
instruments.results1 <- reg.table(instruments1)
instruments.results2 <- reg.table(instruments2)
instruments.results3 <- reg.table(instruments3)
instruments.results4 <- reg.table(instruments4)

cbind(baseline.results, instruments.results1, instruments.results2, instruments.results3, instruments.results4) #Results match Table 2 exactly




#-------------------------------------------------------#
#-----------------FMSC MOMENT SELECTION-----------------#
#-------------------------------------------------------#

#Instrument blocks considered in Table 2 of the paper
attach(clean.data)
baseline <- cbind(lnmort, maleco)
openness <- cbind(trade, coast)
climate <- cbind(frost, humid, latitude)
europe <- cbind(eurfrac, engfrac)
detach(clean.data)


#My functions use matrix notation rather than R's formula objects to specify the 2SLS model, so I need to create a matrix of all complete cases of the instruments and regressors I'll use.
Full.data <- na.omit(data.frame(lngdpc = clean.data$lngdpc, rule = clean.data$rule, malfal = clean.data$malfal, baseline, climate, europe, openness))

#The outcome variable: log real gdp per capita at PPP in 1995 in international dollars
y <- Full.data[,1]

#Include a constant term by appending a column of ones
constant <- rep(1, NROW(Full.data))

#The regressors
X <- data.frame(constant, Full.data[,2:3])
#The columns of X are: constant, rule, malfal
X <- as.matrix(X) #My functions require matrix input

#All exogenous regressors go in the instrument set: here only the constant is exogenous (we instrument both rule and malfal)

#The baseline instruments are assumed to be valid
Z.valid <- data.frame(constant, Full.data[,4:5]) 
#The columns of Z.valid are: constant, lnmort, maleco
Z.valid <- as.matrix(Z.valid) #My functions require matrix input

#The additional instruments may be invalid
Z.additional <- data.frame(Full.data[,6:12])
#The columns of Z.additional are: 
#frost, humid, latitude, eurfrac, engfrac, trade, coast
Z.additional <- as.matrix(Z.additional) #My functions require matrix input

#Encode the instrument blocks used in Table 2: each row of this matrix corresponds to the columns of Z.additional (the potentially invalid instruments that we're considering including). A zero in a given column means that the instrument in that column of Z.additional is not included; a one means that it is included.

#Recall that the columns of Z.additional are:
#frost, humid, latitude, eurfrac, engfrac, trade, coast
#Thus we have:
#climate block: columns 1-3 of Z.additional
#Europe block: columns 4-5 of Z.additional
#openness block: columns 6-7 of Z.additional

climate.instruments <- c(1, 1, 1, 0, 0, 0, 0) #include climate block
europe.instruments <- c(0, 0, 0, 1, 1, 0, 0) #include Europe block
openness.instruments <- c(0, 0, 0, 0, 0, 1, 1) #include openness block

additional.instruments <- as.matrix(rbind(climate.instruments,
                                                europe.instruments,
                                                openness.instruments)) 
#We don't need to specify the baseline model or the model including all instruments, as these are included automatically by my function.

#In the paper, the stated target is the disease-income link, that is the coefficient on malfal. Since the columns of X are constant, rule, malfal this means we're interested in the third coefficient
mu.malfal <- function(x){return(x[3])}
nabla.malfal <- function(x){return(c(0, 0, 1))}
FMSC.malfal <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional, additional.instruments, mu.malfal, nabla.malfal)


#What if the target were rule?
mu.rule <- function(x){return(x[2])}
nabla.rule <- function(x){return(c(0, 1, 0))}
FMSC.rule <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional, additional.instruments, mu.rule, nabla.rule)

#In each of these cases, the FMSC tells us to use the full model.



#Variation #1: what if we consider including some of the instrument blocks together?
climate.only <- c(1, 1, 1, 0, 0, 0, 0) 
europe.only <- c(0, 0, 0, 1, 1, 0, 0) 
openness.only <- c(0, 0, 0, 0, 0, 1, 1)
climate.europe <- c(1, 1, 1, 1, 1, 0, 0)
climate.openness <- c(1, 1, 1, 0, 0, 1, 1)
europe.openness <- c(0, 0, 0, 1, 1, 1, 1)

additional.instruments.var1 <- as.matrix(rbind(climate.only,
                                               europe.only,
                                               openness.only,
                                               climate.europe,
                                               climate.openness,
                                               europe.openness))

FMSC.malfal.var1 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional, additional.instruments.var1, mu.malfal, nabla.malfal)

FMSC.rule.var1 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional, additional.instruments.var1, mu.rule, nabla.rule)
#Again, the FMSC chooses the FMSC in both cases, but when the target is rule, the full instrument is practically indistinguishable from climate plus europe in terms of FMSC.


#Variation #2: What about at about including squares of the endogenous regressors? 
#Recall that the columns of X are: constant, rule, malfal
rule.sq <- X[,2]^2
malfal.sq <- X[,3]^2

Z.additional.var2 <- cbind(rule.sq, malfal.sq)
#The columns of Z.additional.var2 are: rule.sq, malfal.sq
Z.additional.var2 <- as.matrix(Z.additional.var2)

rule.squared <- c(1, 0)
malfal.squared <- c(0, 1)

additional.instruments.var2 <- as.matrix(rbind(rule.squared,
                                               malfal.squared))

FMSC.malfal.var2 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional.var2, additional.instruments.var2, mu.malfal, nabla.malfal)

FMSC.rule.var2 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional.var2, additional.instruments.var2, mu.rule, nabla.rule)
#Here, the FMSC says to use the Full model if the target is rule, but to add malfal.sq if the target is malfal






#Variation #3: What if we also consider the instrument blocks from Table 2?

Z.additional.var3 <- cbind(Z.additional, Z.additional.var2)
#The columns of Z.additional are: 
#frost, humid, latitude, eurfrac, engfrac, trade, coast
#The columns of Z.additional.var2 are: rule.sq, malfal.sq
#Hence, the columns of Z.additional.var3 are:
#frost, humid, latitude, eurfrac, engfrac, trade, coast, rule.sq, malfal.sq

climate.block <-    c(1, 1, 1, 0, 0, 0, 0, 0, 0) 
europe.block <-     c(0, 0, 0, 1, 1, 0, 0, 0, 0) 
openness.block <-   c(0, 0, 0, 0, 0, 1, 1, 0, 0)
rule.sq <-          c(0, 0, 0, 0, 0, 0, 0, 1, 0)
malfal.sq <-        c(0, 0, 0, 0, 0, 0, 0, 0, 1)
endog.sq <-         c(0, 0, 0, 0, 0, 0, 0, 1, 1)

additional.instruments.var3 <- as.matrix(rbind(climate.block,
                                               europe.block,
                                               openness.block,
                                               rule.sq,
                                               malfal.sq,
                                               endog.sq))


FMSC.malfal.var3 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional.var3, additional.instruments.var3, mu.malfal, nabla.malfal)

FMSC.rule.var3 <- FMSC.2SLS.conf.naive(X, y, Z.valid, Z.additional.var3, additional.instruments.var3, mu.rule, nabla.rule)
#Interesting: once again, chooses the full model in each case! I wonder how wide the confidence intervals would be...

#Need to make a table of the 2SLS results in each case. Presumably, they're all giving more or less the same point estimates. Since the result is not sensitive to the instrument choice, the FMSC is telling us to include them all.

#Append squared endogenous regressors to the cleaned dataset 
augmented.data <- cbind(clean.data, rule.sq = clean.data$rule^2, malfal.sq = clean.data$malfal^2)

#Fit all of the 2SLS models over which we apply the FMSC above.
baseline.IV <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco, augmented.data)

add.openness <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + trade + coast, augmented.data)

add.climate <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + frost + humid + latitude, augmented.data)

add.europe <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + eurfrac + engfrac, augmented.data)

climate.europe.IV <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + frost + humid + latitude + eurfrac + engfrac, augmented.data)
  
climate.openness.IV <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + frost + humid + latitude + trade + coast, augmented.data)

europe.openness.IV <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + eurfrac + engfrac + trade + coast, augmented.data)

europe.climate.openness <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + eurfrac + engfrac + frost + humid + latitude + trade + coast, augmented.data)

add.rule.sq <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + rule.sq, augmented.data)

add.malfal.sq <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + malfal.sq, augmented.data)

add.endog.sq <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + rule.sq + malfal.sq, augmented.data)

all.instruments <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + eurfrac + engfrac + frost + humid + latitude + trade + coast + rule.sq + malfal.sq, augmented.data)


#How similar are the results?
reg.table(baseline.IV)
reg.table(add.climate)
reg.table(add.openness)
reg.table(add.europe)
reg.table(climate.europe.IV)
reg.table(climate.openness.IV)
reg.table(europe.openness.IV)
reg.table(europe.climate.openness)
reg.table(add.malfal.sq)
reg.table(add.rule.sq)
reg.table(add.endog.sq)
reg.table(all.instruments)
#The results are all incredibly similar. The differences are small relative to the standard errors. 

#Still need to do the Fully Valid confidence interval correction, but I'm having trouble getting the optimization routine to work.

#Make summary tables of this information for the paper
instrument.names <- c('Baseline', 'Climate', 'Openness', 'Europe', 'Malfalsq', 'Rulesq')


baseline.include <- c(1, 0, 0, 0, 0, 0)
add.climate.include <- c(1, 1, 0, 0, 0, 0)
add.openness.include <- c(1, 0, 1, 0, 0, 0)
add.europe.include <- c(1, 0, 0, 1, 0, 0)
climate.europe.include <- c(1, 1, 0, 1, 0, 0)
climate.openness.include <- c(1, 1, 1, 0, 0, 0) 
europe.openness.include <- c(1, 0, 1, 1, 0, 0)
europe.climate.openness.include <- c(1, 1, 1, 1, 0, 0)
malfal.sq.include <- c(1, 0, 0, 0, 1, 0)
rule.sq.include <- c(1, 0, 0, 0, 0, 1)
endog.sq.include <- c(1, 0, 0, 0, 1, 1)
all.include <- c(1, 1, 1, 1, 1, 1)

instrument.names[as.logical(baseline.include)]
instrument.names[as.logical(add.climate.include)]
instrument.names[as.logical(add.openness.include)]
instrument.names[as.logical(add.europe.include)]
instrument.names[as.logical(climate.europe.include)]
instrument.names[as.logical(climate.openness.include)]
instrument.names[as.logical(europe.openness.include)]
instrument.names[as.logical(europe.climate.openness.include)]
instrument.names[as.logical(malfal.sq.include)]
instrument.names[as.logical(rule.sq.include)]
instrument.names[as.logical(endog.sq.include)]
instrument.names[as.logical(all.include)]


results1 <- reg.table(baseline.IV)
instruments1 <- cbind(baseline.include, baseline.include)
row.names(instruments1) <- instrument.names
results1 <- rbind(results1, instruments1)


results2 <- reg.table(add.climate)
instruments2 <- cbind(add.climate.include, add.climate.include)
row.names(instruments2) <- instrument.names
results2 <- rbind(results2, instruments2)

results3 <- reg.table(add.openness)
instruments3 <- cbind(add.openness.include, add.openness.include)
row.names(instruments3) <- instrument.names
results3 <- rbind(results3, instruments3)

results4 <- reg.table(add.europe)
instruments4 <- cbind(add.europe.include, add.europe.include)
row.names(instruments4) <- instrument.names
results4 <- rbind(results4, instruments4)

results5 <- reg.table(climate.europe.IV)
instruments5 <- cbind(climate.europe.include, climate.europe.include)
row.names(instruments5) <- instrument.names
results5 <- rbind(results5, instruments5)

results6 <- reg.table(climate.openness.IV)
instruments6 <- cbind(climate.openness.include, climate.openness.include)
row.names(instruments6) <- instrument.names
results6 <- rbind(results6, instruments6)

results7 <- reg.table(europe.openness.IV)
instruments7 <- cbind(europe.openness.include, europe.openness.include)
row.names(instruments7) <- instrument.names
results7 <- rbind(results7, instruments7)

results8 <- reg.table(europe.climate.openness)
instruments8 <- cbind(europe.climate.openness.include, europe.climate.openness.include)
row.names(instruments8) <- instrument.names
results8 <- rbind(results8, instruments8)

results9 <- reg.table(add.malfal.sq)
instruments9 <- cbind(malfal.sq.include, malfal.sq.include)
row.names(instruments9) <- instrument.names
results9 <- rbind(results9, instruments9)


results10 <- reg.table(add.rule.sq)
instruments10 <- cbind(rule.sq.include, rule.sq.include)
row.names(instruments10) <- instrument.names
results10 <- rbind(results10, instruments10)

results11 <- reg.table(add.endog.sq)
instruments11 <- cbind(endog.sq.include, endog.sq.include)
row.names(instruments11) <- instrument.names
results11 <- rbind(results11, instruments11)

results12 <- reg.table(all.instruments)
instruments12 <- cbind(all.include, all.include)
row.names(instruments12) <- instrument.names
results12 <- rbind(results12, instruments12)



Part1 <- cbind(results1, results2, results3, results4, results5, results6)
Part2 <- cbind(results7, results8, results9, results10, results11, results12)




paste.below <- function(FMSC.object){
  
  FMSC <- round(FMSC.object$FMSC,2)
  CI <- round(FMSC.object$CI,2)
  
  CI <- cbind(rep(NA, 3), CI)
  CI <- data.frame(FMSC = CI[,1], Estimate = CI[,2])

  out <- rbind(FMSC, CI)
  
  return(out)
  
}#END paste.below


malfal1 <- paste.below(FMSC.malfal)
malfal2 <- paste.below(FMSC.malfal.var1)
malfal3 <- paste.below(FMSC.malfal.var2)
malfal4 <- paste.below(FMSC.malfal.var3)

rule1 <- paste.below(FMSC.rule)
rule2 <- paste.below(FMSC.rule.var1)
rule3 <- paste.below(FMSC.rule.var2)
rule4 <- paste.below(FMSC.rule.var3)

#Dump everything into LaTeX format
library(Hmisc)

#Assign to a dummy variable name to avoid having R try to compile the LaTeX
out <- latex(Part1, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = FALSE)

#Append remaining tables
out <- latex(Part2, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

out <- latex(malfal1, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

out <- latex(malfal2, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

out <- latex(malfal3, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

out <- latex(malfal4, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

out <- latex(rule1, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

out <- latex(rule2, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

out <- latex(rule3, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

out <- latex(rule4, file = 'Output_Tables.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)



