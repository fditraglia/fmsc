#Frank DiTraglia
#September 22nd, 2011

#-------------------------------------------------------#
# Describe the data! Summary stats! What are the units? #
#-------------------------------------------------------#


#This script replicates the main results from Carstensen and Gundlach (2006) using their dataset.


#Work computer
root.dir <- "C:/Documents and Settings/rt356/My Documents/My Dropbox/PhD Dissertation/Job Market Paper/Empirical Examples/Carstensen and Gundlach (Development)"


#Home computer
#root.dir <- "/Users/Frank/Dropbox/PhD Dissertation/Job Market Paper/Empirical Examples/Carstensen and Gundlach (Development)"


#Load packages
library(sem) #contains tsls routine
library(gmm) #to get the J-stat for Andrews' (1999) criteria

#Read the data. The value -999.999 is used to indicate missing data.
setwd(root.dir)
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
#                  Replicate Table 1                    #
#-------------------------------------------------------#


#We have values of lnmort for the following 46 countries:
countries <- clean.data$country[!is.na(clean.data$lnmort)]
cbind(1:length(countries), countries)
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
#[46,] "Vietnam"

#Now, of these countries there are a few missing values for the data used in the baseline models. 
baseline.data <- subset(clean.data, !is.na(lnmort))
baseline.data$country[which(is.na(baseline.data$lngdpc))] #None missing
baseline.data$country[which(is.na(baseline.data$rule))] #Mauritania missing
baseline.data$country[which(is.na(baseline.data$malfal))] #None missing
baseline.data$country[which(is.na(baseline.data$gadp))] #Vietnam missing
baseline.data$country[which(is.na(baseline.data$exprop))] #Mauritania missing
baseline.data$country[which(is.na(baseline.data$lnmort))] #None missing
baseline.data$country[which(is.na(baseline.data$maleco))] #None missing
#We see that Mauritania is missing a value for rule and exprop, while Vietnam is missing a value for gadp. 

#In the paper, Carstensen and Gundlach report use 45 observations for each of the baseline models. For baseline models 1, 2, and 4 the country excluded from the above list is Mauritania. For baseline model 3 the excluded country is Vietname. This is because Mauritania is missing values for rul and exprop, while Vietnam is missing an observation for gadp.


#Formulas for baseline models from Table 1 -- Different measures of institutions and malaria risk as regressors, baseline instrument set
baseline1 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco, clean.data)
baseline2 <- tsls(lngdpc ~ rule + malrisk, ~ lnmort + maleco, clean.data)
baseline3 <- tsls(lngdpc ~ gadp + malfal, ~ lnmort + maleco, clean.data)
baseline4 <- tsls(lngdpc ~ exprop + malfal, ~ lnmort + maleco, clean.data)

#Function to extract quantities needed to replicate Table 1 of the paper from the tsls objects given above
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

#Summarize the results of the Baseline Models to replicate Table 1
baseline.results1 <- reg.table(baseline1)
baseline.results2 <- reg.table(baseline2)
baseline.results3 <- reg.table(baseline3)
baseline.results4 <- reg.table(baseline4)

cbind(baseline.results1, baseline.results2, baseline.results3, baseline.results4) #This matches Table 1 exactly!


#-------------------------------------------------------#
#                  Replicate Table 2                    #
#-------------------------------------------------------#

#This table uses rule and malfal as regressors and looks at how adding additional instruments changes the results. 

#The baseline set of instruments to which others are added is lnmort and maleco.

#The additional blocks of instruments considered here are:
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

#So for this Table we're working with the same countries as in the First Baseline Model (namely the 46 for which we have lnmort excluding Mauritania for which we don't have an observation of rule), and Vietnam. 

#Baseline instruments plus climate block: frost, humid, latitude
instruments1 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + frost + humid + latitude, clean.data)

#Baseline instruments plus Europe block: eurfrac, engfrac
instruments2 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + eurfrac + engfrac, clean.data)

#Baseline instruments plus openness block: coast, trade
instruments3 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + coast + trade, clean.data)

#All instruments
instruments4 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + frost + humid + latitude + eurfrac + engfrac + coast + trade, clean.data)


#Summarize the results and replicate Table 2
instruments.results1 <- reg.table(instruments1)
instruments.results2 <- reg.table(instruments2)
instruments.results3 <- reg.table(instruments3)
instruments.results4 <- reg.table(instruments4)

cbind(instruments.results1, instruments.results2, instruments.results3, instruments.results4) #Matches Table 2 Perfectly!

#At the bottom of the Table are Andrews' (1999) GMM-BIC and GMM-HQ criteria for each of the instrument sets. To replicate these, we need the J-statistic which requires that we carry out efficient GMM rather than tsls. 

#Under homoscedasticity, tsls is efficient GMM, but I'm not sure if that's the assumption they're using here. Also, for the Andrews Criteria we need to used a centered covariance matrix.

#Might want to use some of the code from my simulation study here...
#Alternatively, could extract the objects generated by gmm to create a centered covariance matrix manually
#Probably can assume iid covariance matrix...

#Create blocks of instruments to feed to the function gmm
attach(clean.data)
baseline <- cbind(lnmort, maleco)
openness <- cbind(trade, coast)
climate <- cbind(frost, humid, latitude)
europe <- cbind(eurfrac, engfrac)
detach(clean.data)


#Full.data <- na.omit(data.frame(lngdpc = clean.data$lngdpc, rule = clean.data$rule, malfal = clean.data$malfal, baseline, europe))
#constant <- rep(1, NROW(Full.data))
#y <- Full.data[,1]
#X <- data.frame(constant, Full.data[,2:3])
#All exogenous regressors go in the instrument set: here only the constant
#Z.valid <- data.frame(constant, Full.data[,4:5])
#Z.additional <- data.frame(Full.data[,6:7])
#X <- as.matrix(X)
#Z.valid <- as.matrix(Z.valid)
#Z.additional <- as.matrix(Z.additional)



#-------------------------------------------------------#
#                  Replicate Table 3                    #
#-------------------------------------------------------#

#This table aims to assess the validity of the instruments lnmort and maleco

#Part of this involves using other instruments to identify the model and seeing the result. Are these instruments considered more likely to be valid? Which are they?



#-------------------------------------------------------#
#                     My Additions                      #
#-------------------------------------------------------#

#What about using squared terms or interactions of the regressors themselves as instruments? These will be valid instruments under the assumption that there are no nonlinear effects of the regressors. 

#Perhaps work with the baseline specification since there's no evidence of weak instruments in this case according to the Cragg-Donald statistic from the paper.

#Need to decide on a set of instruments to assume valid. Should it be the baseline instruments? The others? Need to make sure that my baseline instruments aren't weak since my method doesn't address weak instruments.


#Is rule^2 a strong instrument for rule?
attach(clean.data)

x <- rule[!is.na(rule) & !is.na(malfal)]
cor(x, x^2) #0.556 very strong!
cor(x, x^3) #0.891 extremely strong!

#How about malfal?
z <- malfal[!is.na(rule) & !is.na(malfal)]
cor(z, z^2) #0.978 extremely strong! Probably too strong in fact. Ah, this is because so many of the values are zero or one, and zero and one both equal their squares! In contrast, x is a continuous variable

#Interactions?
cor(x, z*x) #0.61 
cor(z, z*x) #-0.72


my.reg1 <- tsls(lngdpc ~ rule + malfal, ~ lnmort + maleco + rule^2 + rule * malfal, clean.data)

my.reg2 <- tsls(lngdpc ~ rule + malfal, ~ rule^2 + rule * malfal, clean.data)

reg.table(my.reg1)
reg.table(my.reg2)


summary(lm(lngdpc ~ rule + malfal))

detach(clean.data)
