#Frank DiTraglia
#February 22st 2014

setwd("~/fmsc/EmpiricalExample/")

#Load packages
library(sem) #contains tsls routine

#Read the data. The value -999.999 is used to indicate missingness.
CGdata <- read.csv("Carstensen_Gundlach.csv", header = TRUE, na.strings = "-999.999", stringsAsFactors = FALSE)


#Rename columns that don't match the variable names from the paper
paper <- c("rule", "malfal", "exprop", "lngdpc", "trade", "latitude", "coast")
dataset <- c("kaufman", "mfalrisk", "exprop2", "lngdpc95", "frarom", "lat", "landsea")
key <- data.frame(paper, dataset, stringsAsFactors = FALSE)
rm(paper, dataset)

for(i in 1:nrow(key)){
  col <- which(names(CGdata) == key$dataset[i])
  names(CGdata)[col] <- key$paper[i]
}


#I use the FMSC to examine the instrument selection exercise in Table 2 of the paper. This table uses lngdpc as the dependent variable and rule and malfal as the independent variables throughout. The various entries in the table 

#This table uses rule and malfal as RHS variables and lngdpc as the outcome variable. The panels of the table look at how adding additional instruments changes the results.

#The instrument blocks considered are:
#     Baseline: lnmort, maleco
#     Climate:  frost, humid, latitude
#     Europe:   eurfrac, engfrac
#     Openness: coast, trade

#Note that there are some missing values in the dataset, so Table 2 is based on a subset including only 44 observations. When estimating with the baseline instruments we could use 45 observations, since lnmort and maleco are observed for Vietnam.

CGdata <- subset(CGdata, 
                    !is.na(lngdpc) &
                    !is.na(rule) & !is.na(malfal) &
                    !is.na(lnmort) & !is.na(maleco) &   
                    !is.na(frost) & !is.na(humid) & 
                    !is.na(eurfrac) & !is.na(engfrac) & 
                    !is.na(coast) & !is.na(trade))


#For convenience, drop the columns we won't be using
keep <- c("lngdpc", "rule", "malfal", "maleco", "lnmort", "frost", "humid", "latitude", "eurfrac", "engfrac", "coast", "trade")
keep.cols <- which(names(CGdata) %in% keep)
CGdata <- CGdata[,keep.cols]

rm(keep, keep.cols)


#Function to extract quantities needed to replicate Table 1 of the paper from the tsls objects given above
reg.table <- function(tsls.object){
  coef <- tsls.object$coefficients
  SE <- sqrt(diag(tsls.object$V))
  N <- tsls.object$n
  z <- qt(0.975, N - length(coef))
  lower <- coef - z * SE
  upper <- coef + z * SE
  results <- round(rbind(coef, SE, lower, upper), 2) 
  return(results[,-1]) #Don't report intercept
}


#The second stage is the same for all instrument sets
model <- formula(lngdpc ~ rule + malfal)

#Set up instrument blocks to match terminology in the paper
baseline <- "~lnmort + maleco"
climate <- "frost + humid + latitude"
europe <- "eurfrac + engfrac"
openness <- "coast + trade"

#Instrument sets numbered to correspond to panels in table
set1 <- as.formula(paste(baseline))

set2 <- as.formula(paste(baseline, "+", 
                         climate))
set3 <- as.formula(paste(baseline, "+", 
                         openness))

set4 <- as.formula(paste(baseline, "+", 
                         europe))

set5 <- as.formula(paste(baseline, "+", 
                         climate, "+", 
                         europe))

set6 <- as.formula(paste(baseline, "+", 
                         climate, "+", 
                         openness))

set7 <- as.formula(paste(baseline, "+", 
                         openness, "+", 
                         europe))

set8 <- as.formula(paste(baseline, "+", 
                         climate, "+", 
                         openness, "+", 
                         europe))

#Fit 2sls for each instrument set
fit1 <- tsls(model, set1, CGdata)
fit2 <- tsls(model, set2, CGdata)
fit3 <- tsls(model, set3, CGdata)
fit4 <- tsls(model, set4, CGdata)
fit5 <- tsls(model, set5, CGdata)
fit6 <- tsls(model, set6, CGdata)
fit7 <- tsls(model, set7, CGdata)
fit8 <- tsls(model, set8, CGdata)


#Summarize the results and replicate Table 2
results1 <- reg.table(fit1)
results2 <- reg.table(fit2)
results3 <- reg.table(fit3)
results4 <- reg.table(fit4)
results5 <- reg.table(fit5)
results6 <- reg.table(fit6)
results7 <- reg.table(fit7)
results8 <- reg.table(fit8)


#Note that there is a very slight discrepancy (0.01) in the standard errors for fit1 compared to those in C&G's paper. This comes from the fact that they include Vietnam when estimating with the baseline instruments. I exclude it since we don't observe the values of any of the other instruments for Vietnam and I want to hold everything constant except the instruments.
cbind(results1, results2, results3, results4) 
cbind(results5, results6, results7, results8)
