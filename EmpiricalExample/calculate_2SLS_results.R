#Frank DiTraglia
#Last Updated: June 13th, 2014

#This script calculates 2SLS Results for all instrument sets considered in my empirical example. It should not be run on its own: it is called by run_empirical_example.R


#-----------------------------------------------
#         DEFINE SECOND STAGE MODEL
#-----------------------------------------------
model <- formula(lngdpc ~ rule + malfal)


#-----------------------------------------------
#    DEFINE INSTRUMENT BLOCKS TO MATCH C&G
#-----------------------------------------------
baseline <- "~lnmort + maleco"
climate <- "frost + humid + latitude"
europe <- "eurfrac + engfrac"
openness <- "coast + trade"

#-----------------------------------------------
#  SET UP INSTRUMENT SETS AS FORMULA OBJECTS 
#-----------------------------------------------

#Numbering indicates how they will be arranged in the table
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

set9 <- as.formula(paste(baseline, "+",
                         "malfal.sq"))

set10 <-as.formula(paste(baseline, "+",
                         "rule.sq"))

set11 <- as.formula(paste(baseline, "+",
                          "malfal.sq", "+",
                          "rule.sq"))

set12 <- as.formula(paste(baseline, "+", 
                          climate, "+", 
                          openness, "+", 
                          europe, "+",
                          "malfal.sq" , "+",
                          "rule.sq"))

rm(baseline, climate, europe, openness)
 
#-----------------------------------------------
#  STORE INSTRUMENT SETS AS A LIST OF FORMULAS 
#-----------------------------------------------
sets.text <- paste('set', 1:12, sep = '', collapse = ', ')
  
instrument.sets <- eval(parse(text = 
                        paste('list(', sets.text, ')', 
                        sep = '')))

eval(parse( text = paste("rm(", sets.text, ")", sep = "")))
rm(sets.text)

#-----------------------------------------------
#      FIT 2SLS FOR EACH INSTRUMENT SET 
#-----------------------------------------------
tsls.fits <- lapply(instrument.sets,
                    function(first.stage){
                      tsls(model, first.stage, CGdata)
                    })
rm(model)

#-----------------------------------------------
# SUMMARIZE 2SLS RESULTS FOR EACH INSTRUMENT SET 
#-----------------------------------------------
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

tsls.results <- lapply(tsls.fits, reg.table)
rm(reg.table)

#There is a very slight discrepancy (a difference of 0.01 in a few places) between these results and those in the original version of my paper. This comes from the way I've chosen to handle missing observations: I now exclude Vietnam from the dataset completely, since we only observe the baseline instruments for this country and I want to hold everything constant except the instruments in this exercise. In constrast CG use Vietnam when it is available, as did I in the original version of the empirical example. Again, the difference in results is miniscule.