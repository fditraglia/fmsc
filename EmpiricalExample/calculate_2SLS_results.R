#Frank DiTraglia
#Last Updated: June 14th, 2014

#This script calculates 2SLS Results for all instrument sets considered in my empirical example. It should not be run on its own: it is called by run_empirical_example.R


#-----------------------------------------------
#         DEFINE SECOND STAGE MODEL
#-----------------------------------------------
model <- formula(lngdpc ~ rule + malfal)

#-----------------------------------------------
#        DEFINE INSTRUMENT BLOCKS
#-----------------------------------------------
Baseline <- c("lnmort", "maleco")
Climate <- c("frost", "humid", "latitude")
Europe <- c("eurfrac", "engfrac")
Openness <- c("coast", "trade")
MalfalSq <- c("malfal.sq")
RuleSq <- c("rule.sq")

instrument.blocks <- list(
  c("Baseline"),
  c("Baseline", "Climate"),
  c("Baseline", "Openness"),
  c("Baseline", "Europe"),
  c("Baseline", "Climate", "Europe"),
  c("Baseline", "Climate", "Openness"),
  c("Baseline", "Openness", "Europe"),
  c("Baseline","Climate", "Openness", "Europe"),
  c("Baseline", "MalfalSq"),
  c("Baseline", "RuleSq"),
  c("Baseline", "MalfalSq", "RuleSq"),
  c("Baseline", "Climate", "Openness", "Europe", "MalfalSq", "RuleSq"))


#-----------------------------------------------
#          LIST OF INSTRUMENT SETS
#-----------------------------------------------
f <- function(x){
  command <- paste("c(", paste(x, collapse = ", "), ")", sep = "")
  eval(parse(text = command))
}

instrument.sets <- lapply(instrument.blocks, f)
rm(f, Baseline, Climate, Europe, Openness, MalfalSq, RuleSq)

#-----------------------------------------------
#      FIT 2SLS FOR EACH INSTRUMENT SET 
#-----------------------------------------------
g <- function(instrument.set){
  first.stage <- paste(instrument.set, collapse = " + ")
  first.stage <- paste("~", first.stage)
  first.stage <- as.formula(first.stage)
  tsls(model, first.stage, CGdata)
}

tsls.fits <- lapply(instrument.sets, g)

rm(g, model)

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

tsls.summaries <- lapply(tsls.fits, reg.table)
rm(reg.table)

#There is a very slight discrepancy (a difference of up to 0.02 in a few places) between these results and those in the original version of my paper. This comes from the way I've chosen to handle missing observations: I now exclude Vietnam from the dataset completely, since we only observe the baseline instruments for this country and I want to hold everything constant except the instruments in this exercise. In constrast CG use Vietnam when it is available, as did I in the original version of the empirical example. Again, the difference in results is miniscule.