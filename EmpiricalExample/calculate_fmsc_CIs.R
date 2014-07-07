# Coverage probabilities
alpha <- 0.05
delta <- 0.05

#------------------------------------------------------------
# Constraint Function - only consider values of tau.star
# in (1 - delta) * 100% Confidence Region
#------------------------------------------------------------
# Note that NOMAD requires constraints of the form c(x) <= 0
# so we write this such that a *negative* value indicates
# that tau.star lies within the confidence region and a 
# positive value indicates that it lies outside this region
q <- length(tau.hat)
chisq.crit <- qchisq(1 - delta, q)

constraint <- function(x){
  stopifnot(length(x) == q)
  diff <- tau.hat - x
  drop(t(diff) %*% tau.var.inv %*% diff - chisq.crit)
}

# Element-wise CIs for tau serve as the lower and upper
# bounds for tau.star in the optimization routine
norm.crit <- qnorm(1 - delta / 2)
tau.upper <- tau.hat + sqrt(diag(tau.var)) * norm.crit
tau.lower <- tau.hat - sqrt(diag(tau.var)) * norm.crit

#Sanity Checks
stopifnot(tau.hat < tau.upper)
stopifnot(tau.hat > tau.lower)
stopifnot(constraint(tau.hat) == -chisq.crit)
stopifnot(constraint(tau.lower) == constraint(tau.upper))

# Use closure to construct ``function factories'' for
# different criteria and target parameters
LambdaQuantile <- function(target, criterion, q){
  function(tau.star){
    row <- which(reg.names == target)
    L.draws <- Lambda(tau.star, row, criterion)
    quantile(L.draws, q)
  }
}


#Set up parameters for optimization routine
snomadr.default <- function(f){
  snomadr(eval.f = f, n = q, x0 = tau.hat, bbin = rep(0,q),
          bbout = c(0,1), lb = tau.lower, ub = tau.upper)
}


#FMSC and posFMSC selected estimators of each target parameter
fmsc.values <- lapply(fmsc.values, as.data.frame)
selected <- lapply(fmsc.values, function(x) apply(x[,1:2], 2, which.min))
mu.est <- lapply(fmsc.values, 
                 function(x) x$est[apply(x[,1:2], 2, which.min)])
names(mu.est$rule) <- c("FMSC", "posFMSC")
names(mu.est$malfal) <- c("FMSC", "posFMSC")

# Number of Observations for CI construction
n <- nrow(CGdata)


#Function to calculate various CIs
target.name <- "malfal" 
criterion.name <- "FMSC"
  
  criterion <- match.fun(criterion.name)
  mu.hat <- mu.est[[target.name]][[criterion.name]]
  S.hat <- selected[[target.name]][[criterion.name]]
  
  g <- LambdaQuantile(target.name, criterion, 1 - alpha / 2)
  h <- LambdaQuantile(target.name, criterion, alpha / 2)
  
  Upper <- snomadr.default(function(x) c(-g(x), constraint(x)))
  Lower <- snomadr.default(function(x) c(h(x), constraint(x)))
  
  constraint(Upper$solution)
  constraint(Lower$solution)
  g(Upper$solution)
  h(Lower$solution)
  g(tau.hat)
  h(tau.hat)
  
  one.step <- mu.hat - c(g(tau.hat), 
                         h(tau.hat)) / sqrt(n)
  two.step <- mu.hat - c(g(Upper$solution), 
                         h(Lower$solution)) / sqrt(n)

  SE.naive <- sqrt(diag(tsls.fits[[S.hat]]$V)[[target.name]])
  z.naive <- qnorm(1 - alpha /2)
  naive <- mu.hat + z.naive * SE.naive * c(-1, 1)
  
  out <- rbind(naive, one.step, two.step)
  row.names(out) <- c("Naive", "1-Step", "2-Step")
  colnames(out) <- c("Lower", "Upper")
  


