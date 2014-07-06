# Constraint Function - only consider values of tau.star
# in (1 - delta) * 100% Confidence Region
#
# Note that NOMAD requires constraints of the form c(x) <= 0
# so we write this such that a *negative* value indicates
# that tau.star lies within the confidence region and a 
# positive value indicates that it lies outside this region
delta <- 0.05
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


g <- LambdaQuantile(target = "rule", criterion = FMSC, q = 0.95)


result <- snomadr(eval.f = function(x) c(-g(x), constraint(x)), 
                  n = q,
                  x0 = tau.hat,
                  bbin = rep(0,q),
                  bbout = c(0,1),
                  lb = tau.lower,
                  ub = tau.upper)

result$solution
constraint(result$solution)
g(result$solution)
g(tau.hat)






