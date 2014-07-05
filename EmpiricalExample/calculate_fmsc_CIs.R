#Test out for rule as target parameter, FMSC selection
reg.names[2]
r <- 2
crit <- FMSC
delta <- 0.05

lower <- function(tau.star){
  L.draws <- Lambda(tau.star, r, crit)
  quantile(L.draws, delta / 2)
}


lower(tau.hat)

