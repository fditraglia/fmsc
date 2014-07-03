library(MASS)
n.sims <- 1000
set.seed(7382)

#Calculate asymptotic variance matrix for tau.hat and its inverse
tau.hat <- drop(fmsc.ingredients$tau)
Omega.full <- fmsc.ingredients$Omega[[8]]
Psi <- fmsc.ingredients$Psi
tau.var <- Psi %*% Omega.full %*% t(Psi)
tau.var.inv <- chol2inv(chol(tau.var))

#Function that returns TRUE if a vector x lies in the
#(1 - delta) * 100% confidence region for tau based on 
#the estimate tau.hat and its asymptotic variance matrix 
#tau.var. Defaults to a 95% confidence region.
R.tau <- function(x, delta = 0.05){
  q <- length(tau.hat)
  stopifnot(length(x) == q)  
  crit <- qchisq(1 - delta, q)
  diff <- tau.hat - x
  t(diff) %*% tau.var.inv %*% diff < crit
}

M <- t(mvrnorm(n = n.sims,
               mu = rep(0, ncol(Omega.full)),
               Sigma = Omega.full))

tau.star <- tau.hat

K <- fmsc.ingredients$K 
K.suspect <- lapply(K, function(mat) mat[,-c(1:3)])

cand <- fmsc.ingredients$cand
cand <- lapply(seq_len(ncol(cand)), 
               function(i) as.logical(cand[,i]))

Psi.Omega.Psi <- lapply(cand, function(x) tau.var[x,x])
tt <- lapply(cand, function(x) tau.hat[x] %o% tau.hat[x])
mapply(function(K,x,y) K%*% (x - y) %*% t(K), K.suspect, tt, Psi.Omega.Psi)


B2 <- lapply(seq_along(K.suspect), function(i) 
      K.suspect[[i]] %*% Psi.Omega.Psi[[i]] %*% t(K.suspect[[i]])) 

K.suspect[[8]] %*%(tt[[8]] - Psi.Omega.Psi[[8]]) %*% t(K.suspect[[8]])
fmsc.ingredients$sqbias[,,8]


Psi.M  <- Psi %*% M

n <- nrow(CGdata)
