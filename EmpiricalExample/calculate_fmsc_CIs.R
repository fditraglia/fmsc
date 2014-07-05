library(MASS)
n.sims <- 10
set.seed(7382)

#Calculate asymptotic variance matrix for tau.hat and its inverse
tau.hat <- drop(fmsc.ingredients$tau)
Omega.full <- fmsc.ingredients$Omega[[8]]
Psi <- fmsc.ingredients$Psi
tau.var <- Psi %*% Omega.full %*% t(Psi)
tau.var.inv <- chol2inv(chol(tau.var))

#Draw the simulations in advance;
M <- t(mvrnorm(n = n.sims,
               mu = rep(0, ncol(Omega.full)),
               Sigma = Omega.full))

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

# To construct valid confidence intervals, we first need code
# to evaluate FMSC_S(tau.star, M) for each specification S
# at simulation draw M and specified value tau.star for the 
# bias parameter. 
#
# We calculate FMSC from an estimate of the ABIAS and AVAR of
# of a given estimator of theta. The AVAR estimate for each
# specification is:
#
#     AVAR = K.S * Xi.S * Omega * Xi.S' * K.S
#
# This doesn't depend on tau.star or M and we've already 
# computed it
avar.S <- fmsc.ingredients$avar


# The ABIAS.sq estimate depends on the matrix
#                 [ 0  0 ]
#     BIAS_MAT =  [ 0  B ]
#
# where
#
#     B = B1 - B2
#     B1 = (Psi * M + tau.star) * (Psi * M + tau.star)'
#     B2 =  Psi * Omega * Psi
#
# Note that *only* B1 depends on tau.star or M. Moreover,
# we have already computed B2: it is simply the variance 
# matrix of tau.hat
B2 <- tau.var

# Now, BIAS_MAT enters the ABIAS formula as follows:
#
#   ABIAS.sq = K.S * Xi.S * BIAS_MAT * Xi.S' * K.S
#
# Although written as a matrix multiplication, the 
# combination of the Xi.S matrices with the zero 
# elements of BIAS_MAT is really just a pair of 
# subsetting operations that serve to remove the
# the elements of K.S that correspond to the valid
# instruments and remove the elements of B that 
# correspond to instruments from z2 that are not
# used in specification S. Rather than explicitly
# forming BIAS_MAT, we'll just take the appropiate
# subsets.
#
# First, we'll get the bias contribution that comes
# from B2 since this doesn't depend on tau or M.
#
# In this particular application there are three
# instruments in z1: a constant, maleco and lnmort.
# This means that the *first three* columns of K.S
# are "zeroed out" by BIAS_MAT for each S. The 
# remaining columns correspond to the "suspect
# instruments:
K <- fmsc.ingredients$K 
K.suspect <- lapply(K, function(mat) mat[,-c(1:3)])

# To subset B, well use logical indicators constructed
# as follows:
cand <- fmsc.ingredients$cand
cand <- lapply(seq_len(ncol(cand)), 
               function(i) as.logical(cand[,i]))

# Now we can calculate the contribution to the squared
# asymptotic bias that comes from B2:
#         
#                        [ 0  0  ]
#     B2.S = K.S * Xi.S  [ 0  B2 ] * Xi.S' * K.S
#
# If we define B1.S analogously, then:
#
#     ABIAS.sq(S) = B1.S - B2.S
# 
# That is, B2.S enters *negatively*
S.contrib <- function(i, Inner){
  S <- cand[[i]]
  K.S <- K.suspect[[i]]
  K.S %*% Inner[S,S] %*% t(K.S)
}
B2.S <- sapply(seq_along(cand), S.contrib, Inner = B2,
               simplify = "array")

# We can check that this works by re-computing the squared
# asymptotic bias matrices that we already have from the 
# C++ code used to calculate the FMSC in the empirical 
# example
sqbias.cpp <- fmsc.ingredients$sqbias
foo <- sapply(seq_along(cand), S.contrib, 
              Inner = tau.hat %o% tau.hat, simplify = "array")
sqbias.R <- foo - B2.S
all.equal(sqbias.R, sqbias.cpp)  
rm(foo, sqbias.cpp, sqbias.R)

# Now, as it happens there we can exploit a substantial 
# simplification in this example because the target
# parameters we consider are simply *elements* of the
# underlying parameter vector beta. This means that we 
# only need to store elements on the *main diagonal* of
# B2.S and avar.S
avar.S <- apply(avar.S, 3, diag)
B2.S <- apply(B2.S, 3, diag)

# While doesn't save us particularly many computations for
# avar.S and B2.S it's *incredibly* helpful for dealing with 
# B1.S(M, tau.star)
#
#                        [ 0  0  ]
#     B1.S = K.S * Xi.S  [ 0  B1 ] * Xi.S' * K.S
#
#     B1 = (Psi * M + tau.star) * (Psi * M + tau.star)'
#
# since it means we *never* have to explicitly form and
# store the outer product matrices. First, we pre-compute
# Psi * M
Psi.M <- Psi %*% M


B1.diag <- lapply(seq_along(cand), function(i)
          (K.suspect[[i]] %*% (Psi.M + tau.star)[cand[[i]],])^2)
                    
sqbias.diag.sim <- lapply(seq_along(B1.diag), function(i)
                           B1.diag[[i]] - B2.diag[,i])

FMSC.diag.sim <- sapply(seq_along(sqbias.diag.sim), function(i)
                        sqbias.diag.sim[[i]] + avar.diag[,i],
                        simplify = "array")

posFMSC.diag.sim <- sapply(seq_along(sqbias.diag.sim), function(i)
                pmax(sqbias.diag.sim[[i]], 0) + avar.diag[,i],
                simplify = "array")


apply(FMSC.diag.sim[3,,], 1, which.min)
bar <- apply(posFMSC.diag.sim[3,,], 1, which.min)
sapply(bar, function(i) K[[bar[i]]][3,] %*% M[cand[[bar[i]]],i])



#Sanity Check
testy <- apply(K.tt.K, 3, diag) - B2.diag + avar.diag
all.equal(testy[2,], fmsc.values$rule[,1])
all.equal(testy[3,], fmsc.values$malfal[,1])


n <- nrow(CGdata)


