# Note: the algorithm given here is designed to be as efficient
# as possible for the *particular* target parameters used in
# the empirical example, but it does *not generalize* to other
# settings. See the comments below for more details. 

#Calculate asymptotic variance matrix for tau.hat and its inverse
tau.hat <- drop(fmsc.ingredients$tau)
Omega.full <- fmsc.ingredients$Omega[[8]]
Psi <- fmsc.ingredients$Psi
tau.var <- Psi %*% Omega.full %*% t(Psi)
tau.var.inv <- chol2inv(chol(tau.var))

#Names of the regressors
reg.names <- rownames(fmsc.ingredients$est)

#Draw the simulations in advance
n.sims <- 10000
set.seed(7382)
M <- t(mvrnorm(n = n.sims,
               mu = rep(0, ncol(Omega.full)),
               Sigma = Omega.full))

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
# Psi * M where each *column* of the result is a sim draw
Psi.M <- Psi %*% M

# Now rather than computing the whole outer product as a 
# function of tau.star and M, we can get the diagonal elements
# for *all draws of M simultaneously* as follows
B1.S <- function(tau.star, row){
  sapply(seq_along(cand), function(i)
    (K.suspect[[i]][row,] %*% (Psi.M + tau.star)[cand[[i]],])^2)
}

# Now we're ready to start putting the pieces back together.
sqbias <- function(tau.star, row){
  B1 <- B1.S(tau.star, row)
  B2 <- B2.S[row,]
  apply(B1, 1, function(x) x - B2)
}


sqbias.pos <- function(tau.star, row){
  pmax(sqbias(tau.star, row), 0)
}

FMSC <- function(tau.star, row){
  criterion <- apply(sqbias(tau.star, row), 2, 
                     function(x) x + avar.S[row,])
  apply(criterion, 2, which.min)
}

posFMSC <- function(tau.star, row){
  criterion <- apply(sqbias.pos(tau.star, row), 2, 
                     function(x) x + avar.S[row,])
  apply(criterion, 2, which.min)
}


Lambda <- function(tau.star, row, criterion){
  index <- criterion(tau.star, row)
  g <- function(i){
    K.row <- K[[index[i]]][row,,drop = FALSE]
    M.ind <- c(TRUE, TRUE, TRUE, cand[[index[i]]])
    -1 * K.row %*% M[M.ind, i, drop = FALSE]
  }
  sapply(seq_along(index), g)
}


