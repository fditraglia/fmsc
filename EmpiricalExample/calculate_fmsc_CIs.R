library(MASS)
n.sims <- 1000
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


# The first step is contructing a function that can evaluate
# FMSC_S(tau.star, M). The "meat" of the sandwich doesn't
# depend on the target parameter so we'll start here. There 
# are two components, one of which depends on M and tau.star
# and one of which doesn't. 

# The AVAR component for each candidate is *fixed* and we've 
# already computed it: these are our estimates of the product
# K.S * XI.S * Omega * Xi.S' * K.S'
avar <- fmsc.ingredients$avar

# The Bias component *does* depend on M and tau.star:
#     K.S * Xi.S * BIAS.MAT * Xi.S' * K.S
# where BIAS.MAT is a 2x2 block matrix whose only
# non-zero block is the lower diagonal block, which is:
#     B = B1 - B2
# where:
#     B1 = (Psi * M + tau.star)(Psi * M + tau.star)'
#     B2 = Psi * Omega * Psi'
# Notice that only B1 depends on M and tau.star. Moreover,
# B2 is simply the variance matrix of tau.hat which we
# already computed above
B2 <- tau.var

# Now, we do not in fact need to form BIAS.MAT. It is
# a waste of storage and flops to explicitly multiply
# all of the zeros with the zeros in Xi.S. The only 
# part of K.S that plays a role in the ABIAS is the 
# part that multiplies the *non-zero* elements of the
# matrix BIAS.MAT. So what we really want is to subset
# K.S. Remember that K.S has as many rows as there are
# regressors, and as many columns as there are instruments
# in specification S. Since we always include the valid
# instruments, we can partition K.S as follows:
#
#   K.S = [K.S(valid)  K.S(suspect)]
#
# where K.S(valid) corresponds to the valid instruments 
# and K.S(suspect) corresponds to the potentially invalid
# instruments included in specification S. Accordingly, 
# we'll first extract K.S(suspect) for each S:
K <- fmsc.ingredients$K 
K.suspect <- lapply(K, function(mat) mat[,-c(1:3)])

# All that pre-multiplying by Xi.S and post-multiplying
# by Xi.S' does is subset the non-zero block of BIAS.MAT
# to correspond to the instruments used in specification
# S. We can do this without explicitly forming the Xi
# matrices by using the candidate indicators that we 
# used in the FMSC calculation. These correspond to the
# columns of z2 used in estimation. We'll extract them
# and convert them from a matrix ofzeros and ones to a
# list of TRUE and FALSE. 
cand <- fmsc.ingredients$cand
cand <- lapply(seq_len(ncol(cand)), 
               function(i) as.logical(cand[,i]))

# Now we can use cand to do the subsetting that Xi
# accomplishes via matrix multiplication. First, 
# we'll do this for B2, turning it into a list where
# each element refers to a particular candidate
B2 <- lapply(cand, function(x) B2[x,x])

# And now pre-multiply by the *relevent* part of K
# and post-multiple by the *relevant* part of K.S
# the result will be a list of matrices each of which
# has the same dimensions
B2 <- lapply(seq_along(K.suspect), function(i)
              K.suspect[[i]] %*% B2[[i]] %*% t(K.suspect[[i]]))

# As a sanity check, let's make sure we can get the same
# squared asymptotic bias matrix as the C++ code gave us
# First we need to extract it and convert it to a list
sqbias.cpp <- fmsc.ingredients$sqbias
sqbias.cpp <- lapply(1:(dim(sqbias.cpp)[3]),
                     function(i) sqbias.cpp[,,i])
tt <- lapply(cand, function(x) tau.hat[x] %o% tau.hat[x])
K.tt.K <- lapply(seq_along(tt), function(i)
                K.suspect[[i]] %*% tt[[i]] %*% t(K.suspect[[i]]))
sqbias.R <- lapply(seq_along(B2), function(i)
                K.tt.K[[i]] - B2[[i]])
all.equal(sqbias.R, sqbias.cpp)
rm(sqbias.cpp, tt, K.tt.K, sqbias.R)

# Now that we know this works, we turn our attention to B1

tau.star <- tau.hat
Psi.M  <- Psi %*% M
n <- nrow(CGdata)
