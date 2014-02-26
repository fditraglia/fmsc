#This script tests the R's built-in function for canonical correlation analysis, and creates a "minimialist" version to port to C++


#Minimialist R code following the source of cancor and using the defaults. Basically I've just removed all the error-checking, renamed some variables, and made the code clearer at the expense of speed by adding more temporaries than are strictly necessary.
cancor_r <- function(x, y){
  n_obs <- nrow(x)
  n_x <- ncol(x)
  n_y <- ncol(y)
  x <- x - rep(colMeans(x), rep(n_obs, n_x))
  y <- y - rep(colMeans(y), rep(n_obs, n_y))
  QR_x <- qr(x)
  QR_y <- qr(y)
  dx <- QR_x$rank
  dy <- QR_y$rank
  
  Qx <- qr.Q(QR_x)
  Rx <- qr.R(QR_x)
  Qy <- qr.Q(QR_y)
  Ry <- qr.R(QR_y) 
  
  z <- svd(t(Qx) %*% Qy, dx, dy)
  xcoef <- backsolve(Rx, z$u)
  ycoef <- backsolve(Ry, z$v)
  cor <- z$d
  out <- list(cor, xcoef = xcoef, ycoef = ycoef)
  return(out)
}


#Check if the results are the same using the example from the help file for cancor
pop <- LifeCycleSavings[, 2:3]
oec <- LifeCycleSavings[, -(2:3)]
foo <- cancor(pop, oec)[1:3]
bar <- cancor_r(pop, oec)
foo
bar
all.equal(foo, bar, check.attributes = FALSE)

