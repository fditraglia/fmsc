#This script tests the R's built-in function for canonical correlation analysis, and creates a "minimialist" version to port to C++


#Minimialist R code following the source of cancor and using the defaults. Basically I've just removed all the error-checking and renamed some variables.
cancor_r <- function(x, y){
  n_obs <- nrow(x)
  n_x <- ncol(x)
  n_y <- ncol(y)
  x <- x - rep(colMeans(x), rep(n_obs, n_x))
  y <- y - rep(colMeans(y), rep(n_obs, n_y))
  qx <- qr(x)
  qy <- qr(y)
  dx <- qx$rank
  dy <- qy$rank
  z <- svd(qr.qty(qx, qr.qy(qy, diag(1, n, dy)))[1:dx,], dx, dy)
  xcoef <- backsolve((qx$qr)[1:dx, 1:dx], z$u)
  ycoef <- backsolve((qy$qr)[1:dy, 1:dy], z$v)
  cor <- z$d
  out <- list(xcoef = xcoef, ycoef = ycoef, cor = cor)
  return(out)
}


#Check if the results are the same using the example from the help file for cancor
pop <- LifeCycleSavings[, 2:3]
oec <- LifeCycleSavings[, -(2:3)]
foo <- cancor(pop, oec)
bar <- cancor_r(pop, oec)
foo$ycoef - bar$ycoef
foo$xcoef - bar$xcoef
foo$cor - bar$cor
