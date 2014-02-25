#This script tests the R's built-in function for canonical correlation analysis, and creates a "minimialist" version to port to C++

#The example from ?cancor
pop <- LifeCycleSavings[, 2:3]
oec <- LifeCycleSavings[, -(2:3)]
cancor(pop, oec)


#Minimialist R code following the source of cancor and using the defaults. Basically I've just removed all the error-checking
x <- as.matrix(pop)
y <- as.matrix(oec)
nr <- nrow(x)
ncx <- ncol(x)
ncy <- ncol(y)
xcenter <- colMeans(x, )
x <- x - rep(xcenter, rep(nr, ncx))
ycenter <- colMeans(y)
y <- y - rep(ycenter, rep(nr, ncy))
qx <- qr(x)
qy <- qr(y)
dx <- qx$rank
dy <- qy$rank
z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1:dx,], dx, dy)
xcoef <- backsolve((qx$qr)[1:dx, 1:dx], z$u)
ycoef <- backsolve((qy$qr)[1:dy, 1:dy], z$v)
cor <- z$d

ycoef - cancor(pop, oec)$ycoef
xcoef - cancor(pop, oec)$xcoef
all.equal(cor, cancor(pop, oec)$cor)