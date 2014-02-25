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
x <- x - rep(xcenter, rep.int(nr, ncx))
ycenter <- colMeans(y)
y <- y - rep(ycenter, rep.int(nr, ncy))
qx <- qr(x)
qy <- qr(y)
dx <- qx$rank
dy <- qy$rank
z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx,,drop = FALSE], dx, dy)
xcoef <- backsolve((qx$qr)[1L:dx, 1L:dx, drop = FALSE], z$u)
rownames(xcoef) <- colnames(x)[qx$pivot][1L:dx]
ycoef <- backsolve((qy$qr)[1L:dy, 1L:dy, drop = FALSE], z$v)
rownames(ycoef) <- colnames(y)[qy$pivot][1L:dy]
cor <- z$d
