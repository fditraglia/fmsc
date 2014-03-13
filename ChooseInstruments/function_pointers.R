setwd("~/fmsc/ChooseInstruments")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("function_pointers.cpp")


#Simple example with scalars
foo(1, 1) #Passes Plus to ResultPlusOne -> a + b + 1
bar(1, 1) #Passes Minus to ResultPlusOne -> a - b + 1

#More complicated example with matrix inner products
k <- 4
b <- 1:k / 10
V <- matrix(rnorm(k^2), k, k)
baz1(b, V) - V[1,1]
baz2(b, V) - V[2,2]
baz3(b, V) -  (t(b^2) %*% V %*% b^2)

#Passing a function pointer *through* one function and *into* another
PassThrough(1, 1) #Should give same result as foo
foo(1, 1)
