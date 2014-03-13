setwd("~/fmsc/ChooseInstruments")
library(Rcpp)
sourceCpp("function_pointers.cpp")

foo(1, 1) #Passes Plus to ResultPlusOne -> a + b + 1
bar(1, 1) #Passes Minus to ResultPlusOne -> a - b + 1
