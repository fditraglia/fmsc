#Depends on fmsc_2sls_cpp, contained in functions.cpp

#----------------------------------------------------------------
#Set up data matrices to input to fmsc_2SLS_cpp()
#----------------------------------------------------------------
valid.instruments <- instrument.sets[[1]]
suspect.instruments <- setdiff(instrument.sets[[8]], valid.instruments)
y <- CGdata$lngdpc
constant <-  rep(1, nrow(CGdata))
x <- as.matrix(cbind(constant, CGdata[, c("rule", "malfal")]))
z1 <- as.matrix(cbind(constant, CGdata[, valid.instruments]))
z2 <- as.matrix(CGdata[, suspect.instruments])

#----------------------------------------------------------------
#Set up matrix of candidates to input to fmsc_2SLS_cpp()
#----------------------------------------------------------------
#    Each column is a vector of ones and zeros that indexes the
#    columns of z2 included in that specification. The valid and
#    full specifications added automatically by fmsc_2SLS_cpp()
#----------------------------------------------------------------
z2.cols <- lapply(instrument.sets, match, suspect.instruments)
z2.cols <- lapply(z2.cols, function(x) x[!is.na(x)])
f <- function(x){
  out <- rep(0, ncol(z2))
  out[x] <- 1
  return(out)
}
z2.ind <- sapply(lapply(z2.cols, f), rbind)
candidates <- z2.ind[,-c(1,8)]

#----------------------------------------------------------------
#Calculate the basic quantities needed for the FMSC
#----------------------------------------------------------------
fmsc.ingredients <- fmsc_2SLS_cpp(x, y, z1, z2, candidates)
rownames(fmsc.ingredients$est) <- colnames(x)

#----------------------------------------------------------------
#Clean Up
#----------------------------------------------------------------
rm(z2.ind, z2.cols, candidates, x, y, z1, z2, constant)
rm(valid.instruments, suspect.instruments, f)