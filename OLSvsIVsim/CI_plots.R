#Load and unpack results
load("CI_results.Rdata")
coverage.prob <- results$coverage.prob
median.width <- results$median.width
rm(results)

params <- coverage.prob[,c("n", "p", "r")]
keep.cover <- c("TSLS", "FMSC_naive", "FMSC_correct", "AVG_correct")
keep.width <- setdiff(keep.cover, "FMSC_naive")

coverage.prob <- coverage.prob[,keep.cover]
coverage.prob <- round(100 * coverage.prob)
coverage.prob <- cbind(params, coverage.prob)

median.width <- median.width[,keep.width]
width.rel.tsls <- (median.width[,setdiff(keep.width, "TSLS")] - median.width$TSLS) / median.width$TSLS
width.rel.tsls <- round(100 * width.rel.tsls)
width.rel.tsls <- cbind(params, width.rel.tsls)

#Clean up
rm(params, median.width)

wide.p.by.r <- function(long.fix.n){
  n <- unique(long.fix.n$n)
  data <- long.fix.n[,setdiff(names(long.fix.n), "n")]
  out <- reshape(data, idvar = "p", timevar = "r", , direction = "wide")
  names(out) <- c(paste("N =", n), unique(long.fix.n$r))
  return(out)
} 

foo <- subset(width.rel.tsls, n == 500)[,c("n", "p", "r", "AVG_correct")]
wide.p.by.r(foo)

data <- coverage.prob[,c("n", "p", "r", "AVG_correct")]
lapply(unique(data$n), function(x) wide.p.by.r(subset(data, n == x)))

