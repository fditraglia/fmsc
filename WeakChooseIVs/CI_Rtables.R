#Load and unpack results
load("./Results/CI_results.Rdata")
coverage.prob <- results$coverage.prob
median.width <- results$median.width
rm(results)

params <- coverage.prob[,c("n", "g", "r")]
keep.cover <- c("Valid", "FMSC_naive", "FMSC_2step", "posFMSC_2step")
keep.width <- setdiff(keep.cover, "FMSC_naive")

coverage.prob <- coverage.prob[,keep.cover]
coverage.prob <- round(100 * coverage.prob)
coverage.prob <- cbind(params, coverage.prob)

median.width <- median.width[,keep.width]
width.rel.valid <- (median.width[,setdiff(keep.width, "Valid")] - median.width$Valid) / median.width$Valid
width.rel.valid <- round(100 * width.rel.valid)
width.rel.valid <- cbind(params, width.rel.valid)

#Clean up
rm(params, median.width)

wide.g.by.r <- function(long.fix.n){
  n <- unique(long.fix.n$n)
  data <- long.fix.n[,setdiff(names(long.fix.n), "n")]
  out <- reshape(data, idvar = "g", timevar = "r", , direction = "wide")
  names(out) <- c(paste("N =", n), unique(long.fix.n$r))
  return(out)
} 


n.by.g.by.r <- function(dataframe, colName){
  data <- dataframe[,c("n", "g", "r", colName)]
  lapply(unique(data$n), function(x) wide.g.by.r(subset(data, n == x)))
}

cover.list <- lapply(keep.cover, function(colName) n.by.g.by.r(coverage.prob, colName))
names(cover.list) <- keep.cover

width.list <- lapply(setdiff(keep.width, "Valid"), function(colName) n.by.g.by.r(width.rel.valid, colName))
names(width.list) <- setdiff(keep.width, "Valid")

rm(width.rel.valid, keep.cover, keep.width, coverage.prob)
rm(n.by.g.by.r)
rm(wide.g.by.r)

