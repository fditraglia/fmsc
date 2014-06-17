n.reps <- 1000
rho.fine <- seq(0, 0.6, 0.1)
pi.fine <- seq(0.1, 0.6, 0.1)
rho.coarse <- seq(0.1, 0.2, 0.1)
pi.coarse <- seq(0.1, 0.2, 0.1)

sample.size.grid <- c(50, 100, 250)

#Run the simulation over all values in rho.fine for a fixed value of pi and sample size
fix.pi <- function(pi.value, sample.size){
  out <- lapply(rho.fine, function(r) mse_compare_default_cpp(pi.value, r, sample.size, n.reps))
  out <- do.call(rbind, out) #Convert list to matrix
  nRows <- length(rho.fine)
  cbind(n = rep(sample.size, nRows), p = rep(pi.value, nRows), r = rho.fine, out)
}

#Run the simulation over all values in pi.fine for a fixed value of rho and sample size
fix.rho <- function(rho.value, sample.size){
  out <- lapply(pi.fine, function(p) mse_compare_default_cpp(p, rho.value, sample.size, n.reps))
  out <- do.call(rbind, out) #Convert list to matrix
  nRows <- length(pi.fine)
  cbind(n = rep(sample.size, nRows), p = pi.fine, r = rep(rho.value, nRows), out)
}

#Run fix.pi over all values in pi.coarse for a fixed sample size
panels.coarse.pi <- function(sample.size){
  lapply(pi.coarse, function(p) fix.pi(p, sample.size))
} 

#Run fix.rho over all values in rho.coarse for a fixed sample size
panels.coarse.rho <- function(sample.size){
  lapply(rho.coarse, function(r) fix.rho(r, sample.size))
}

results.coarse.rho <- lapply(sample.size.grid, panels.coarse.rho)
results.coarse.pi <- lapply(sample.size.grid, panels.coarse.pi)

#Output results as R object
results <- list(coarse.rho = results.coarse.rho,
                coarse.pi = results.coarse.pi)
save(results, file = "mse_results.Rdata")

#Clean up
rm(n.reps, sample.size.grid)
rm(pi.fine, pi.coarse)
rm(rho.fine, rho.coarse)
rm(fix.pi, fix.rho)
rm(panels.coarse.pi, panels.coarse.rho)
rm(results.coarse.pi, results.coarse.rho)
rm(results)