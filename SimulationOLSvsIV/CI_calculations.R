#Define simulation parameters
n.reps <- 10000
rho.grid <- seq(from = 0, to = 0.5, by = 0.1)
pi.grid <- seq(from = 0.1, to = 0.6, by = 0.1)
n.grid <- c(50, 100, 500)

params <- expand.grid(n = n.grid, p = pi.grid, r = rho.grid)

CI.results <- mcmapply(CIs_compare_default_cpp, 
                     p = params$p, r = params$r, n = params$n, 
                     n_reps = n.reps, SIMPLIFY = FALSE,
                     mc.cores = nCores)

coverage.prob <- cbind(params, t(sapply(CI.results, function(x) x$coverage.prob)))
median.width <- cbind(params, t(sapply(CI.results, function(x) x$median.width)))

results <- list(coverage.prob = coverage.prob, median.width = median.width)

save(results, file = "./Results/CI_results.Rdata")

#Clean up
rm(n.reps, rho.grid, pi.grid, n.grid)
rm(params, CI.results, coverage.prob, median.width)
rm(results)
