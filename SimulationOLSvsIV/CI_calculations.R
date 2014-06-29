n.reps <- 5000

rho.grid <- seq(from = 0, to = 0.5, by = 0.1)
pi.grid <- seq(from = 0.1, to = 0.6, by = 0.1)
sample.size.grid <- c(250, 500, 1000)

params <- expand.grid(n = sample.size.grid, p = pi.grid, r = rho.grid)

CI.results <- mapply(CIs_compare_default_cpp, 
                     p = params$p, r = params$r, n = params$n, 
                     n_reps = n.reps, SIMPLIFY = FALSE)

coverage.prob <- cbind(params, t(sapply(CI.results, function(x) x$coverage.prob)))
median.width <- cbind(params, t(sapply(CI.results, function(x) x$median.width)))

results <- list(coverage.prob = coverage.prob, median.width = median.width)

save(results, file = "CI_results.Rdata")

#Clean up
rm(n.reps, rho.grid, pi.grid, sample.size.grid)
rm(params, CI.results, coverage.prob, median.width)
rm(results)
