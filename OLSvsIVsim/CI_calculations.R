n.reps <- 1000

rho.grid <- seq(from = 0, to = 0.1, by = 0.1)
pi.grid <- seq(from = 0.1, to = 0.2, by = 0.1)
sample.size.grid <- c(250, 500)

params <- expand.grid(n = sample.size.grid, p = pi.grid, r = rho.grid)

CI.results <- mapply(CIs_compare_default_cpp, 
                     p = params$p, r = params$r, n = params$n, 
                     n_reps = n.reps, SIMPLIFY = FALSE)

cover <- do.call(function(x) {x$coverage.prob}, CI.results)
i