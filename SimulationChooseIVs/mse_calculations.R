reps <- 20000
n.grid <- 500
r.grid <- c(0, 0.2, 0.4)
g.grid <- c(0, 0.2, 0.4)
params <- expand.grid(g.grid, r.grid, n.grid)
names(params) <- c("g", "r", "n")

out <- mcmapply(mse_compare_default_cpp, g = params$g, r = params$r, n = params$n, n_reps = reps, mc.cores = 5)
out <- cbind(params, t(out))
save(out, file = "rmse_results.Rdata")
rm(reps, n.grid, r.grid, g.grid, params, out)
