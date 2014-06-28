# reps <- 20000
# n.grid <- 500
# r.grid <- c(0, 0.2, 0.4)
# g.grid <- c(0, 0.2, 0.4)
# params <- expand.grid(g.grid, r.grid, n.grid)
# names(params) <- c("g", "r", "n")
# 
# system.time(out <- t(mapply(mse_compare_default_cpp, g = params$g, r = params$r, n = params$n, n_reps = reps)))
# out <- cbind(params, out)
# save(out, file = "rmse_results.Rdata")

