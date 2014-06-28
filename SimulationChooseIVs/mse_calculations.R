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
set.seed(1234)
g <- 0.4
r <- 0.2
foo <- dgp_test(g, r, 50000)
var(with(foo, cbind(e,v,w))) - foo$V
all.equal(foo$y, foo$x * 0.5 + foo$e)
all.equal(foo$x, foo$z %*% rep(0.1, 3) + foo$w * g + foo$v)
