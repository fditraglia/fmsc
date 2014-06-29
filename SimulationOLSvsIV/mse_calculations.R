n.reps <- 100
rho.fine <- seq(0, 0.6, 0.1)
pi.fine <- seq(0.1, 0.6, 0.1)
rho.coarse <- c(0, 0.25, 0.5)
pi.coarse <- c(0.1, 0.3, 0.5)
n.grid<- c(50, 100, 250)

params.pi.coarse <- expand.grid(n.grid, pi.coarse, rho.fine)
params.rho.coarse <- expand.grid(n.grid, pi.fine, rho.coarse)
names(params.rho.coarse) <- names(params.pi.coarse) <- c("n", "p", "r")

results.rho.coarse <- mapply(mse_compare_default_cpp, 
                                p = params.rho.coarse$p,
                                r = params.rho.coarse$r,
                                n = params.rho.coarse$n,
                                n_reps = n.reps)
results.rho.coarse <- cbind(params.rho.coarse, 
                            t(results.rho.coarse))

results.pi.coarse <- mapply(mse_compare_default_cpp, 
                               p = params.pi.coarse$p,
                               r = params.pi.coarse$r,
                               n = params.pi.coarse$n,
                               n_reps = n.reps)
results.pi.coarse <- cbind(params.pi.coarse, 
                           t(results.pi.coarse))

#Output results as R object
results <- list(coarse.rho = results.rho.coarse,
                coarse.pi = results.pi.coarse)
save(results, file = "./Results/mse_results.Rdata")

#Clean up
rm(n.reps, n.grid)
rm(pi.fine, pi.coarse)
rm(rho.fine, rho.coarse)
rm(params.pi.coarse, params.rho.coarse)
rm(results.rho.coarse, results.pi.coarse)
rm(results)