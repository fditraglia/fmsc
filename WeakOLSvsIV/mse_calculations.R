n.reps <- 10000
rho.fine <- seq(0, 0.6, 0.01)
pi.fine <- seq(0.1, 0.6, 0.01)
rho.coarse <- c(0, 0.1, 0.2)
pi.coarse <- c(0.2, 0.4, 0.6)
n.grid<- c(50, 100, 500)

params.pi.coarse <- expand.grid(n = n.grid,
                                p = pi.coarse, 
                                r = rho.fine)
params.rho.coarse <- expand.grid(n = n.grid, 
                                 p = pi.fine, 
                                 r = rho.coarse)

results.rho.coarse <- mcmapply(mse_compare_default_cpp, 
                                p = params.rho.coarse$p,
                                r = params.rho.coarse$r,
                                n = params.rho.coarse$n,
                                n_reps = n.reps,
                                mc.cores = nCores)
results.rho.coarse <- cbind(params.rho.coarse, 
                            t(results.rho.coarse))

results.pi.coarse <- mcmapply(mse_compare_default_cpp, 
                               p = params.pi.coarse$p,
                               r = params.pi.coarse$r,
                               n = params.pi.coarse$n,
                               n_reps = n.reps,
                               mc.cores = nCores)
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