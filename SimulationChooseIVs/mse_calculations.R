n.reps <- 10
gamma.fine <- seq(0, 0.4, 0.05)
gamma.coarse <- c(0, 0.2, 0.4) 
rho.fine <- seq(0, 0.4, 0.05)
rho.coarse <- c(0, 0.2, 0.4)
n.grid <- c(250, 500, 1000)  

params.gamma.coarse <- expand.grid(n = n.grid,
                                   g = gamma.coarse,
                                   r = rho.fine)
params.rho.coarse <- expand.grid(n = n.grid,
                                 g = gamma.fine,
                                 r = rho.coarse)

results.rho.coarse <- mcmapply(mse_compare_default_cpp,
                               g = params.rho.coarse$g,
                               r = params.rho.coarse$r,
                               n = params.rho.coarse$n,
                               n_reps = n.reps,
                               mc.cores = nCores)
results.rho.coarse <- cbind(params.rho.coarse, 
                            t(results.rho.coarse))

results.gamma.coarse <- mcmapply(mse_compare_default_cpp,
                               g = params.gamma.coarse$g,
                               r = params.gamma.coarse$r,
                               n = params.gamma.coarse$n,
                               n_reps = n.reps,
                               mc.cores = nCores)
results.gamma.coarse <- cbind(params.gamma.coarse, 
                            t(results.gamma.coarse))

#Output results as R object
results <- list(coarse.rho = results.rho.coarse,
                coarse.gamma = results.gamma.coarse)
save(results, file = "./Results/mse_results.Rdata")

#Clean up
rm(n.reps, n.grid)
rm(gamma.fine, gamma.coarse)
rm(rho.fine, rho.coarse)
rm(params.gamma.coarse, params.rho.coarse)
rm(results.rho.coarse, results.gamma.coarse)
rm(results)
