n.reps <- 20000
gamma.coarse <- c(0.2, 0.3, 0.4) 
rho.fine <- seq(0, 0.3, 0.01)
n.grid <- c(50, 100, 500)  

params.gamma.coarse <- expand.grid(n = n.grid,
                                   g = gamma.coarse,
                                   r = rho.fine)

results.gamma.coarse <- mcmapply(mse_compare_default_cpp,
                               g = params.gamma.coarse$g,
                               r = params.gamma.coarse$r,
                               n = params.gamma.coarse$n,
                               n_reps = n.reps,
                               mc.cores = nCores)

results.gamma.coarse <- cbind(params.gamma.coarse, 
                            t(results.gamma.coarse))

#Output results as R object
results <- list(coarse.rho = NULL,
                coarse.gamma = results.gamma.coarse)
save(results, file = "./Results/mse_results.Rdata")

#Clean up
rm(n.reps, n.grid)
rm(gamma.coarse)
rm(rho.fine)
rm(params.gamma.coarse)
rm(results.gamma.coarse)
rm(results)
