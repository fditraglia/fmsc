n.reps <- 10000
rho.fine <- seq(0, 0.6, 0.01)
pi.coarse <- c(0.01, 0.05, 0.1)
n.grid<- c(50, 100, 500)

params.pi.coarse <- expand.grid(n = n.grid,
                                p = pi.coarse, 
                                r = rho.fine)

results.pi.coarse <- mcmapply(mse_compare_default_cpp, 
                               p = params.pi.coarse$p,
                               r = params.pi.coarse$r,
                               n = params.pi.coarse$n,
                               n_reps = n.reps,
                               mc.cores = nCores)

results.pi.coarse <- cbind(params.pi.coarse, 
                           t(results.pi.coarse))

#Output results as R object
results <- list(coarse.rho = NULL,
                coarse.pi = results.pi.coarse)
save(results, file = "./Results/mse_results.Rdata")

#Clean up
rm(n.reps, n.grid)
rm(pi.coarse)
rm(rho.fine)
rm(params.pi.coarse)
rm(results.pi.coarse)
rm(results)