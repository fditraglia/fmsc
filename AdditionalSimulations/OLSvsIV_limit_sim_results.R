set.seed(8372)
n_cores <- parallel::detectCores()

tau_seq <- seq(0, 8)
pi_sq_seq <- seq(0.05, 0.4, 0.05)
params <- expand.grid(tau = tau_seq, pi_sq = pi_sq_seq)

#=================== OLS, TSLS and "Naive" CIs

nonsim90 <- Map(function(tau, pi_sq)
  fmscr::OLSvsIV_nonsimCI(tau, pi_sq, size = 0.1, n_sim = 10000),
  params$tau, params$pi_sq)

nonsim95 <- Map(function(tau, pi_sq)
  fmscr::OLSvsIV_nonsimCI(tau, pi_sq, size = 0.05, n_sim = 10000),
  params$tau, params$pi_sq)

#=================== One-step Simulation-based CIs

onestep90 <- parallel::mcMap(function(tau, pi_sq)
  fmscr::OLSvsIV_onestepCI(tau, pi_sq, size = 0.1, inc = 0.01,
                           n_sim_outer = 10000, n_sim_inner = 1000),
  params$tau, params$pi_sq, mc.cores = n_cores)

onestep95 <- parallel::mcMap(function(tau, pi_sq)
  fmscr::OLSvsIV_onestepCI(tau, pi_sq, size = 0.05, inc = 0.005,
                           n_sim_outer = 10000, n_sim_inner = 1000),
  params$tau, params$pi_sq, mc.cores = n_cores)

#==================== Two-step Simulation-based CIs
# Still need to write the C++ code
#   - Uniformly valid two-step CI
#   - Referee's Suggestion

#==================== Collect and Save Results

to_dataframe <- function(results_list){
  # Helper function to put the output of the C++ functions into
  # a more convenient form for later use
  do.call("rbind", lapply(results_list, as.data.frame))
}

limitsim90 <- cbind(params, to_dataframe(nonsim90), to_dataframe(onestep90))
limitsim95 <- cbind(params, to_dataframe(nonsim95), to_dataframe(onestep95))

save(limitsim90, file = "OLSvsIV_limitsim90.Rd")
save(limitsim95, file = "OLSvsIV_limitsim95.Rd")

rm(list = ls())
