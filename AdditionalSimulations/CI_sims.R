set.seed(8372)
n_cores <- parallel::detectCores()

n_reps <- 1000
alpha_seq <- c(0.2, 0.1, 0.05)
pi_sq_seq <- seq(0.1, 0.4, 0.1)
N_seq <- c(50, 100, 500)

#========================= OLS vs TSLS Simulation Example
rho_seq <- seq(0, 0.5, 0.1)
params <- expand.grid(rho = rho_seq, pi_sq = pi_sq_seq, alpha = alpha_seq,
                      N = N_seq)

OLSvsIV_CIs <- parallel::mcMap(function(alpha, rho, pi_sq, N)
  fmscr::CIsim_OLSvsIV(alpha, rho, pi_sq, N, n_reps),
  params$alpha, params$rho, params$pi_sq, params$N)

OLSvsIV_CIs <- do.call(rbind, OLSvsIV_CIs)
OLSvsIV_CIs <- cbind(params, OLSvsIV_CIs)
save(OLSvsIV_CIs, file = "OLSvsIV_CIs.Rd")

#========================= Choose IVs Simulation
g_sq_seq <- seq(0.1, 0.4, 0.1)
params <- expand.grid(rho = rho_seq, g_sq = g_sq_seq, alpha = alpha_seq,
                      N = N_seq)

chooseIVs_CIs <- parallel::mcMap(function(alpha, rho, g_sq, N)
  fmscr::CIsim_chooseIVs(alpha, rho, g_sq, N, n_reps),
  params$alpha, params$rho, params$g_sq, params$N)

chooseIVs_CIs <- do.call(rbind, chooseIVs_CIs)
chooseIVs_CIs <- cbind(params, chooseIVs_CIs)
save(chooseIVs_CIs, file = "chooseIVs_CIs.Rd")

