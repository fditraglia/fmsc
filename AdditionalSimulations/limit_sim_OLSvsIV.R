#========================= OLS vs TSLS Simulation Example
set.seed(8372)
n_cores <- parallel::detectCores()

alpha_seq <- c(0.2, 0.1, 0.05)
pi_sq_seq <- seq(0.1, 0.4, 0.1)
tau_seq <- seq(0, 5, 1)

params <- expand.grid(tau = tau_seq, pi_sq = pi_sq_seq, alpha = alpha_seq)

#========================= Size Distortion of Naive Intervals
cover_naive <- unlist(Map(function(alpha, tau, pi_sq)
  fmscr::cover_naive(alpha, tau, bias_coef = 1,
                     tau_sd = sqrt((1 - pi_sq)/pi_sq), efficient_sd = 1),
  params$alpha, params$tau, params$pi_sq))

cover_naive <- data.frame(params, coverprob = cover_naive)


#========================= Expected Relative Width: Naive vs. Valid
params <- expand.grid(tau = tau_seq, pi_sq = pi_sq_seq)

relwidth_naive <- unlist(Map(function(pi_sq, tau)
  fmscr::expect_rel_width(tau, bias_coef = 1,
                   tau_sd = sqrt((1 - pi_sq)/pi_sq),
                   efficient_sd = 1),
  params$pi_sq, params$tau))

relwidth_naive <- data.frame(params, erelwidth = relwidth_naive)

#========================= Relative Width of Infeasible post-FMSC CI
params <- expand.grid(tau = tau_seq, pi_sq = pi_sq_seq, alpha = alpha_seq)
relwidth_infeas <- unlist(parallel::mcMap(function(alpha, tau, pi_sq)
  fmscr::rel_width_FMSCinfeas(alpha, tau, bias_coef = 1,
                              tau_sd = sqrt((1 - pi_sq)/pi_sq),
                              efficient_sd = 1, equal.tailed = FALSE),
  params$alpha, params$tau, params$pi_sq, mc.cores = n_cores))

relwidth_infeas <- data.frame(params, relwidth = relwidth_infeas)


#========================== Tables of Results
xtabs(I(100 * round(coverprob, 2)) ~ pi_sq + tau + alpha, cover_naive)
xtabs(I(100 * round(erelwidth, 2)) ~ pi_sq + tau, relwidth_naive)
xtabs(I(100 * round(relwidth, 2)) ~ pi_sq + tau + alpha, relwidth_infeas)


