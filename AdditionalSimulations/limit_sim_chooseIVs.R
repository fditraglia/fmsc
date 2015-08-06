#======================== Choosing IVs Simulation Example
set.seed(8372)
n_cores <- parallel::detectCores()

alpha_seq <- c(0.2, 0.1, 0.05)
tau_seq <- seq(0, 5, 1)
g_sq_seq <- seq(0.1, 0.4, 0.1)


#======================== Size Distortion of Naive Intervals
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq, alpha = alpha_seq)

cover_naive <- unlist(Map(function(alpha, tau, g_sq)
  fmscr::cover_naive(alpha, tau,
                     bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                     tau_sd = sqrt(1 + 9 * g_sq),
                     efficient_sd = sqrt(1 / (g_sq + 1/9))),
  params$alpha, params$tau, params$g_sq))

cover_naive <- data.frame(params, coverprob = cover_naive)

#========================= Expected Relative Width: Naive vs. Valid
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq)

relwidth_naive <- unlist(Map(function(g_sq, tau)
  fmscr::expect_rel_width(tau, bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                   tau_sd = sqrt(1 + 9 * g_sq),
                   efficient_sd = sqrt(1 / (g_sq + 1/9))),
  params$g_sq, params$tau))

relwidth_naive <- data.frame(params, erelwidth=relwidth_naive)


#========================= Relative Width of Infeasible post-FMSC CI
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq, alpha = alpha_seq)
relwidth_infeas <- unlist(parallel::mcMap(function(alpha, tau, g_sq)
  fmscr::rel_width_FMSCinfeas(alpha, tau,
                              bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                              tau_sd = sqrt(1 + 9 * g_sq),
                              efficient_sd = sqrt(1 / (g_sq + 1/9)),
                              equal.tailed = FALSE),
  params$alpha, params$tau, params$g_sq, mc.cores = n_cores))

relwidth_infeas <- data.frame(params, relwidth = relwidth_infeas)

#========================= One-step CI: Equal-tailed
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq, alpha = alpha_seq)
onestep_equal <- parallel::mcMap(function(alpha, tau, g_sq)
  fmscr::onestepCI(alpha, tau,
                   bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                   tau_sd = sqrt(1 + 9 * g_sq),
                   efficient_sd = sqrt(1 / (g_sq + 1/9)),
                   equal.tailed = TRUE, ndraws = 5),
  params$alpha, params$tau, params$g_sq, mc.cores = n_cores)

onestep_equal <- data.frame(params, do.call(rbind, onestep_equal))

#========================= One-step CI: Shortest
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq, alpha = alpha_seq)
onestep_short <- parallel::mcMap(function(alpha, tau, g_sq)
  fmscr::onestepCI(alpha, tau,
                   bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                   tau_sd = sqrt(1 + 9 * g_sq),
                   efficient_sd = sqrt(1 / (g_sq + 1/9)),
                   equal.tailed = FALSE, ndraws = 5),
  params$alpha, params$tau, params$g_sq, mc.cores = n_cores)

onestep_short <- data.frame(params, do.call(rbind, onestep_short))

#========================= Two-step CI: Equal
#-------- a1 = 0.5 * alpha
#-------- a2 = 0.5 * alpha
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq, alpha = alpha_seq)
twostep_equal <- parallel::mcMap(function(alpha, tau, g_sq)
  fmscr::twostepCI(alpha, tau,
                   bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                   tau_sd = sqrt(1 + 9 * g_sq),
                   efficient_sd = sqrt(1 / (g_sq + 1/9)),
                   a1 = 0.5 * alpha, ndraws = 5),
  params$alpha, params$tau, params$g_sq, mc.cores = n_cores)

twostep_equal <- data.frame(params, do.call(rbind, twostep_equal))


#========================= Two-step CI: Wide Tau
#-------- a1 = 0.25 * alpha
#-------- a2 = 0.75 * alpha
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq, alpha = alpha_seq)
twostep_widetau <- parallel::mcMap(function(alpha, tau, g_sq)
  fmscr::twostepCI(alpha, tau,
                   bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                   tau_sd = sqrt(1 + 9 * g_sq),
                   efficient_sd = sqrt(1 / (g_sq + 1/9)),
                   a1 = 0.25 * alpha, ndraws = 5),
  params$alpha, params$tau, params$g_sq, mc.cores = n_cores)

twostep_widetau <- data.frame(params, do.call(rbind, twostep_widetau))

#========================= Two-step CI: Narrow Tau
#-------- a1 = 0.75 * alpha
#-------- a2 = 0.25 * alpha
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq, alpha = alpha_seq)
twostep_narrowtau <- parallel::mcMap(function(alpha, tau, g_sq)
  fmscr::twostepCI(alpha, tau,
                   bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                   tau_sd = sqrt(1 + 9 * g_sq),
                   efficient_sd = sqrt(1 / (g_sq + 1/9)),
                   a1 = 0.75 * alpha, ndraws = 5),
  params$alpha, params$tau, params$g_sq, mc.cores = n_cores)

twostep_narrowtau <- data.frame(params, do.call(rbind, twostep_narrowtau))

#========================== Store Results
chooseIVs <- list(cover_naive, relwidth_naive, relwidth_infeas,
                  onestep_equal, onestep_short, twostep_equal, twostep_widetau,
                  twostep_narrowtau)
names(chooseIVs) <- c('cover_naive', 'relwidth_naive', 'relwidth_infeas',
                  'onestep_equal', 'onestep_short', 'twostep_equal',
                  'twostep_widetau', 'twostep_narrowtau')

save(chooseIVs, file = "chooseIVs_limit_sim.Rd")

#========================== Tables of Results
#xtabs(I(100 * round(coverprob, 2)) ~ g_sq + tau + alpha, cover_naive)
#xtabs(I(100 * round(erelwidth, 2)) ~ g_sq + tau, relwidth_naive)
#xtabs(I(100 * round(relwidth, 2)) ~ g_sq + tau + alpha, relwidth_infeas)


