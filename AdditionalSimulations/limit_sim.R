#======================= OLS vs IV Example
pi_sq_seq <- seq(0.1, 0.4, 0.05)
tau_seq <- 0:8
alpha_seq <- c(0.05, 0.1)
params <- expand.grid(tau = tau_seq, pi_sq = pi_sq_seq, alpha = alpha_seq)

#-------- Size Distortion of Naive Intervals
cover_naive_OLSvsIV <- unlist(Map(function(alpha, tau, pi_sq)
  fmscr::cover_naive(alpha, tau, bias_coef = 1,
                     tau_sd = sqrt((1 - pi_sq)/pi_sq), efficient_sd = 1),
  params$alpha, params$tau, params$pi_sq))

cover_naive_OLSvsIV <- data.frame(params, coverprob = cover_naive_OLSvsIV)
xtabs(I(100 * round(coverprob, 2)) ~ pi_sq + tau + alpha, cover_naive_OLSvsIV)


#-------- Expected Relative Width: Naive vs. Valid
params <- expand.grid(tau = tau_seq, pi_sq = pi_sq_seq)

relwidth_naive_OLSvsIV <- unlist(Map(function(pi_sq, tau)
  expect_rel_width(tau, bias_coef = 1,
                   tau_sd = sqrt((1 - pi_sq)/pi_sq),
                   efficient_sd = 1),
  params$pi_sq, params$tau))

relwidth_naive_OLSvsIV <- data.frame(params, erelwidth = relwidth_naive_OLSvsIV)
xtabs(I(100 * round(erelwidth, 2)) ~ pi_sq + tau, relwidth_naive_OLSvsIV)



#======================= Choosing IVs Example
g_sq_seq <- seq(0.05, 0.35, 0.05)
tau_seq <- 0:8
alpha_seq <- c(0.05, 0.1)
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq, alpha = alpha_seq)

#-------- Size Distortion of Naive Intervals
cover_naive_chooseIV <- unlist(Map(function(alpha, tau, g_sq)
  fmscr::cover_naive(alpha, tau,
                     bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                     tau_sd = sqrt(1 + 9 * g_sq),
                     efficient_sd = sqrt(1 / (g_sq + 1/9))),
  params$alpha, params$tau, params$g_sq))

cover_naive_chooseIV <- data.frame(params, coverprob = cover_naive_chooseIV)
xtabs(I(100 * round(coverprob, 2)) ~ g_sq + tau + alpha, cover_naive_chooseIV)

#-------- Expected Relative Width: Naive vs. Valid
params <- expand.grid(tau = tau_seq, g_sq = g_sq_seq)

relwidth_naive_chooseIV <- unlist(Map(function(g_sq, tau)
  expect_rel_width(tau, bias_coef = sqrt(g_sq)/(g_sq + 1/9),
                   tau_sd = sqrt(1 + 9 * g_sq),
                   efficient_sd = sqrt(1 / (g_sq + 1/9))),
  params$g_sq, params$tau))

relwidth_naive_chooseIV <- data.frame(params, erelwidth=relwidth_naive_chooseIV)
xtabs(I(100 * round(erelwidth, 2)) ~ g_sq + tau, relwidth_naive_chooseIV)


