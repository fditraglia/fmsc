load("OLSvsIV_limitsim90.Rd")
load("OLSvsIV_limitsim95.Rd")

#======================== Coverage Tables (Actual - Nominal)
coverage_tables <- list(
  xtabs(I(100 * round(naiveCover, 2) - 90) ~ pi_sq + tau, limitsim90),
  xtabs(I(100 * round(fmsc1Cover, 2) - 90) ~ pi_sq + tau, limitsim90),
  xtabs(I(100 * round(avg1Cover, 2) - 90) ~ pi_sq + tau, limitsim90),
  xtabs(I(100 * round(naiveCover, 2) - 95) ~ pi_sq + tau, limitsim95),
  xtabs(I(100 * round(fmsc1Cover, 2) - 95) ~ pi_sq + tau, limitsim95),
  xtabs(I(100 * round(avg1Cover, 2) - 95) ~ pi_sq + tau, limitsim95))
names(coverage_tables) <- c("Naive 90%", "FMSC1 90%", "AVG1 90%",
                            "Naive 95%", "FMSC1 95%", "AVG1 95%")

#======================== Width Tables (Median Relative to TSLS)
width_tables <- list(
  xtabs(I(100 * round(naiveWidth/tslsWidth, 2)) ~ pi_sq + tau, limitsim95),
  xtabs(I(100 * round(fmsc1Width/tslsWidth, 2)) ~ pi_sq + tau, limitsim95),
  xtabs(I(100 * round(avg1Width/tslsWidth, 2)) ~ pi_sq + tau, limitsim95),
  xtabs(I(100 * round(naiveWidth/tslsWidth, 2)) ~ pi_sq + tau, limitsim95),
  xtabs(I(100 * round(fmsc1Width/tslsWidth, 2)) ~ pi_sq + tau, limitsim95),
  xtabs(I(100 * round(avg1Width/tslsWidth, 2)) ~ pi_sq + tau, limitsim95))
names(width_tables) <- c("Naive 90%", "FMSC1 90%", "AVG1 90%",
                          "Naive 95%", "FMSC1 95%", "AVG1 95%")


