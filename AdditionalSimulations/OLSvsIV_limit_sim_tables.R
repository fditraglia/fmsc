load("OLSvsIV_limitsim90.Rd")
load("OLSvsIV_limitsim95.Rd")

#======================= Replacement for %in%
# 0.3 on work machine != 0.3 on home machine!
# Code for match.numeric is adapted from StackOverflow:
# http://stackoverflow.com/questions/16309750/match-does-not-work
match.numeric <- function(x, table) {
  are.equal <- function(x, y) isTRUE(all.equal(x, y))
  match.one <- function(x, table)
    match(TRUE, vapply(table, are.equal, logical(1L), x = x))
  vapply(x, match.one, integer(1L), table)
}
`%nin%` <- function(x, table) !is.na(match.numeric(x, table))

#======================== Coarsen Simulation Grid for Tables
coarse90 <- subset(limitsim90, (tau %nin% 0:8) & (pi_sq %nin% (2:8/20)))
coarse95 <- subset(limitsim95, (tau %nin% 0:8) & (pi_sq %nin% (2:8/20)))

#======================== Coverage Tables (Actual - Nominal)
coverage_tables <- list(
  xtabs(I(100 * round(naiveCover, 2) - 90) ~ pi_sq + tau, coarse90),
  xtabs(I(100 * round(fmsc1Cover, 2) - 90) ~ pi_sq + tau, coarse90),
  xtabs(I(100 * round(avg1Cover, 2) - 90) ~ pi_sq + tau, coarse90),
  xtabs(I(100 * round(naiveCover, 2) - 95) ~ pi_sq + tau, coarse95),
  xtabs(I(100 * round(fmsc1Cover, 2) - 95) ~ pi_sq + tau, coarse95),
  xtabs(I(100 * round(avg1Cover, 2) - 95) ~ pi_sq + tau, coarse95))
names(coverage_tables) <- c("Naive 90%", "FMSC1 90%", "AVG1 90%",
                            "Naive 95%", "FMSC1 95%", "AVG1 95%")

#======================== Width Tables (Median Relative to TSLS)
width_tables <- list(
  xtabs(I(100 * round(naiveWidth/tslsWidth, 2)) ~ pi_sq + tau, coarse95),
  xtabs(I(100 * round(fmsc1Width/tslsWidth, 2)) ~ pi_sq + tau, coarse95),
  xtabs(I(100 * round(avg1Width/tslsWidth, 2)) ~ pi_sq + tau, coarse95),
  xtabs(I(100 * round(naiveWidth/tslsWidth, 2)) ~ pi_sq + tau, coarse95),
  xtabs(I(100 * round(fmsc1Width/tslsWidth, 2)) ~ pi_sq + tau, coarse95),
  xtabs(I(100 * round(avg1Width/tslsWidth, 2)) ~ pi_sq + tau, coarse95))
names(width_tables) <- c("Naive 90%", "FMSC1 90%", "AVG1 90%",
                          "Naive 95%", "FMSC1 95%", "AVG1 95%")


