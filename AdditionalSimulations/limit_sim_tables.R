setwd("~/fmsc/AdditionalSimulations/")
load("OLSvsIV_limit_sim.Rd")
load("chooseIVs_limit_sim.Rd")


#============================ Append Relative Width Columns

#-------------------------------- OLS vs IV
# For the OLS vs IV example, the width of the "valid" estimator
# is 2 * qnorm(1 - alpha/2) * sqrt(1 / pi_sq)
OLSvsIV$onestep_equal$widthvalid <- with(OLSvsIV$onestep_equal,
                                         2 * qnorm(1 - alpha/2) * sqrt(1/pi_sq))
OLSvsIV$onestep_equal$relwidth <- with(OLSvsIV$onestep_equal,
                                       avgwidth / widthvalid)

OLSvsIV$onestep_short$widthvalid <- with(OLSvsIV$onestep_short,
                                         2 * qnorm(1 - alpha/2) * sqrt(1/pi_sq))
OLSvsIV$onestep_short$relwidth <- with(OLSvsIV$onestep_short,
                                       avgwidth / widthvalid)

OLSvsIV$twostep_equal$widthvalid <- with(OLSvsIV$twostep_equal,
                                         2 * qnorm(1 - alpha/2) * sqrt(1/pi_sq))
OLSvsIV$twostep_equal$relwidth <- with(OLSvsIV$twostep_equal,
                                         avgwidth / widthvalid)

OLSvsIV$twostep_widetau$widthvalid <- with(OLSvsIV$twostep_widetau,
                                         2 * qnorm(1 - alpha/2) * sqrt(1/pi_sq))
OLSvsIV$twostep_widetau$relwidth <- with(OLSvsIV$twostep_widetau,
                                         avgwidth / widthvalid)

OLSvsIV$twostep_narrowtau$widthvalid <- with(OLSvsIV$twostep_narrowtau,
                                         2 * qnorm(1 - alpha/2) * sqrt(1/pi_sq))
OLSvsIV$twostep_narrowtau$relwidth <- with(OLSvsIV$twostep_narrowtau,
                                         avgwidth / widthvalid)


#--------------------------- Choose IVs
# For the choosing IVs example, the width of the "valid" estimator
# is 2 * qnorm(1 - alpha/2) * 3
chooseIVs$onestep_equal$widthvalid <- with(OLSvsIV$onestep_equal,
                                         2 * qnorm(1 - alpha/2) * 3)
chooseIVs$onestep_equal$relwidth <- with(OLSvsIV$onestep_equal,
                                       avgwidth / widthvalid)

chooseIVs$onestep_short$widthvalid <- with(OLSvsIV$onestep_short,
                                         2 * qnorm(1 - alpha/2) * 3)
chooseIVs$onestep_short$relwidth <- with(OLSvsIV$onestep_short,
                                       avgwidth / widthvalid)

chooseIVs$twostep_equal$widthvalid <- with(OLSvsIV$twostep_equal,
                                         2 * qnorm(1 - alpha/2) * 3)
chooseIVs$twostep_equal$relwidth <- with(OLSvsIV$twostep_equal,
                                         avgwidth / widthvalid)

chooseIVs$twostep_widetau$widthvalid <- with(OLSvsIV$twostep_widetau,
                                           2 * qnorm(1 - alpha/2) * 3)
chooseIVs$twostep_widetau$relwidth <- with(OLSvsIV$twostep_widetau,
                                         avgwidth / widthvalid)

chooseIVs$twostep_narrowtau$widthvalid <- with(OLSvsIV$twostep_narrowtau,
                                             2 * qnorm(1 - alpha/2) * 3)
chooseIVs$twostep_narrowtau$relwidth <- with(OLSvsIV$twostep_narrowtau,
                                         avgwidth / widthvalid)


#=========================== Coverage of Naive Intervals
xtabs(I(100 * round(coverprob,2)) ~ pi_sq + tau + alpha, OLSvsIV$cover_naive)
xtabs(I(100 * round(coverprob,2)) ~ g_sq + tau + alpha, chooseIVs$cover_naive)

#=========================== CI widths: naive and infeasible FMSC
xtabs(I(100 * round(erelwidth,2)) ~ pi_sq + tau, OLSvsIV$relwidth_naive)
xtabs(I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$relwidth_infeas)
xtabs(I(100 * round(erelwidth,2)) ~ g_sq + tau, chooseIVs$relwidth_naive)
xtabs(I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$relwidth_infeas)

#=========================== Coverage of One-step Equal-Tailed Intervals
xtabs(I(100 * round(coverage,2)) ~ pi_sq + tau + alpha, OLSvsIV$onestep_equal)
xtabs(I(100 * round(coverage,2)) ~ g_sq + tau + alpha, chooseIVs$onestep_equal)

#============================ Relative Width of One-step Equal-Tailed Intervals
xtabs(I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$onestep_equal)
xtabs(I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$onestep_equal)

#=========================== Coverage of One-step Shortest Intervals
xtabs(I(100 * round(coverage,2)) ~ pi_sq + tau + alpha, OLSvsIV$onestep_short)
xtabs(I(100 * round(coverage,2)) ~ g_sq + tau + alpha, chooseIVs$onestep_short)

#============================ Relative Width of One-step Shortest Intervals
xtabs(I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$onestep_short)
xtabs(I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$onestep_short)

#=========================== Coverage of Two-step Intervals a1 = a2
xtabs(I(100 * round(coverage,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_equal)
xtabs(I(100 * round(coverage,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_equal)

#============================ Relative Width of Two-step Intervals a1 = a2
xtabs(I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_equal)
xtabs(I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_equal)

#=========================== Coverage of Two-step Intervals a1 < a2
xtabs(I(100 * round(coverage,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_widetau)
xtabs(I(100 * round(coverage,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_widetau)

#============================ Relative Width of Two-step Intervals a1 < a2
xtabs(I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_widetau)
xtabs(I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_widetau)
