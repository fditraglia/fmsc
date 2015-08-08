setwd("~/fmsc/AdditionalSimulations/")
rm(list = ls())

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
#=========================== Helper Function
make_TeXtables <- function(xtab_list, row_lab, col_lab){
  names_list <- paste0("\\", names(xtab_list))
  out <- Map(function(xtab, tab_name)
    fmscr::TeXtable(xtab, tab_name, row_lab, col_lab), xtab_list, names_list)
  out <- paste(unlist(out), collapse = "\n \n \\vspace{2em} \n \n")
  return(out)
}

#=========================== Coverage of Naive Intervals
c_naive_OLSvsIV <- fmscr::rtables(
  I(100 * round(coverprob,2)) ~ pi_sq + tau + alpha, OLSvsIV$cover_naive)
c_naive_OLSvsIV <- make_TeXtables(c_naive_OLSvsIV, "\\pi^2", "\\tau")

c_naive_chooseIV <- fmscr::rtables(
  I(100 * round(coverprob,2)) ~ g_sq + tau + alpha, chooseIVs$cover_naive)
c_naive_chooseIV <- make_TeXtables(c_naive_chooseIV, "\\gamma^2", "\\tau")


#=========================== CI widths: naive and infeasible FMSC
w_naive_OLSvsIV <- xtabs(
  I(100 * round(erelwidth,2)) ~ pi_sq + tau, OLSvsIV$relwidth_naive)
w_naive_OLSvsIV <- fmscr::TeXtable(w_naive_OLSvsIV, row_lab = "\\pi^2", 
                                   col_lab = "\\tau")

w_infeas_OLSvsIV <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$relwidth_infeas)
w_infeas_OLSvsIV <- make_TeXtables(w_infeas_OLSvsIV, "\\pi^2", "\\tau")


w_naive_chooseIVs <- xtabs(
  I(100 * round(erelwidth,2)) ~ g_sq + tau, chooseIVs$relwidth_naive)
w_naive_chooseIVs <- fmscr::TeXtable(w_naive_chooseIVs, row_lab = "\\pi^2", 
                                   col_lab = "\\tau") 

w_infeas_chooseIVs <- fmscr::rtables(I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$relwidth_infeas)
w_infeas_chooseIVs <- make_TeXtables(w_infeas_chooseIVs, "\\gamma^2", "\\tau")


#=========================== Coverage of One-step Equal-Tailed Intervals
c_1equal_OLSvsIV <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ pi_sq + tau + alpha, OLSvsIV$onestep_equal)
c_1equal_OLSvsIV <- make_TeXtables(c_1equal_OLSvsIV, "\\pi^2", "\\tau")

c_1equal_chooseIVs <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ g_sq + tau + alpha, chooseIVs$onestep_equal)
c_1equal_chooseIVs <- make_TeXtables(c_1equal_chooseIVs, "\\gamma^2", "\\tau")

#============================ Relative Width of One-step Equal-Tailed Intervals
w_1equal_OLSvsIV <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$onestep_equal)
w_1equal_OLSvsIV <- make_TeXtables(w_1equal_OLSvsIV, "\\pi^2", "\\tau")

w_1equal_chooseIVs <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$onestep_equal)
w_1equal_chooseIVs <- make_TeXtables(w_1equal_chooseIVs, "\\gamma^2", "\\tau")


#=========================== Coverage of One-step Shortest Intervals
c_1short_OLSvsIV <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ pi_sq + tau + alpha, OLSvsIV$onestep_short)
c_1short_OLSvsIV <- make_TeXtables(c_1short_OLSvsIV, "\\pi^2", "\\tau")

c_1short_chooseIVs <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ g_sq + tau + alpha, chooseIVs$onestep_short)
c_1short_chooseIVs <- make_TeXtables(c_1short_chooseIVs, "\\gamma^2", "\\tau")

#============================ Relative Width of One-step Shortest Intervals
w_1short_OLSvsIV <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$onestep_short)
w_1short_OLSvsIV <- make_TeXtables(w_1short_OLSvsIV, "\\pi^2", "\\tau")

w_1short_chooseIVs <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$onestep_short)
w_1short_chooseIVs <- make_TeXtables(w_1short_chooseIVs, "\\gamma^2", "\\tau")


#=========================== Coverage of Two-step Intervals a1 = a2
c_2equal_OLSvsIV <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_equal)
c_2equal_OLSvsIV <- make_TeXtables(c_2equal_OLSvsIV, "\\pi^2", "\\tau")

c_2equal_chooseIVs <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_equal)
c_2equal_chooseIVs <- make_TeXtables(c_2equal_chooseIVs, "\\gamma^2", "\\tau")

#============================ Relative Width of Two-step Intervals a1 = a2
w_2equal_OLSvsIV <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_equal)
w_2equal_OLSvsIV <- make_TeXtables(w_2equal_OLSvsIV, "\\pi^2", "\\tau")

w_2equal_chooseIVs <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_equal)
w_2equal_chooseIVs <- make_TeXtables(w_2equal_chooseIVs, "\\gamma^2", "\\tau")


#=========================== Coverage of Two-step Intervals a1 < a2
c_2widetau_OLSvsIV <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_widetau)
c_2widetau_OLSvsIV <- make_TeXtables(c_2widetau_OLSvsIV, "\\pi^2", "\\tau")

c_2widetau_chooseIVs <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_widetau)
c_2widetau_chooseIVs <- make_TeXtables(c_2widetau_chooseIVs, "\\gamma^2", "\\tau")

#============================ Relative Width of Two-step Intervals a1 < a2
w_2widetau_OLSvsIV <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_widetau)
w_2widetau_OLSvsIV <- make_TeXtables(w_2widetau_OLSvsIV, "\\pi^2", "\\tau")

w_2widetau_chooseIVs <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_widetau)
w_2widetau_chooseIVs <- make_TeXtables(w_2widetau_chooseIVs, "\\gamma^2", "\\tau")

#=========================== Coverage of Two-step Intervals a1 > a2
c_2narrowtau_OLSvsIV <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_narrowtau)
c_2narrowtau_OLSvsIV <- make_TeXtables(c_2narrowtau_OLSvsIV, "\\pi^2", "\\tau")

c_2narrowtau_chooseIVs <- fmscr::rtables(
  I(100 * round(coverage,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_narrowtau)
c_2narrowtau_chooseIVs <- make_TeXtables(c_2narrowtau_chooseIVs, "\\gamma^2", "\\tau")

#============================ Relative Width of Two-step Intervals a1 > a2
w_2narrowtau_OLSvsIV <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ pi_sq + tau + alpha, OLSvsIV$twostep_narrowtau)
w_2narrowtau_OLSvsIV <- make_TeXtables(w_2narrowtau_OLSvsIV, "\\pi^2", "\\tau")

w_2narrowtau_chooseIVs <- fmscr::rtables(
  I(100 * round(relwidth,2)) ~ g_sq + tau + alpha, chooseIVs$twostep_narrowtau)
w_2narrowtau_chooseIVs <- make_TeXtables(w_2narrowtau_chooseIVs, "\\gamma^2", "\\tau")

#=========================== Clean Up
rm(OLSvsIV, chooseIVs)
rm(make_TeXtables)

#=========================== Output Tables as TeX
setwd("~/fmsc/AdditionalSimulations/LimitSimResults/")
tables_list <- ls()

for(i in 1:length(tables_list)){
  temp_file <- paste0(tables_list[i], ".tex")
  out_table <- get(tables_list[i])
  cat(out_table, file = paste0("./", temp_file))
}  

rm(list = ls())