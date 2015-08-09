setwd("~/fmsc/AdditionalSimulations/")
rm(list = ls())

#=========================== Helper Function
make_TeXtables <- function(xtab_list, row_lab, col_lab){
  names_list <- paste0("\\", names(xtab_list))
  out <- Map(function(xtab, tab_name)
    fmscr::TeXtable(xtab, tab_name, row_lab, col_lab), xtab_list, names_list)
  out <- paste(unlist(out), collapse = "\n \n \\vspace{2em} \n \n")
  return(out)
}

#============================ OLS vs TSLS Example
load("OLSvsIV_CIs.Rd")
OLSvsIV_CIs_50 <- subset(OLSvsIV_CIs, N == 50)
OLSvsIV_CIs_100 <- subset(OLSvsIV_CIs, N == 100)
OLSvsIV_CIs_500 <- subset(OLSvsIV_CIs, N == 500)


#--------------------- Coverage of TSLS Estimator
c_tsls_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(c_tsls,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_tsls_OLSvsIV_50 <- make_TeXtables(c_tsls_OLSvsIV_50, "\\pi^2", "\\rho")

c_tsls_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(c_tsls,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
c_tsls_OLSvsIV_100 <- make_TeXtables(c_tsls_OLSvsIV_100, "\\pi^2", "\\rho")

c_tsls_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(c_tsls,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_tsls_OLSvsIV_500 <- make_TeXtables(c_tsls_OLSvsIV_500, "\\pi^2", "\\rho")



#--------------------- Coverage of Naive CI
c_naive_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(c_naive,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_naive_OLSvsIV_50 <- make_TeXtables(c_naive_OLSvsIV_50, "\\pi^2", "\\rho")

c_naive_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(c_naive,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
c_naive_OLSvsIV_100 <- make_TeXtables(c_naive_OLSvsIV_100, "\\pi^2", "\\rho")

c_naive_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(c_naive,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_naive_OLSvsIV_500 <- make_TeXtables(c_naive_OLSvsIV_500, "\\pi^2", "\\rho")

#--------------------- Relative Width of Naive CI
w_naive_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(w_naive / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
w_naive_OLSvsIV_50 <- make_TeXtables(w_naive_OLSvsIV_50, "\\pi^2", "\\rho")

w_naive_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(w_naive / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
w_naive_OLSvsIV_100 <- make_TeXtables(w_naive_OLSvsIV_100, "\\pi^2", "\\rho")

w_naive_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(w_naive / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_500)
w_naive_OLSvsIV_500 <- make_TeXtables(w_naive_OLSvsIV_500, "\\pi^2", "\\rho")



#--------------------- Coverage of 1-Step Equal-Tailed CI
c_1equal_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(c_1equal,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_1equal_OLSvsIV_50 <- make_TeXtables(c_1equal_OLSvsIV_50, "\\pi^2", "\\rho")

c_1equal_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(c_1equal,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
c_1equal_OLSvsIV_100 <- make_TeXtables(c_1equal_OLSvsIV_100, "\\pi^2", "\\rho")

c_1equal_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(c_1equal,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_1equal_OLSvsIV_500 <- make_TeXtables(c_1equal_OLSvsIV_500, "\\pi^2", "\\rho")

#--------------------- Relative Width of 1-Step Equal-Tailed CI
w_1equal_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(w_1equal / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
w_1equal_OLSvsIV_50 <- make_TeXtables(w_1equal_OLSvsIV_50, "\\pi^2", "\\rho")

w_1equal_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(w_1equal / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
w_1equal_OLSvsIV_100 <- make_TeXtables(w_1equal_OLSvsIV_100, "\\pi^2", "\\rho")

w_1equal_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(w_1equal / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_500)
w_1equal_OLSvsIV_500 <- make_TeXtables(w_1equal_OLSvsIV_500, "\\pi^2", "\\rho")



#--------------------- Coverage of 1-Step Shortest CI
c_1short_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(c_1short,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_1short_OLSvsIV_50 <- make_TeXtables(c_1short_OLSvsIV_50, "\\pi^2", "\\rho")

c_1short_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(c_1short,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
c_1short_OLSvsIV_100 <- make_TeXtables(c_1short_OLSvsIV_100, "\\pi^2", "\\rho")

c_1short_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(c_1short,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_1short_OLSvsIV_500 <- make_TeXtables(c_1short_OLSvsIV_500, "\\pi^2", "\\rho")

#--------------------- Relative Width of 1-Step Shortest CI
w_1short_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(w_1short / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
w_1short_OLSvsIV_50 <- make_TeXtables(w_1short_OLSvsIV_50, "\\pi^2", "\\rho")

w_1short_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(w_1short / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
w_1short_OLSvsIV_100 <- make_TeXtables(w_1short_OLSvsIV_100, "\\pi^2", "\\rho")

w_1short_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(w_1short / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_500)
w_1short_OLSvsIV_500 <- make_TeXtables(w_1short_OLSvsIV_500, "\\pi^2", "\\rho")



#--------------------- Coverage of 2-Step CI a1 = a2
c_2equal_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(c_2equal,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_2equal_OLSvsIV_50 <- make_TeXtables(c_2equal_OLSvsIV_50, "\\pi^2", "\\rho")

c_2equal_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(c_2equal,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
c_2equal_OLSvsIV_100 <- make_TeXtables(c_2equal_OLSvsIV_100, "\\pi^2", "\\rho")

c_2equal_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(c_2equal,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_2equal_OLSvsIV_500 <- make_TeXtables(c_2equal_OLSvsIV_500, "\\pi^2", "\\rho")

#--------------------- Relative Width of 2-Step CI a1 = a2
w_2equal_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(w_2equal / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
w_2equal_OLSvsIV_50 <- make_TeXtables(w_2equal_OLSvsIV_50, "\\pi^2", "\\rho")

w_2equal_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(w_2equal / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
w_2equal_OLSvsIV_100 <- make_TeXtables(w_2equal_OLSvsIV_100, "\\pi^2", "\\rho")

w_2equal_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(w_2equal / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_500)
w_2equal_OLSvsIV_500 <- make_TeXtables(w_2equal_OLSvsIV_500, "\\pi^2", "\\rho")



#--------------------- Coverage of 2-Step CI a1 < a2
c_2tauwide_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(c_2tauwide,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_2tauwide_OLSvsIV_50 <- make_TeXtables(c_2tauwide_OLSvsIV_50, "\\pi^2", "\\rho")

c_2tauwide_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(c_2tauwide,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
c_2tauwide_OLSvsIV_100 <- make_TeXtables(c_2tauwide_OLSvsIV_100, "\\pi^2", "\\rho")

c_2tauwide_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(c_2tauwide,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_2tauwide_OLSvsIV_500 <- make_TeXtables(c_2tauwide_OLSvsIV_500, "\\pi^2", "\\rho")

#--------------------- Relative Width of 2-Step CI a1 < a2
w_2tauwide_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(w_2tauwide / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
w_2tauwide_OLSvsIV_50 <- make_TeXtables(w_2tauwide_OLSvsIV_50, "\\pi^2", "\\rho")

w_2tauwide_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(w_2tauwide / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
w_2tauwide_OLSvsIV_100 <- make_TeXtables(w_2tauwide_OLSvsIV_100, "\\pi^2", "\\rho")

w_2tauwide_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(w_2tauwide / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_500)
w_2tauwide_OLSvsIV_500 <- make_TeXtables(w_2tauwide_OLSvsIV_500, "\\pi^2", "\\rho")



#--------------------- Coverage of 2-Step CI a1 > a2
c_2taunarrow_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(c_2taunarrow,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_2taunarrow_OLSvsIV_50 <- make_TeXtables(c_2taunarrow_OLSvsIV_50, "\\pi^2", "\\rho")

c_2taunarrow_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(c_2taunarrow,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
c_2taunarrow_OLSvsIV_100 <- make_TeXtables(c_2taunarrow_OLSvsIV_100, "\\pi^2", "\\rho")

c_2taunarrow_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(c_2taunarrow,2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
c_2taunarrow_OLSvsIV_500 <- make_TeXtables(c_2taunarrow_OLSvsIV_500, "\\pi^2", "\\rho")

#--------------------- Relative Width of 2-Step CI a1 > a2
w_2taunarrow_OLSvsIV_50 <- fmscr::rtables(
  I(100 * round(w_2taunarrow / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_50)
w_2taunarrow_OLSvsIV_50 <- make_TeXtables(w_2taunarrow_OLSvsIV_50, "\\pi^2", "\\rho")

w_2taunarrow_OLSvsIV_100 <- fmscr::rtables(
  I(100 * round(w_2taunarrow / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_100)
w_2taunarrow_OLSvsIV_100 <- make_TeXtables(w_2taunarrow_OLSvsIV_100, "\\pi^2", "\\rho")

w_2taunarrow_OLSvsIV_500 <- fmscr::rtables(
  I(100 * round(w_2taunarrow / w_tsls, 2)) ~ pi_sq + rho + alpha, OLSvsIV_CIs_500)
w_2taunarrow_OLSvsIV_500 <- make_TeXtables(w_2taunarrow_OLSvsIV_500, "\\pi^2", "\\rho")



#============================ Clean Up
rm(OLSvsIV_CIs, OLSvsIV_CIs_50, OLSvsIV_CIs_100, OLSvsIV_CIs_500)




#============================ Choosing IVs Example
load("chooseIVs_CIs.Rd")

chooseIVs_CIs_50 <- subset(chooseIVs_CIs, N == 50)
chooseIVs_CIs_100 <- subset(chooseIVs_CIs, N == 100)
chooseIVs_CIs_500 <- subset(chooseIVs_CIs, N == 500)

#--------------------- Coverage of Valid Estimator
c_valid_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(c_valid,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_valid_chooseIVs_50 <- make_TeXtables(c_valid_chooseIVs_50, "\\gamma^2", "\\rho")

c_valid_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(c_valid,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
c_valid_chooseIVs_100 <- make_TeXtables(c_valid_chooseIVs_100, "\\gamma^2", "\\rho")

c_valid_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(c_valid,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_valid_chooseIVs_500 <- make_TeXtables(c_valid_chooseIVs_500, "\\gamma^2", "\\rho")



#--------------------- Coverage of Naive CI
c_naive_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(c_naive,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_naive_chooseIVs_50 <- make_TeXtables(c_naive_chooseIVs_50, "\\gamma^2", "\\rho")

c_naive_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(c_naive,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
c_naive_chooseIVs_100 <- make_TeXtables(c_naive_chooseIVs_100, "\\gamma^2", "\\rho")

c_naive_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(c_naive,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_naive_chooseIVs_500 <- make_TeXtables(c_naive_chooseIVs_500, "\\gamma^2", "\\rho")

#--------------------- Relative Width of Naive CI
w_naive_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(w_naive / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
w_naive_chooseIVs_50 <- make_TeXtables(w_naive_chooseIVs_50, "\\gamma^2", "\\rho")

w_naive_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(w_naive / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
w_naive_chooseIVs_100 <- make_TeXtables(w_naive_chooseIVs_100, "\\gamma^2", "\\rho")

w_naive_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(w_naive / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_500)
w_naive_chooseIVs_500 <- make_TeXtables(w_naive_chooseIVs_500, "\\gamma^2", "\\rho")



#--------------------- Coverage of 1-Step Equal-Tailed CI
c_1equal_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(c_1equal,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_1equal_chooseIVs_50 <- make_TeXtables(c_1equal_chooseIVs_50, "\\gamma^2", "\\rho")

c_1equal_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(c_1equal,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
c_1equal_chooseIVs_100 <- make_TeXtables(c_1equal_chooseIVs_100, "\\gamma^2", "\\rho")

c_1equal_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(c_1equal,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_1equal_chooseIVs_500 <- make_TeXtables(c_1equal_chooseIVs_500, "\\gamma^2", "\\rho")

#--------------------- Relative Width of 1-Step Equal-Tailed CI
w_1equal_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(w_1equal / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
w_1equal_chooseIVs_50 <- make_TeXtables(w_1equal_chooseIVs_50, "\\gamma^2", "\\rho")

w_1equal_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(w_1equal / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
w_1equal_chooseIVs_100 <- make_TeXtables(w_1equal_chooseIVs_100, "\\gamma^2", "\\rho")

w_1equal_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(w_1equal / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_500)
w_1equal_chooseIVs_500 <- make_TeXtables(w_1equal_chooseIVs_500, "\\gamma^2", "\\rho")



#--------------------- Coverage of 1-Step Shortest CI
c_1short_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(c_1short,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_1short_chooseIVs_50 <- make_TeXtables(c_1short_chooseIVs_50, "\\gamma^2", "\\rho")

c_1short_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(c_1short,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
c_1short_chooseIVs_100 <- make_TeXtables(c_1short_chooseIVs_100, "\\gamma^2", "\\rho")

c_1short_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(c_1short,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_1short_chooseIVs_500 <- make_TeXtables(c_1short_chooseIVs_500, "\\gamma^2", "\\rho")

#--------------------- Relative Width of 1-Step Shortest CI
w_1short_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(w_1short / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
w_1short_chooseIVs_50 <- make_TeXtables(w_1short_chooseIVs_50, "\\gamma^2", "\\rho")

w_1short_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(w_1short / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
w_1short_chooseIVs_100 <- make_TeXtables(w_1short_chooseIVs_100, "\\gamma^2", "\\rho")

w_1short_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(w_1short / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_500)
w_1short_chooseIVs_500 <- make_TeXtables(w_1short_chooseIVs_500, "\\gamma^2", "\\rho")



#--------------------- Coverage of 2-Step CI a1 = a2
c_2equal_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(c_2equal,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_2equal_chooseIVs_50 <- make_TeXtables(c_2equal_chooseIVs_50, "\\gamma^2", "\\rho")

c_2equal_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(c_2equal,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
c_2equal_chooseIVs_100 <- make_TeXtables(c_2equal_chooseIVs_100, "\\gamma^2", "\\rho")

c_2equal_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(c_2equal,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_2equal_chooseIVs_500 <- make_TeXtables(c_2equal_chooseIVs_500, "\\gamma^2", "\\rho")

#--------------------- Relative Width of 2-Step CI a1 = a2
w_2equal_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(w_2equal / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
w_2equal_chooseIVs_50 <- make_TeXtables(w_2equal_chooseIVs_50, "\\gamma^2", "\\rho")

w_2equal_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(w_2equal / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
w_2equal_chooseIVs_100 <- make_TeXtables(w_2equal_chooseIVs_100, "\\gamma^2", "\\rho")

w_2equal_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(w_2equal / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_500)
w_2equal_chooseIVs_500 <- make_TeXtables(w_2equal_chooseIVs_500, "\\gamma^2", "\\rho")



#--------------------- Coverage of 2-Step CI a1 < a2
c_2tauwide_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(c_2tauwide,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_2tauwide_chooseIVs_50 <- make_TeXtables(c_2tauwide_chooseIVs_50, "\\gamma^2", "\\rho")

c_2tauwide_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(c_2tauwide,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
c_2tauwide_chooseIVs_100 <- make_TeXtables(c_2tauwide_chooseIVs_100, "\\gamma^2", "\\rho")

c_2tauwide_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(c_2tauwide,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_2tauwide_chooseIVs_500 <- make_TeXtables(c_2tauwide_chooseIVs_500, "\\gamma^2", "\\rho")

#--------------------- Relative Width of 2-Step CI a1 < a2
w_2tauwide_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(w_2tauwide / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
w_2tauwide_chooseIVs_50 <- make_TeXtables(w_2tauwide_chooseIVs_50, "\\gamma^2", "\\rho")

w_2tauwide_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(w_2tauwide / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
w_2tauwide_chooseIVs_100 <- make_TeXtables(w_2tauwide_chooseIVs_100, "\\gamma^2", "\\rho")

w_2tauwide_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(w_2tauwide / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_500)
w_2tauwide_chooseIVs_500 <- make_TeXtables(w_2tauwide_chooseIVs_500, "\\gamma^2", "\\rho")



#--------------------- Coverage of 2-Step CI a1 > a2
c_2taunarrow_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(c_2taunarrow,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_2taunarrow_chooseIVs_50 <- make_TeXtables(c_2taunarrow_chooseIVs_50, "\\gamma^2", "\\rho")

c_2taunarrow_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(c_2taunarrow,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
c_2taunarrow_chooseIVs_100 <- make_TeXtables(c_2taunarrow_chooseIVs_100, "\\gamma^2", "\\rho")

c_2taunarrow_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(c_2taunarrow,2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
c_2taunarrow_chooseIVs_500 <- make_TeXtables(c_2taunarrow_chooseIVs_500, "\\gamma^2", "\\rho")

#--------------------- Relative Width of 2-Step CI a1 > a2
w_2taunarrow_chooseIVs_50 <- fmscr::rtables(
  I(100 * round(w_2taunarrow / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_50)
w_2taunarrow_chooseIVs_50 <- make_TeXtables(w_2taunarrow_chooseIVs_50, "\\gamma^2", "\\rho")

w_2taunarrow_chooseIVs_100 <- fmscr::rtables(
  I(100 * round(w_2taunarrow / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_100)
w_2taunarrow_chooseIVs_100 <- make_TeXtables(w_2taunarrow_chooseIVs_100, "\\gamma^2", "\\rho")

w_2taunarrow_chooseIVs_500 <- fmscr::rtables(
  I(100 * round(w_2taunarrow / w_valid, 2)) ~ g_sq + rho + alpha, chooseIVs_CIs_500)
w_2taunarrow_chooseIVs_500 <- make_TeXtables(w_2taunarrow_chooseIVs_500, "\\gamma^2", "\\rho")

#================================= Clean Up
rm(chooseIVs_CIs_500, chooseIVs_CIs_100, chooseIVs_CIs_50, chooseIVs_CIs)
rm(make_TeXtables)

#=========================== Output Tables as TeX
setwd("~/fmsc/AdditionalSimulations/CISimResults/")
tables_list <- ls()

for(i in 1:length(tables_list)){
  temp_file <- paste0(tables_list[i], ".tex")
  out_table <- get(tables_list[i])
  cat(out_table, file = paste0("./", temp_file))
}  

rm(list = ls())