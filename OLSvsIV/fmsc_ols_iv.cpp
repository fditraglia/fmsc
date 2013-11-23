/*------------------------------------------------------------
Filename:        fmsc_ols_iv.cpp
Author:          Frank DiTraglia
First Version:   2013-23-11
This Version:    2013-23-11
--------------------------------------------------------------
Carries out FMSC calculations for IV versus 2SLS example with
a single endogenous regressor.
------------------------------------------------------------*/
  
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List fmsc_ols_iv_cpp(x, y, z, DHW_levels){
  /*------------------------------------------------------------------
  NOTE: This function assumes that x is a column of observations for 
        a single endogenous regressor. In other words, it assumes 
        that all exogenous regressors, including the constant, have 
        been "projected out" of the system.
  ------------------------------------------------------------------
  Arguments:
   x               matrix of observations of single endog. regressor
   y               matrix of observations for outcome variable
   z               matrix of observations for the instruments
   DHW.levels      optional vector of significance levels for a 
                     Durbin-Hausman-Wu tests. If specified, the 
                     function returns, in addition to the OLS, IV 
                     and FMSC-selected estimators, the corresponding 
                     DHW pretest estimators.
  ------------------------------------------------------------------*/
  return;
  
}