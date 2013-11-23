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
List fmsc_ols_iv_cpp(NumericVector x_r, NumericVector y_r, 
                      NumericMatrix z_r){ 
                        
                      //Add DHW and averaging stuff later
  
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
//  zz.inv <- chol2inv(qr.R(qr(z))) 
//  x.hat <- z %*% (zz.inv %*% zx)
//  b.tsls <- crossprod(x.hat, y) / crossprod(x.hat)
//  tsls.residuals <- y - x %*% b.tsls
//  s.e.squared <- crossprod(tsls.residuals) / n
//  tau <- crossprod(x, tsls.residuals) / sqrt(n)
//  s.x.squared <- xx / n
//  g.squared <- t(zx) %*% zz.inv %*% zx / n
//  s.v.squared <- s.x.squared - g.squared
//  Tfmsc <- tau^2 * g.squared / (s.v.squared * s.e.squared * s.x.squared)
//  b.fmsc <- ifelse(Tfmsc < 2, b.ols, b.tsls)

  int n = y_r.size();
  int n_z = z_r.ncol();
  
  //Initialize armadillo matrices and vectors corresponding to 
  //R input and reuse original memory
  arma::colvec x(x_r.begin(), n, false);
  arma::colvec y(y_r.begin(), n, false);
  arma::mat z(z_r.begin(), n, n_z, false);
  
  //Calculations for OLS estimator
  double xx = dot(x, x);
  double b_ols = dot(x, y) / xx
  
  //Calculations for 2SLS Estimator
  arma::colvec zx = trans(z) * x;
  
  
  
  return;
  
}