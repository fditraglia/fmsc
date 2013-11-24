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
  ------------------------------------------------------------------*/
  int n = y_r.size();
  int n_z = z_r.ncol();
  
  //Initialize armadillo matrices and vectors corresponding to 
  //R input and reuse original memory
  arma::colvec x(x_r.begin(), n, false);
  arma::colvec y(y_r.begin(), n, false);
  arma::mat z(z_r.begin(), n, n_z, false);
  
  //Calculations for OLS estimator
  double xx = arma::dot(x, x);
  double b_ols = arma::dot(x, y) / xx;
  
  //Calculations for 2SLS Estimator
  arma::colvec first_stage_coefs = arma::solve(z,x);
  double b_tsls = arma::as_scalar(arma::solve(z * first_stage_coefs, y)); 
  arma::colvec tsls_resid = y - x * b_tsls; 
  
  //2SLS Weighting Matrix and "gamma-squared"
  arma::mat Qz, Rz;
  arma::qr_econ(Qz, Rz, z);
  arma::mat Rzinv = arma::inv(arma::trimatu(Rz));
  arma::mat zz_inv = Rzinv * arma::trans(Rzinv);
  arma::colvec zx = arma::trans(z) * x;
  double g_squared = arma::as_scalar(arma::trans(zx) * zz_inv * zx) / n;

  //Remaining quantities needed for FMSC
  double s_e_squared = arma::dot(tsls_resid, tsls_resid) / n;
  double s_x_squared = xx / n;
  double s_v_squared = s_x_squared - g_squared;
  double tau = arma::dot(x, tsls_resid) / sqrt(n);
  
  //Calculate FMSC "test statistic" and selected estimator
  double Tfmsc = pow(tau, 2) * g_squared / (s_v_squared * 
                                      s_e_squared * s_x_squared);
  double b_fmsc;
  if(Tfmsc < 2){
    b_fmsc = b_ols;
  } else {
    b_fmsc = b_tsls;
  }
  
  return List::create(Named("b.ols") = b_ols, Named("b.tsls") = b_tsls, 
                                      Named("b.fmsc") = b_fmsc);
  
}