/*------------------------------------------------------------
Filename:        simulation_functions_OLSvsIV.cpp
Author:          Frank DiTraglia
First Version:   2014-01-22
This Version:    2014-01-22
--------------------------------------------------------------
C++ code to carry out the OLS vs. IV simulation experiments.
All functions here correspond to R versions contained in the 
file simulation_functions_OLSvsIV.R.
------------------------------------------------------------*/
  
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List dgp_cpp(double b, NumericVector PI_r, NumericMatrix Ve_r, 
              NumericMatrix Vz_r, int n){
/*------------------------------------------------------------------
# Arguments:
#  b       coefficient on x in the second stage
#  PI      vector of coefficients on z in the second stage
#  V.e     variance-covariance matrix of the errors epsilon and v
#  V.z     variance-covariance matrix of the instruments z
#  n       sample size
#------------------------------------------------------------------- 
# Returns: list of matrices x, y and z containing simulated dataset
-------------------------------------------------------------------*/ 
    
    RNGScope scope; // also done by sourceCpp()
    
    //Number of instruments
    int n_z = Vz_r.ncol();
    
    //Initialize armadillo matrices corresponding to R input matrices
    //that reuse the original memory
    arma::mat Ve(Ve_r.begin(), 2, 2, false); //Error terms 
    arma::mat Vz(Vz_r.begin(), n_z, n_z, false); //Instruments
    arma::colvec p(PI_r.begin(), PI_r.size(), false); //Can't call it PI!
    
    //Generate errors
    arma::colvec stdnorm_errors = rnorm(n * 2);
    arma::mat errors = trans(arma::chol(Ve) * reshape(stdnorm_errors, 2, n));
    
    //Generate instruments
    arma::colvec stdnorm_z = rnorm(n * n_z);
    arma::mat z = trans(arma::chol(Vz) * reshape(stdnorm_z, n_z, n));

    //Generate endogenous variables: x and y
    arma::colvec x = z * p + errors.col(1); //Remember: zero indexing!
    arma::colvec y = b * x + errors.col(0); 
    
    return List::create(Named("x") = x, Named("y") = y, Named("z") = z);
                        
}

// [[Rcpp::export]]
NumericVector fmsc_ols_iv_cpp(NumericVector x_r, NumericVector y_r, 
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
  
  //Calculate Optimal Weight for OLS estimator
  //Plug in asymptotically unbiased estimator of squared A-Bias
  //If A-bias estimate is negative, set to zero. Make sure weight
  //lies in [0,1].
  double tau_squared_est = pow(tau, 2) - (s_e_squared * s_x_squared 
                                        * s_v_squared / g_squared);
  double bias_est;
  if((tau_squared_est / pow(s_x_squared, 2)) >= 0){
    bias_est = tau_squared_est / pow(s_x_squared, 2);
  } else {
    bias_est = 0;
  }
  double var_diff = s_e_squared * (1 / g_squared - 1 / s_x_squared);
  double omega_star = 1 / (1 + (bias_est / var_diff));
  if(omega_star > 1){
    omega_star = 1;
  }
  if(omega_star < 0){
    omega_star = 0;
  }
  
  //Optimal Averaging Estimator
  double b_star = omega_star * b_ols + (1 - omega_star) * b_tsls;
  
  //Create and return vector of results
  NumericVector out = NumericVector::create(b_ols, b_tsls, b_fmsc, b_star);
  out.names() = CharacterVector::create("b.ols", "b.tsls", "b.fmsc", "b.star");
  return out;  
  
}

