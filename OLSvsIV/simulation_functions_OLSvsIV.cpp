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
arma::mat dgp_cpp(double b, arma::colvec p, arma::mat Ve, 
              arma::mat Vz, int n){
/*------------------------------------------------------------------
# Arguments:
#  b       coefficient on x in the second stage
#  PI      vector of coefficients on z in the second stage
#  V.e     variance-covariance matrix of the errors epsilon and v
#  V.z     variance-covariance matrix of the instruments z
#  n       sample size
#------------------------------------------------------------------- 
# Returns: a matrix of simulated data. The first column is x, the
#          single endogenous regressor, the second column is y, the
#          outcome, and all remaining columns are z, the matrix of
#          instrumental variables.
-------------------------------------------------------------------*/ 
    
    RNGScope scope;
    
    //Number of instruments
    int n_z = Vz.n_cols;
    
    //Generate errors
    arma::colvec stdnorm_errors = rnorm(n * 2);
    arma::mat errors = trans(arma::chol(Ve) * reshape(stdnorm_errors, 2, n));
    
    //Generate instruments
    arma::colvec stdnorm_z = rnorm(n * n_z);
    arma::mat z = trans(arma::chol(Vz) * reshape(stdnorm_z, n_z, n));

    //Generate endogenous variables: x and y
    arma::colvec x = z * p + errors.col(1); //Remember: zero indexing!
    arma::colvec y = b * x + errors.col(0); 
    
    
    arma::mat sim_data = arma::mat(n, 2 + n_z);
    sim_data.col(0) = x;
    sim_data.col(1) = y;
    sim_data.cols(2, n_z + 1) = z;
    
    return (sim_data);
    
}





// [[Rcpp::export]]
NumericVector fmsc_ols_iv_cpp(arma::mat data){ 
                        
/*------------------------------------------------------------------
# Estimates OLS and IV estimators (without constant terms) using
# the supplied dataset, calculates the FMSC criterion to choose
# between them, estimates the AMSE-optimal weights on OLS and IV
# and uses these to construct an optimal averaging estimator.
#-------------------------------------------------------------------
# Arguments:
#   data           a matrix of data for OLS/IV. The first column is
#                  x, the single endogenous regressor, the second 
#                  column is y, the outcome, and all remaining
#                  columns are the instruments. 
------------------------------------------------------------------*/

  int n_z = data.n_cols - 2;
  int n = data.n_rows;
  
  arma::colvec x = data.col(0);
  arma::colvec y = data.col(1);
  arma::mat z = data.cols(2, n_z + 1);

  //OLS estimator
  double xx = arma::dot(x, x);
  double b_ols = arma::dot(x, y) / xx;
  
  //2SLS Estimator
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


//Run the simulation once with "default" values
//for "uninteresting" parameters.

// [[Rcpp::export]]
NumericVector simple_sim_cpp(double p , double r, int n){
  
  //default parameter values
  double b = 1;
  arma::colvec p_vec = p * arma::ones(3);
  arma::mat Vz = arma::eye(3, 3);
  
  arma::mat Ve;
  Ve << 1<< r << arma::endr
     << r << 1 << arma::endr;

  
  arma::mat sim_data = dgp_cpp(b, p_vec, Ve, Vz, n);
  
  NumericVector results = fmsc_ols_iv_cpp(sim_data);
  return(results);

  
}



//Next write mse_compare_cpp as well as an mse helper function