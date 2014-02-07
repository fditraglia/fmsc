/*------------------------------------------------------------
Filename:        simulation_functions_OLSvsIV_classes.cpp
Author:          Frank DiTraglia
First Version:   2014-02-07
This Version:    2014-02-07
--------------------------------------------------------------
C++ code to carry out the OLS vs. IV simulation experiments.
Unlike the versions in simulation_functions_OLSvsIV.cpp, the
code here uses classes.
------------------------------------------------------------*/
  
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;




class dgp_OLS_IV {
  private:
    int n_z; 
    arma::colvec stdnorm_errors, stdnorm_z; 
    arma::mat errors;
  public:
    arma::colvec x; 
    arma::colvec y; 
    arma::mat z;
    //Class constructor
    dgp_OLS_IV(double, arma::colvec, arma::mat, arma::mat, int);
};


dgp_OLS_IV::dgp_OLS_IV(double b, arma::colvec p, arma::mat Ve, 
              arma::mat Vz, int n){
//Constructor for class dgp_OLS_IV
//Generates simulations from the dgp

  RNGScope scope;
  n_z = Vz.n_cols;
      
  //Generate errors
  stdnorm_errors = rnorm(n * 2);
  errors = trans(arma::chol(Ve) * reshape(stdnorm_errors, 2, n));
      
  //Generate instruments
  stdnorm_z = rnorm(n * n_z);
  z = trans(arma::chol(Vz) * reshape(stdnorm_z, n_z, n));

  //Generate endogenous variables: x and y
  x = z * p + errors.col(1); //Remember: zero indexing!
  y = b * x + errors.col(0);   
}






class fmsc_OLS_IV {
  private:
    int n_z, n;
    double xx, g_squared, s_e_squared, s_x_squared, 
        s_v_squared, tau, ols_estimate, tsls_estimate;
    arma::colvec tsls_resid, first_stage_coefs, zx;
    arma::mat Qz, Rz, Rzinv, zz_inv;
  public:
    //class constructor
    fmsc_OLS_IV(arma::colvec, arma::colvec, arma::mat);
    //member functions
    double Tfmsc();
    double b_ols();
    double b_tsls();
    double b_fmsc();
    double b_DHW(double);
    double b_AVG();
};
  

fmsc_OLS_IV::fmsc_OLS_IV(arma::colvec x, arma::colvec y, 
                          arma::mat z){
//Class constructor for fmsc_OLS_IV
//Calculates all "basic" quantities 

  n_z = z.n_cols;
  n = z.n_rows;
      
  //OLS estimator
  xx = arma::dot(x, x);
  ols_estimate = arma::dot(x, y) / xx;
  
  //2SLS Estimator
  first_stage_coefs = arma::solve(z,x);
  tsls_estimate = arma::as_scalar(arma::solve(z * first_stage_coefs, y)); 
  tsls_resid = y - x * tsls_estimate; 
  
  //2SLS Weighting Matrix and "gamma-squared"
  arma::qr_econ(Qz, Rz, z);
  Rzinv = arma::inv(arma::trimatu(Rz));
  zz_inv = Rzinv * arma::trans(Rzinv);
  zx = arma::trans(z) * x;
  g_squared = arma::as_scalar(arma::trans(zx) * zz_inv * zx) / n;

  //Remaining quantities needed for FMSC
  s_e_squared = arma::dot(tsls_resid, tsls_resid) / n;
  s_x_squared = xx / n;
  s_v_squared = s_x_squared - g_squared;
  tau = arma::dot(x, tsls_resid) / sqrt(n);
}

double fmsc_OLS_IV::b_ols(){
//Member function of class fmsc_OLS_IV
//Extracts OLS estimate
  return(ols_estimate);
}


double fmsc_OLS_IV::b_tsls(){
//Member function of class fmsc_OLS_IV
//Extracts 2SLS estimate
  return(tsls_estimate);
}
    
double fmsc_OLS_IV::Tfmsc(){
//Member function of class fmsc_OLS_IV
//Calculates FMSC "test statistic"
  return(pow(tau, 2) * g_squared / (s_v_squared * 
                                  s_e_squared * s_x_squared));      
}
    
double fmsc_OLS_IV::b_fmsc(){
//Member function of class fmsc_OLS_IV
//Calculates estimator selected by FMS
  double out;
  if(Tfmsc() < 2){
    out = ols_estimate;
  } else {
    out = tsls_estimate;
  }
  return(out);
}
    

double fmsc_OLS_IV::b_DHW(double level){
//Member function of class fmsc_OLS_IV
//Calculates DHW pre-test estimator at a given confidence level
  double out;
  double DHW_crit = R::qchisq(level, 1, 1, 0);
  if(Tfmsc() > DHW_crit){
    out = tsls_estimate;
  }else{
    out = ols_estimate;
  }
  return(out);
}


double fmsc_OLS_IV::b_AVG(){
//Member function of class fmsc_OLS_IV
//Calculates feasible version of AMSE-optimal averaging estimator
  double tau_squared_est = pow(tau, 2) - (s_e_squared * s_x_squared 
                                        * s_v_squared / g_squared);
  double sq_bias_est;
  //If the squared bias estimate is negative, set it to zero
  //so weight lies in [0,1]
  if((tau_squared_est / pow(s_x_squared, 2)) >= 0){
    sq_bias_est = tau_squared_est / pow(s_x_squared, 2);
  } else {
    sq_bias_est = 0;
  }
  double var_diff = s_e_squared * (1 / g_squared - 1 / s_x_squared);
  double omega = 1 / (1 + (sq_bias_est / var_diff));  
  double out = omega * ols_estimate + (1 - omega) * tsls_estimate;
  return(out);
}





double MSE_trim(arma::colvec x, double truth, double trim){
/*-------------------------------------------------------
# Calculates trimmed mean-squared error.
#--------------------------------------------------------
#  x        vector of estimates
#  truth    true value of the parameter
#  trim     fraction of estimates to discard (half from
#             each tail) before calculating MSE
#-------------------------------------------------------*/

  int k = x.n_elem;
  int tail_drop = ceil(k * trim / 2);
  
  arma::colvec x_trimmed = arma::sort(x);
  x_trimmed = x_trimmed(arma::span(tail_drop, k - tail_drop - 1));
  
  arma::colvec truth_vec = truth * arma::ones(x_trimmed.n_elem);
  arma::colvec errors = x_trimmed - truth_vec;
  double MSE = arma::dot(errors, errors) / errors.n_elem;
  
  return(MSE);
  
}



// [[Rcpp::export]]
NumericVector mse_compare_classes(double b, arma::colvec p, arma::mat Ve, 
              arma::mat Vz, int n, int n_reps){
//Function to run n_reps of the simulation study and calculate the MSE
//of the various estimators
  
  arma::colvec ols(n_reps);
  arma::colvec tsls(n_reps);
  arma::colvec fmsc(n_reps);
  arma::colvec DHW90(n_reps);
  arma::colvec DHW95(n_reps);
  arma::colvec AVG(n_reps);

  for(int i = 0; i < n_reps; i++){
    
    dgp_OLS_IV sim(b, p, Ve, Vz, n);
    fmsc_OLS_IV est(sim.x, sim.y, sim.z);

    ols(i) = est.b_ols();
    tsls(i) = est.b_tsls();
    fmsc(i) = est.b_fmsc();
    DHW90(i) = est.b_DHW(0.90);
    DHW95(i) = est.b_DHW(0.95);
    AVG(i) = est.b_AVG();
    
  }
  
  double MSE_ols = MSE_trim(ols, b, 0);
  double MSE_tsls = MSE_trim(tsls, b, 0);
  double MSE_fmsc = MSE_trim(fmsc, b, 0);
  double MSE_DHW90 = MSE_trim(DHW90, b, 0);
  double MSE_DHW95 = MSE_trim(DHW95, b, 0);
  double MSE_star = MSE_trim(AVG, b, 0);
  
  //Create and return vector of results
  NumericVector out = NumericVector::create(MSE_ols, MSE_tsls, MSE_fmsc, 
                        MSE_star, MSE_DHW90, MSE_DHW95);
  out.names() = CharacterVector::create("b.ols", "b.tsls", "b.fmsc", 
                        "b.star", "b.DHW90", "b.DHW95");
  return(out);  

  
}







// [[Rcpp::export]]
NumericVector mse_compare_default_classes(double p , double r, int n, 
                                        int n_reps){
//Runs the simulation once with "default" values
//for "uninteresting" parameters.

  double b = 1;
  arma::colvec p_vec = p * arma::ones(3);
  arma::mat Vz = arma::eye(3, 3);
  
  arma::mat Ve;
  Ve << 1<< r << arma::endr
     << r << 1 << arma::endr;
  
  NumericVector out = mse_compare_classes(b, p_vec, Ve, Vz, n, n_reps);
  return(out);
  
}
