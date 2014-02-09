/*------------------------------------------------------------
# Filename:        simulation_functions_OLSvsIV.cpp
# Author:          Frank DiTraglia
# First Version:   2014-01-22
# This Version:    2014-02-08
#--------------------------------------------------------------
# C++ code to carry out the OLS vs. IV simulation experiments.
#------------------------------------------------------------*/
  
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;




class dgp_OLS_IV {
/*--------------------------------------------------
# Class for simulating from dgp in OLS vs IV example
#-------------------------------------------------*/
  public:
    arma::colvec x;  //sims for endogenous regressor
    arma::colvec y;  //sims for outcome
    arma::mat z;     //sims for instrumental vars
    //Class constructor
    dgp_OLS_IV(double, arma::colvec, arma::mat, arma::mat, int);
};

dgp_OLS_IV::dgp_OLS_IV(double b, arma::colvec p, arma::mat Ve, 
              arma::mat Vz, int n){
/*--------------------------------------------------
# Constructor: generates simulations from dgp
#---------------------------------------------------
# Arguments:
#  b       coefficient on x in the 2nd stage
#  PI      vector of coeffs on z in the 1st stage
#  V.e     var-covar matrix of epsilon and v 
#  V.z     var-covar matrix z
#  n       sample size
#---------------------------------------------------
# Initializes:
#     x, y, z
#-------------------------------------------------*/
  RNGScope scope;
  int n_z = Vz.n_cols;
  
  arma::colvec stdnorm_errors = rnorm(n * 2);
  arma::mat errors = arma::trans(arma::chol(Ve) * reshape(stdnorm_errors, 2, n));
  arma::colvec stdnorm_z = rnorm(n * n_z);
  
  z = trans(arma::chol(Vz) * reshape(stdnorm_z, n_z, n));
  x = z * p + errors.col(1); //Remember: zero indexing!
  y = b * x + errors.col(0);   
}






class fmsc_OLS_IV {
/*--------------------------------------------------
# Class for FMSC calculations in OLS vs IV example
#-------------------------------------------------*/
  private:
    int n_z, n;
    double xx, g_sq, s_e_sq_ols, s_e_sq_tsls, s_x_sq, 
        s_v_sq, tau, tau_var, ols_estimate, tsls_estimate;
    arma::colvec ols_resid, tsls_resid, first_stage_coefs, zx;
    arma::mat zz_inv, zz, CI_sims;
    void draw_CI_sims(int); //Initialize CI_sims for later use
  public:
    //class constructor
    fmsc_OLS_IV(arma::colvec, arma::colvec, arma::mat);
    //member functions
    double Tfmsc();       //return FMSC "test statistic"
    double b_ols();       //return OLS estimate
    double b_tsls();      //return tsls estimate
    double b_fmsc();      //return FMSC-selected estimate
    double b_DHW(double); //return DHW pre-test estimate
    double b_AVG();       //return feasible averaging estimate
    arma::rowvec CI_ols(double);  //Confidence interval (CI) for ols
    arma::rowvec CI_tsls(double); //CI for tsls
    arma::rowvec CI_fmsc_naive(double); //Naive CI post-fmsc
    arma::rowvec CI_tau(double); //CI for bias parameter tau
//    arma::rowvec CI_fmsc_correct(double, int); //Corrected CI post-fmsc
};
  

fmsc_OLS_IV::fmsc_OLS_IV(arma::colvec x, arma::colvec y, 
                          arma::mat z){
/*--------------------------------------------------
# Constructor: calculates all "basic quantities"
#---------------------------------------------------
# Arguments:
#  x          vector of obs. for endog regressor
#  y          vector of obs. for outcome
#  z          matrix of obs. for instruments
#---------------------------------------------------
# Initializes:
#    n_z, n, xx, g_sq, s_e_sq_ols, s_e_tsls, 
#    s_x_sq, s_v_sq, tau, tau_var, ols_estimate, 
#    ols_resid, tsls_estimate, tsls_resid, 
#    first_stage_coefs, zx, zz_inv, zz
#-------------------------------------------------*/
  n_z = z.n_cols;
  n = z.n_rows;

  xx = arma::dot(x, x);
  ols_estimate = arma::dot(x, y) / xx;
  ols_resid = y - x * ols_estimate;
  
  first_stage_coefs = arma::solve(z,x);
  tsls_estimate = arma::as_scalar(arma::solve(z * first_stage_coefs, y)); 
  tsls_resid = y - x * tsls_estimate; 
  
  arma::mat Qz, Rz, Rzinv;
  arma::qr_econ(Qz, Rz, z);
  Rzinv = arma::inv(arma::trimatu(Rz));
  zz_inv = Rzinv * arma::trans(Rzinv);
  zz = trans(trimatu(Rz)) * trimatu(Rz);
  zx = arma::trans(z) * x;
  g_sq = arma::as_scalar(arma::trans(zx) * zz_inv * zx) / n;
  
  s_e_sq_ols = arma::dot(ols_resid, ols_resid) / n;
  s_e_sq_tsls = arma::dot(tsls_resid, tsls_resid) / n;
  s_x_sq = xx / n;
  s_v_sq = s_x_sq - g_sq;
  tau = arma::dot(x, tsls_resid) / sqrt(n);
  tau_var = s_e_sq_tsls * s_x_sq * (s_x_sq / g_sq - 1);
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
  return(pow(tau, 2) * g_sq / (s_v_sq * 
                                  s_e_sq_tsls * s_x_sq));      
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
    
double fmsc_OLS_IV::b_DHW(double alpha){
//Member function of class fmsc_OLS_IV
//Calculates DHW pre-test estimator
//Arguments: alpha = significance level for test
  double out;
  double DHW_crit = R::qchisq(1 - alpha, 1, 1, 0);
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
  double tau_squared_est = pow(tau, 2) - (s_e_sq_tsls * s_x_sq 
                                        * s_v_sq / g_sq);
  double sq_bias_est;
  //If the squared bias estimate is negative, set it to zero
  //so weight lies in [0,1]
  if((tau_squared_est / pow(s_x_sq, 2)) >= 0){
    sq_bias_est = tau_squared_est / pow(s_x_sq, 2);
  } else {
    sq_bias_est = 0;
  }
  double var_diff = s_e_sq_tsls * (1 / g_sq - 1 / s_x_sq);
  double omega = 1 / (1 + (sq_bias_est / var_diff));  
  double out = omega * ols_estimate + (1 - omega) * tsls_estimate;
  return(out);
}

arma::rowvec fmsc_OLS_IV::CI_ols(double alpha){
//Member function of class fmsc_OLS_IV
//Returns confidence interval (lower, upper) for OLS estimator
//Arguments: alpha = significance level (0.05 = 95% CI)
  double z_quantile = R::qnorm(1 - alpha/2, 0, 1, 1, 0);
  double SE_ols = sqrt(s_e_sq_ols / (n * s_x_sq));
  double lower = ols_estimate - z_quantile * SE_ols;
  double upper = ols_estimate + z_quantile * SE_ols;
  
  arma::rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}

arma::rowvec fmsc_OLS_IV::CI_tsls(double alpha){
//Member function of class fmsc_OLS_IV
//Returns confidence interval (lower, upper) for TSLS estimator
//Arguments: alpha = significance level (0.05 = 95% CI)
  double z_quantile = R::qnorm(1 - alpha/2, 0, 1, 1, 0);
  double SE_tsls = sqrt(s_e_sq_tsls / (n * g_sq));
  double lower = tsls_estimate - z_quantile * SE_tsls;
  double upper = tsls_estimate + z_quantile * SE_tsls;
  
  arma::rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}

arma::rowvec fmsc_OLS_IV::CI_fmsc_naive(double alpha){
//Member function of class fmsc_OLS_IV
//Returns naive confidence interval (lower, upper)
//for post-FMSC estimator
//Arguments: alpha = significance level (0.05 = 95% CI)
    arma::rowvec out(2);
    if(Tfmsc() < 2){
      out = CI_ols(alpha);
    } else {
      out = CI_tsls(alpha);
    }
  return(out);
}

arma::rowvec fmsc_OLS_IV::CI_tau(double delta){
//Member function of class fmsc_OLS_IV
//Returns CI (lower, upper) for bias parameter tau
//Arguments: delta = significance level (0.05 = 95% CI)
  double z_quantile = R::qnorm(1 - delta/2, 0, 1, 1, 0);
  double SE_tau = sqrt(tau_var);
  double lower = tau - z_quantile * SE_tau;
  double upper = tau + z_quantile * SE_tau;
  
  arma::rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}


void fmsc_OLS_IV::draw_CI_sims(int n_sims){
//Member function of class fmsc_OLS_IV
//Initializes the matrix CI_sims for use by other member functions
  RNGScope scope;
  
  //Construct variance matrix of M
  arma::mat Omega(n_z + 1, n_z + 1);
  Omega(0, 0) = s_x_sq;
  Omega(0, arma::span(1, n_z)) = arma::trans(zx) / n;
  Omega(arma::span(0, n_z), 0) = zx / n;
  Omega(arma::span(1, n_z), arma::span(1, n_z)) = zz / n;
  Omega = s_e_sq_tsls * Omega;
  
  //Construct matrix D that multiplies M
  arma::mat D(3, n_z + 1);
  D(0, 0) = 1 / s_x_sq;
  D(0, arma::span(1, n_z)) = arma::zeros<arma::rowvec>(n_z);
  D(1, 0) = 0;
  D(1, arma::span(1, n_z)) = arma::trans(first_stage_coefs) / g_sq;
  D(2, 0) = 1;
  D(2, arma::span(1, n_z)) = -1 * s_x_sq * arma::trans(first_stage_coefs) / g_sq;
  
  //Draw simulations
  arma::colvec stdnorm = rnorm(n_sims * Omega.n_rows);
  CI_sims = D * arma::chol(Omega) * reshape(stdnorm, Omega.n_rows, n);
}


//arma::rowvec fmsc_OLS_IV::CI_fmsc_correct(double level, int n_sims){
////Member function of class fmsc_OLS_IV
////Returns corrected confidence interval (lower, upper)
////for post-FMSC estimator
////Arguments: level = confidence level (0.95 is a 95% CI)
//}





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

double MAD(arma::colvec x, double truth){
/*-------------------------------------------------------
# Calculates median absolute deviation.
#--------------------------------------------------------
#  x        vector of estimates
#  truth    true value of the parameter
#-------------------------------------------------------*/
  arma::colvec truth_vec = truth * arma::ones(x.n_rows);
  arma::colvec abs_dev = abs(x - truth_vec);
  return(arma::median(abs_dev)); 
}


double coverage_prob(arma::mat conf_intervals, double truth){
/*-------------------------------------------------------
# Calculates the coverage probability of a matrix of
# confidence intervals.
#--------------------------------------------------------
#  conf_intervals   matrix of confidence intervals in 
#                     which each row is a CI, the 1st
#                     column is the lower limit, and the
#                     2nd column is the upper limit 
#                           
#  truth            true value of the parameter for which
#                       the CIs were constructed
#-------------------------------------------------------*/
  arma::colvec truth_vec = truth * arma::ones(conf_intervals.n_rows);
  arma::colvec cover_lower = arma::conv_to<arma::colvec>
                    ::from(conf_intervals.col(0) < truth_vec);
  arma::colvec cover_upper = arma::conv_to<arma::colvec>
                    ::from(conf_intervals.col(1) > truth_vec);
  arma::colvec cover = cover_lower % cover_upper;
  return(arma::sum(cover) / cover.n_elem);
}


double median_width(arma::mat conf_intervals){
/*-------------------------------------------------------
# Calculates the median width of a matrix of confidence
# intervals.
#--------------------------------------------------------
#  conf_intervals   matrix of confidence intervals in 
#                     which each row is a CI, the 1st
#                     column is the lower limit, and the
#                     2nd column is the upper limit 
#-------------------------------------------------------*/
  arma::colvec width = conf_intervals.col(1) - conf_intervals.col(0);
  return(arma::median(width));
}



// [[Rcpp::export]]
NumericVector mse_compare_cpp(double b, arma::colvec p, arma::mat Ve, 
              arma::mat Vz, int n, int n_reps){
//Function to run n_reps of the simulation study and calculate the MSE
//of various estimators

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
    DHW90(i) = est.b_DHW(0.1);
    DHW95(i) = est.b_DHW(0.05);
    AVG(i) = est.b_AVG();
    
  }
  
  double const trim_frac = 0; //Change this if you want trimmed MSE
  
  double MSE_ols = MSE_trim(ols, b, trim_frac);
  double MSE_tsls = MSE_trim(tsls, b, trim_frac);
  double MSE_fmsc = MSE_trim(fmsc, b, trim_frac);
  double MSE_DHW90 = MSE_trim(DHW90, b, trim_frac);
  double MSE_DHW95 = MSE_trim(DHW95, b, trim_frac);
  double MSE_star = MSE_trim(AVG, b, trim_frac);
  
  //Create and return vector of results
  NumericVector out = NumericVector::create(MSE_ols, MSE_tsls, MSE_fmsc, 
                        MSE_star, MSE_DHW90, MSE_DHW95);
  out.names() = CharacterVector::create("b.ols", "b.tsls", "b.fmsc", 
                        "b.star", "b.DHW90", "b.DHW95");
  return(out);  

  
}




// [[Rcpp::export]]
NumericVector mse_compare_default_cpp(double p , double r, int n, 
                                        int n_reps){
//Wrapper for mse_compare_cpp
//Runs simulation once with default values for "uninteresting" params

  double b = 1;
  arma::colvec p_vec = p * arma::ones(3);
  arma::mat Vz = arma::eye(3, 3);
  
  arma::mat Ve;
  Ve << 1<< r << arma::endr
     << r << 1 << arma::endr;
  
  NumericVector out = mse_compare_cpp(b, p_vec, Ve, Vz, n, n_reps);
  return(out);  
}




// [[Rcpp::export]]
arma::mat test_CIs_cpp(double p , double r, int n, 
                                        int n_reps){
//Function to test the confidence interval code
  double b = 1;
  arma::colvec p_vec = p * arma::ones(3);
  arma::mat Vz = arma::eye(3, 3);
  
  arma::mat Ve;
  Ve << 1<< r << arma::endr
     << r << 1 << arma::endr;

  arma::mat out(n_reps, 2, arma::fill::zeros);
  
  for(int i = 0; i < n_reps; i++){
    dgp_OLS_IV sim(b, p_vec, Ve, Vz, n);
    fmsc_OLS_IV est(sim.x, sim.y, sim.z);
    out.row(i) = est.CI_fmsc_naive(0.05); 
  }
  return(out);  
}
