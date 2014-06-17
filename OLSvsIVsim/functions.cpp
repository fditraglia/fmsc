/*------------------------------------------------------------
# Filename:        functions.cpp
# Author:          Frank DiTraglia
# First Version:   2014-01-22
# Last Updated:    2014-06-17
#--------------------------------------------------------------
# C++ code to carry out the OLS vs. IV simulation experiments.
#------------------------------------------------------------*/
  
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


double sample_quantile(colvec x, double p){
/*-------------------------------------------------------
# Calculates a sample quantile
#--------------------------------------------------------
#  x        vector of data
#  p        probability for desired quantile (e.g. 0.25
#             gives the first quartile, 0.5 the median)
#--------------------------------------------------------
# Details:
#           There are many competing definitions of
#           sample quantiles (see Hyndman & Fan, 1996).
#           Here we simply use the R default definition,
#           which corresponds to Definition 7 in Hyndman
#           & Fan. See ?quantile in R for more details.
#-------------------------------------------------------*/
  int n = x.n_elem;
  double m = 1 - p;
  int j = floor(n * p + m);
  double g = n * p + m - j;
  colvec x_sort = sort(x);
  return((1 - g) * x_sort(j - 1) + g * x_sort(j));
}

double MSE_trim(colvec x, double truth, double trim){
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
  
  colvec x_trimmed = sort(x);
  x_trimmed = x_trimmed(span(tail_drop, k - tail_drop - 1));
  
  colvec truth_vec = truth * ones(x_trimmed.n_elem);
  colvec errors = x_trimmed - truth_vec;
  double MSE = dot(errors, errors) / errors.n_elem;
  
  return(MSE);  
}

double MAD(colvec x, double truth){
/*-------------------------------------------------------
# Calculates median absolute deviation.
#--------------------------------------------------------
#  x        vector of estimates
#  truth    true value of the parameter
#-------------------------------------------------------*/
  colvec truth_vec = truth * ones(x.n_rows);
  colvec abs_dev = abs(x - truth_vec);
  return(median(abs_dev)); 
}


double coverage_prob(mat conf_intervals, double truth){
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
  colvec truth_vec = truth * ones(conf_intervals.n_rows);
  colvec cover_lower = conv_to<colvec>
                    ::from(conf_intervals.col(0) < truth_vec);
  colvec cover_upper = conv_to<colvec>
                    ::from(conf_intervals.col(1) > truth_vec);
  colvec cover = cover_lower % cover_upper;
  return(sum(cover) / cover.n_elem);
}


double median_width(mat conf_intervals){
/*-------------------------------------------------------
# Calculates the median width of a matrix of confidence
# intervals.
#--------------------------------------------------------
#  conf_intervals   matrix of confidence intervals in 
#                     which each row is a CI, the 1st
#                     column is the lower limit, and the
#                     2nd column is the upper limit 
#-------------------------------------------------------*/
  colvec width = conf_intervals.col(1) - conf_intervals.col(0);
  return(median(width));
}




class dgp_OLS_IV {
/*--------------------------------------------------
# Class for simulating from dgp in OLS vs IV example
#-------------------------------------------------*/
  public:
    colvec x;  //sims for endogenous regressor
    colvec y;  //sims for outcome
    mat z;     //sims for instrumental vars
    //Class constructor
    dgp_OLS_IV(double, colvec, mat, mat, int);
};

dgp_OLS_IV::dgp_OLS_IV(double b, colvec p, mat Ve, 
              mat Vz, int n){
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
  
  colvec stdnorm_errors = rnorm(n * 2);
  mat errors = trans(chol(Ve) * reshape(stdnorm_errors, 2, n));
  colvec stdnorm_z = rnorm(n * n_z);
  
  z = trans(chol(Vz) * reshape(stdnorm_z, n_z, n));
  x = z * p + errors.col(1); //Remember: zero indexing!
  y = b * x + errors.col(0);   
}






class fmsc_OLS_IV {
/*--------------------------------------------------
# Class for FMSC calculations in OLS vs IV example
#-------------------------------------------------*/
  public:
    //class constructor
    fmsc_OLS_IV(const colvec&, const colvec&, 
          const mat&);
    //member functions
    double Tfmsc();       //return FMSC "test statistic"
    double b_ols();       //return OLS estimate
    double b_tsls();      //return tsls estimate
    double b_fmsc();      //return FMSC-selected estimate
    double b_DHW(double); //return DHW pre-test estimate
    double b_AVG();       //return feasible averaging estimate
    rowvec CI_ols(double);  //Confidence interval (CI) for ols
    rowvec CI_tsls(double); //CI for tsls
    rowvec CI_fmsc_naive(double); //Naive CI post-fmsc
    void draw_CI_sims(int); //Initialize CI_sims for use by other members
    rowvec CI_fmsc_1step(double); //Simulation-based CI for 
                    //b_fmsc that simply plugs in tau-hat
    rowvec CI_fmsc_correct(double, double, int); //Simulation-based CI for
                    //b_fmsc based on the two-step procedure using CI for tau
    rowvec CI_AVG_1step(double); //Simulation-based CI for 
                    //b_AVG that simply plugs in tau-hat 
    rowvec CI_AVG_correct(double, double, int); //Simulation-based CI for
                   //b_AVG based on the two-step procedure using CI for tau
  private:
    int n_z, n;
    double xx, g_sq, s_e_sq_ols, s_e_sq_tsls, s_x_sq, 
        s_v_sq, tau, tau_var, ols_estimate, tsls_estimate;
    colvec ols_resid, tsls_resid, first_stage_coefs, zx;
    mat zz_inv, zz, CI_sims;
    //private member functions
    rowvec CI_tau(double); //CI for bias parameter tau
    rowvec CI_Lambda_fmsc(double, double); //Simulation-based CI for
                    //Lambda post-FMSC evaluated at particular value of tau
    rowvec CI_Lambda_AVG(double, double); //Simulation-based CI for
                      //Lambda based on the averaging estimator evaluated 
                      //at particular value of tau
};
  

fmsc_OLS_IV::fmsc_OLS_IV(const colvec& x, const colvec& y, 
                          const mat& z){
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

  xx = dot(x, x);
  ols_estimate = dot(x, y) / xx;
  ols_resid = y - x * ols_estimate;
  
  first_stage_coefs = solve(z,x);
  tsls_estimate = as_scalar(solve(z * first_stage_coefs, y)); 
  tsls_resid = y - x * tsls_estimate; 
  
  mat Qz, Rz, Rzinv;
  qr_econ(Qz, Rz, z);
  Rzinv = inv(trimatu(Rz));
  zz_inv = Rzinv * trans(Rzinv);
  zz = trans(trimatu(Rz)) * trimatu(Rz);
  zx = trans(z) * x;
  g_sq = as_scalar(trans(zx) * zz_inv * zx) / n;
  
  s_e_sq_ols = dot(ols_resid, ols_resid) / n;
  s_e_sq_tsls = dot(tsls_resid, tsls_resid) / n;
  s_x_sq = xx / n;
  s_v_sq = s_x_sq - g_sq;
  tau = dot(x, tsls_resid) / sqrt(n);
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
  //The arguments in the preceding are (q, df, lower.tail, log.p)
  //For details of the C internals, enter qchisq at the R command line
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

rowvec fmsc_OLS_IV::CI_ols(double alpha){
//Member function of class fmsc_OLS_IV
//Returns confidence interval (lower, upper) for OLS estimator
//Arguments: alpha = significance level (0.05 = 95% CI)
  double z_quantile = R::qnorm(1 - alpha/2, 0, 1, 1, 0);
  double SE_ols = sqrt(s_e_sq_ols / (n * s_x_sq));
  double lower = ols_estimate - z_quantile * SE_ols;
  double upper = ols_estimate + z_quantile * SE_ols;
  
  rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}

rowvec fmsc_OLS_IV::CI_tsls(double alpha){
//Member function of class fmsc_OLS_IV
//Returns confidence interval (lower, upper) for TSLS estimator
//Arguments: alpha = significance level (0.05 = 95% CI)
  double z_quantile = R::qnorm(1 - alpha/2, 0, 1, 1, 0);
  double SE_tsls = sqrt(s_e_sq_tsls / (n * g_sq));
  double lower = tsls_estimate - z_quantile * SE_tsls;
  double upper = tsls_estimate + z_quantile * SE_tsls;
  
  rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}

rowvec fmsc_OLS_IV::CI_fmsc_naive(double alpha){
//Member function of class fmsc_OLS_IV
//Returns naive confidence interval (lower, upper)
//for post-FMSC estimator
//Arguments: alpha = significance level (0.05 = 95% CI)
    rowvec out(2);
    if(Tfmsc() < 2){
      out = CI_ols(alpha);
    } else {
      out = CI_tsls(alpha);
    }
  return(out);
}

rowvec fmsc_OLS_IV::CI_tau(double delta){
//Member function of class fmsc_OLS_IV
//Returns CI (lower, upper) for bias parameter tau
//Arguments: delta = significance level (0.05 = 95% CI)
  double z_quantile = R::qnorm(1 - delta/2, 0, 1, 1, 0);
  double SE_tau = sqrt(tau_var);
  double lower = tau - z_quantile * SE_tau;
  double upper = tau + z_quantile * SE_tau;
  
  rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}


void fmsc_OLS_IV::draw_CI_sims(int n_sims){
//Member function of class fmsc_OLS_IV
//Initializes the matrix CI_sims for use by other member functions
  RNGScope scope;
  
  //Construct variance matrix of M
  mat Omega(n_z + 1, n_z + 1);
  Omega(0, 0) = s_x_sq;
  Omega(0, span(1, n_z)) = trans(zx) / n;
  Omega(span(1, n_z), 0) = zx / n;
  Omega(span(1, n_z), span(1, n_z)) = zz / n;
  Omega = s_e_sq_tsls * Omega;
  
  //Construct matrix D that multiplies M
  mat D(3, n_z + 1);
  D(0, 0) = 1 / s_x_sq;
  D(0, span(1, n_z)) = zeros<rowvec>(n_z);
  D(1, 0) = 0;
  D(1, span(1, n_z)) = trans(first_stage_coefs) / g_sq;
  D(2, 0) = 1;
  D(2, span(1, n_z)) = -1 * s_x_sq * 
                      trans(first_stage_coefs) / g_sq;
  
  //Draw simulations
  colvec stdnorm = rnorm(n_sims * Omega.n_rows);
  CI_sims = D * chol(Omega) * reshape(stdnorm, Omega.n_rows, n);
}


rowvec fmsc_OLS_IV::CI_Lambda_fmsc(double alpha, double tau_star){
//Member function of class fmsc_OLS_IV
//Constructs a (1 - alpha) * 100% CI for Lambda at a given value of tau 
//based on FMSC selection. This function assumes that draw_CI_sims has
//already been called so that the data member CI_sims is available.
  rowvec one_vec = ones<rowvec>(CI_sims.n_cols);
  rowvec A = (tau_star / s_x_sq) * one_vec + CI_sims.row(0);
  rowvec B = CI_sims.row(1);
  rowvec C = tau_star * one_vec + CI_sims.row(2);
  
  rowvec w = conv_to<rowvec>
                    ::from(pow(C, 2) < (2 * tau_var * one_vec));
                    
  colvec Lambda = conv_to<colvec>
                          ::from(w % A + (one_vec - w) % B);
  double lower = sample_quantile(Lambda, alpha/2);
  double upper = sample_quantile(Lambda, 1 - alpha/2);
  rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}


rowvec fmsc_OLS_IV::CI_fmsc_1step(double alpha){
//Member function of class fmsc_OLS_IV
//Returns 1-step corrected confidence interval (lower, upper)
//This function assumes that draw_CI_sims has already been called
  rowvec Lambda_interval = CI_Lambda_fmsc(alpha, tau);
  double Lambda_lower = Lambda_interval(0);
  double Lambda_upper = Lambda_interval(1);
  double lower = b_fmsc() - Lambda_upper / sqrt(n);
  double upper = b_fmsc() - Lambda_lower / sqrt(n);
  rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}

rowvec fmsc_OLS_IV::CI_fmsc_correct(double alpha, 
                                  double delta, int n_grid){
//Member function of class fmsc_OLS_IV
//Returns corrected confidence interval (lower, upper)
//for post-FMSC estimator with asymptotic coverage probability
//of at least 1 - (alpha + delta)
//This function assumes that draw_CI_sims has already been called
//Arguments: 
//  alpha       significance level for Lambda CI conditional on tau
//  delta       significance level for tau confidence interval
  rowvec tau_interval = CI_tau(delta);
  double tau_lower = tau_interval(0);
  double tau_upper = tau_interval(1);
  colvec tau_star = linspace(tau_lower, tau_upper, n_grid);
  mat Lambda_CIs(n_grid, 2);
  
  for(int i = 0; i < n_grid; i++){
    Lambda_CIs.row(i) = CI_Lambda_fmsc(alpha, tau_star(i));
  }
  
  double Lambda_lower_min = min(Lambda_CIs.col(0));
  double Lambda_upper_max = max(Lambda_CIs.col(1));
  
  double lower = b_fmsc() - Lambda_upper_max / sqrt(n);
  double upper = b_fmsc() - Lambda_lower_min / sqrt(n);
  rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}


rowvec fmsc_OLS_IV::CI_Lambda_AVG(double alpha, double tau_star){
//Member function of class fmsc_OLS_IV
//Constructs a (1 - alpha) * 100% CI for Lambda at a given value of tau 
//based on feasible AMSE averaging. This function assumes draw_CI_sims has
//already been called so that the data member CI_sims is available.
  rowvec one_vec = ones<rowvec>(CI_sims.n_cols);
  rowvec A = (tau_star / s_x_sq) * one_vec + CI_sims.row(0);
  rowvec B = CI_sims.row(1);
  rowvec C = tau_star * one_vec + CI_sims.row(2);
  
  rowvec sq_bias_est = (pow(C, 2)  - tau_var * one_vec) 
                                / pow(s_x_sq, 2);
  rowvec zero_vec = zeros<rowvec>(CI_sims.n_cols);
  sq_bias_est = max(zero_vec, sq_bias_est);
  
  double var_diff = s_e_sq_tsls * (1 / g_sq - 1 / s_x_sq);
  rowvec w = one_vec / (one_vec + sq_bias_est / var_diff);
  colvec Lambda = conv_to<colvec>
                          ::from(w % A + (one_vec - w) % B);
                          
  double lower = sample_quantile(Lambda, alpha/2);
  double upper = sample_quantile(Lambda, 1 - alpha/2);
  rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}

rowvec fmsc_OLS_IV::CI_AVG_1step(double alpha){
//Member function of class fmsc_OLS_IV
//Returns 1-step corrected confidence interval (lower, upper)
//This function assumes that draw_CI_sims has already been called
  rowvec Lambda_interval = CI_Lambda_AVG(alpha, tau);
  double Lambda_lower = Lambda_interval(0);
  double Lambda_upper = Lambda_interval(1);
  double lower = b_AVG() - Lambda_upper / sqrt(n);
  double upper = b_AVG() - Lambda_lower / sqrt(n);
  rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}


rowvec fmsc_OLS_IV::CI_AVG_correct(double alpha, 
                                double delta, int n_grid){
//Member function of class fmsc_OLS_IV
//Returns corrected confidence interval (lower, upper)
//for feasible averaging estimator with asymptotic coverage probability
//of at least 1 - (alpha + delta)
//This function assumes that draw_CI_sims has already been called
//Arguments: 
//  alpha       significance level for Lambda CI conditional on tau
//  delta       significance level for tau confidence interval
  rowvec tau_interval = CI_tau(delta);
  double tau_lower = tau_interval(0);
  double tau_upper = tau_interval(1);
  colvec tau_star = linspace(tau_lower, tau_upper, n_grid);
  mat Lambda_CIs(n_grid, 2);
  
  for(int i = 0; i < n_grid; i++){
    Lambda_CIs.row(i) = CI_Lambda_AVG(alpha, tau_star(i));
  }
  
  double Lambda_lower_min = min(Lambda_CIs.col(0));
  double Lambda_upper_max = max(Lambda_CIs.col(1));
  
  double lower = b_AVG() - Lambda_upper_max / sqrt(n);
  double upper = b_AVG() - Lambda_lower_min / sqrt(n);
  rowvec out(2);
  out(0) = lower;
  out(1) = upper;
  return(out);
}







// [[Rcpp::export]]
NumericVector mse_compare_cpp(double b, colvec p, mat Ve, 
              mat Vz, int n, int n_reps){
//Function to run n_reps of the simulation study and calculate the MSE
//of various estimators

  colvec ols(n_reps);
  colvec tsls(n_reps);
  colvec fmsc(n_reps);
  colvec DHW90(n_reps);
  colvec DHW95(n_reps);
  colvec AVG(n_reps);

  for(int i = 0; i < n_reps; i++){
    
    dgp_OLS_IV sim_data(b, p, Ve, Vz, n);
    fmsc_OLS_IV sim_results(sim_data.x, sim_data.y, sim_data.z);

    ols(i) = sim_results.b_ols();
    tsls(i) = sim_results.b_tsls();
    fmsc(i) = sim_results.b_fmsc();
    DHW90(i) = sim_results.b_DHW(0.1);
    DHW95(i) = sim_results.b_DHW(0.05);
    AVG(i) = sim_results.b_AVG();
    
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
//Runs simulation once with default values for "uninteresting" params.

//In this specification there are three instruments (z1, z2, z3)
//which ensures that the finite-sample RMSE of the 2SLS estimator exists.
//Further, the design sets
//  var(x) = 1
//  corr(e,v) = r
//  cor(x, z1 + z2 + z3) = p
 
  double b = 0.5;
  colvec p_vec = p * ones(3);
  mat Vz = eye(3, 3) / 3;
  
  mat Ve;
  Ve << 1 << r * sqrt(1 - pow(p,2)) << endr
     << r * sqrt(1 - pow(p,2)) << 1 - pow(p, 2) << endr;
  
  NumericVector out = mse_compare_cpp(b, p_vec, Ve, Vz, n, n_reps);
  return(out);  
}




// [[Rcpp::export]]
mat test_CIs_cpp(double p , double r, int n, 
                                        int n_reps){
//Function to test the confidence interval code with the same parameter
//values as in mse_compare_default_cpp
//In this specification there are three instruments (z1, z2, z3)
//which ensures that the finite-sample RMSE of the 2SLS estimator exists.
//Further, the design sets
//  var(x) = 1
//  corr(e,v) = r
//  cor(x, z1 + z2 + z3) = p
 
  double b = 0.5;
  colvec p_vec = p * ones(3);
  mat Vz = eye(3, 3) / 3;
  
  mat Ve;
  Ve << 1 << r * sqrt(1 - pow(p,2)) << endr
     << r * sqrt(1 - pow(p,2)) << 1 - pow(p, 2) << endr;

  mat out(n_reps, 2, fill::zeros);
  
  for(int i = 0; i < n_reps; i++){
    
    dgp_OLS_IV sim_data(b, p_vec, Ve, Vz, n);
    fmsc_OLS_IV sim_results(sim_data.x, sim_data.y, sim_data.z);
    
    sim_results.draw_CI_sims(1000);
    out.row(i) = sim_results.CI_AVG_correct(0.02, 0.08, 100);
    //out.row(i) = sim_results.CI_AVG_1step(0.05);
 
  }
  return(out);  
}
