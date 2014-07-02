/*-------------------------------------------------------
# Filename:        functions.cpp
# Author:          Frank DiTraglia
# First Version:   2014-02-12
# Last Updated:    2014-07-02 
#--------------------------------------------------------
# C++ code to calculate FMSC for the empirical example 
# and interface with R.
#-------------------------------------------------------*/


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//Class for carrying out two-stage least squares with various options for
//variance matrix estimation
class tsls_fit {
  public:
    tsls_fit(const mat&, const colvec&, const mat&);
    colvec est() {return(b);} 
    colvec SE_textbook() {return(sqrt(diagvec(V_textbook() / n)));}
    colvec SE_robust() {return(sqrt(diagvec(V_robust() / n)));}
    colvec SE_center() {return(sqrt(diagvec(V_center() / n)));}
    colvec resid() {return(residuals);}
    mat Omega_textbook() {return(s_sq * Rz.t() * Rz / n);}
    mat Omega_robust() {return(Z_copy.t() * D * Z_copy / n);}
    mat Omega_center() {return(Z_copy.t() * (D / n -  (residuals * residuals.t()) / (n * n)) * Z_copy);}
    mat V_textbook() {return(n * s_sq * Rtilde_inv * Rtilde_inv.t());}
    mat V_robust() {return(n * n * C * Omega_robust() * C.t());}
    mat V_center() {return(n * n * C * Omega_center() * C.t());}
    mat C;
  private:
    int n, k;
    double s_sq;
    colvec b, residuals;
    mat Z_copy, Qz, Rz, Qtilde, Rtilde, Rtilde_inv, D;
};
//Class constructor
tsls_fit::tsls_fit(const mat& X, const colvec& y, const mat& Z){
  n = y.n_elem;
  k = X.n_cols;
  Z_copy = Z; 
  qr_econ(Qz, Rz, Z);
  qr_econ(Qtilde, Rtilde, Qz.t() * X);
  b = solve(trimatu(Rtilde), Qtilde.t() * Qz.t() * y);
  residuals = y - X * b;
  D =  diagmat(pow(residuals, 2));
  s_sq = dot(residuals, residuals) / (n - k);
  Rtilde_inv = inv(trimatu(Rtilde)); 
  C = Rtilde_inv * Qtilde.t() * inv(trimatl(Rz.t()));
}


//Class for FMSC calculations to choose IVs for 2SLS estimation.
//As implemented, it only allows target parameters that are weighted 
//averages of the elements of beta. A more general solution would be 
//to pass a function pointer instead of a vector of weights.
class fmsc_chooseIV {
  public:
    fmsc_chooseIV(const mat&, const colvec&, const mat&, const mat&, umat); 
    //Extract Target parameter estimators 
    //  All candidate specifications:
    colvec mu(colvec weights){return(estimates.t() * weights);}
    //  Valid model only:
    double mu_valid(colvec weights){
      colvec mu_estimates = mu(weights);
      return(mu_estimates(0));
    }
    //  Full model only
    double mu_full(colvec weights){
      colvec mu_estimates = mu(weights);
      return(mu_estimates(mu_estimates.n_elem - 1));
    }
    //Squared Asymptotic Bias of mu estimators 
    //  (all candidate specifications)
    colvec abias_sq(colvec weights){
      colvec out(z2_indicators.n_cols);
      for(int i = 0; i < K.n_elem; i++){
        out(i) = as_scalar(weights.t() * sqbias_inner.slice(i) * weights);
      }
      return(out);
    }
    //Asymptotic Variance of mu estimators
    //  (all candidate specifications)
    colvec avar(colvec weights){
      colvec out(z2_indicators.n_cols);
      for(int i = 0; i < K.n_elem; i++){
        out(i) = as_scalar(weights.t() * avar_inner.slice(i) * weights);
      }
      return(out);
    }
    //Plain vanilla FMSC
    //  (allows negative squared bias estimate)
    colvec fmsc(colvec weights){
      return(abias_sq(weights) + avar(weights));
    }
    //Positive-Part FMSC
    //  (set negative squared bias estimates to zero)
    colvec fmsc_pos(colvec weights){
      colvec first_term = abias_sq(weights);
      first_term = max(first_term, zeros<colvec>(first_term.n_elem));
      colvec second_term = avar(weights);
      return(first_term + second_term);
    }
    //FMSC-selected estimator of mu
    double mu_fmsc(colvec weights){
      colvec criterion_values = fmsc(weights);
      uword which_min;
      criterion_values.min(which_min);
      colvec b_fmsc = estimates.col(which_min);
      return(dot(weights, b_fmsc));
    }    
    //positive-part FMSC-selected estimator of mu
    double mu_fmsc_pos(colvec weights){
      colvec criterion_values = fmsc_pos(weights);
      uword which_min;
      criterion_values.min(which_min);
      colvec b_fmsc = estimates.col(which_min);
      return(dot(weights, b_fmsc));
    }   
    tsls_fit valid, full;
    colvec tau;
    mat Psi, tau_outer_est, Bias_mat, estimates;
    umat z2_indicators; 
    field<mat> K, Omega;
    cube sqbias_inner, avar_inner; 
    int n_obs, n_z1, n_z2, n_z, n_params;
};
//Class constructor - initialization list ensures tsls_fit constructor 
//is called before entering body of the present constuctor 
fmsc_chooseIV::fmsc_chooseIV(const mat& x, const colvec& y, const mat& z1, 
           const mat& z2, umat candidates = zeros<umat>(1,1)): 
                valid(x, y, z1), full(x, y, join_rows(z1, z2)){
    //If specified, each column of candidates is an indicator vector 
    //that corresponds to the columns of z2 used in estimation for that 
    //candidate. We always calculate the Valid and Full. The default
    //is to calculate *only* these.
    n_z1 = z1.n_cols;
    n_z2 = z2.n_cols;
    n_z = n_z1 + n_z2;
    n_obs = y.n_elem;
    n_params = x.n_cols;
    uvec valid_indicator = zeros<uvec>(n_z2);
    uvec full_indicator = ones<uvec>(n_z2);
    uvec z1_indicator = ones<uvec>(n_z1); //Always include z1
    
    mat K_valid = n_obs * valid.C;
    mat Omega_valid = valid.Omega_robust(); //no centering needed
    mat K_full = n_obs * full.C;
    mat Omega_full = full.Omega_center();
    Psi =  join_rows(-1 * z2.t() * x * K_valid / n_obs, eye(n_z2, n_z2));
    tau = z2.t() * valid.resid() / sqrt(n_obs);
    tau_outer_est = tau * tau.t() - Psi * Omega_full * Psi.t();
    Bias_mat = mat(n_z, n_z, fill::zeros);
    Bias_mat(span(n_z1, n_z - 1), span(n_z1, n_z - 1)) = tau_outer_est;
    
    //Are there any additional candidates besides valid and full?
    int n_add_cand;
    if(all(vectorise(candidates) == 0)){
      n_add_cand = 0;
      z2_indicators = join_rows(valid_indicator, full_indicator);
    }else{
      n_add_cand = candidates.n_cols;
      z2_indicators = join_rows(valid_indicator, candidates);
      z2_indicators = join_rows(z2_indicators, full_indicator);
    }

    //Vector to store current candidate specification
    //Used to subset *full* instrument matrix, *not* z2
    uvec candidate(n_z);    
    //Matrix to store estimate from each candidate
    mat estimates_temp(n_params, n_add_cand + 2);
    //Temporary storage for K(S) and Omega(S) 
    //Dims vary by candidate so store in a field
    field<mat> K_temp(n_add_cand + 2);
    field<mat> Omega_temp(n_add_cand + 2);
    //Matrix to store Xi(S) * Bias_mat * Xi(S)'
    //as we loop over candidate moment sets
    mat inner;
    //Temp storage for K(S)Xi(S) * Bias_mat * Xi(S)'K(S)'
    //and K(S)Xi(S) * Omega * K(S)'Xi(S)'. Since dims are
    //the same for each candidate, store results in a cube
    cube sqbias_inner_temp(n_params, n_params, n_add_cand + 2);
    cube avar_inner_temp(n_params, n_params, n_add_cand + 2);
    
    //Results for valid estimator - outside loop since already fitted
    K_temp(0) = K_valid;
    Omega_temp(0) = Omega_valid;
    estimates_temp.col(0) = valid.est();
    //Use find to convert vector of 0-1 indicators into
    //vector of *positions* of the non-zeros
    candidate = find(join_cols(z1_indicator, valid_indicator));
    inner = Bias_mat.submat(candidate, candidate);
    sqbias_inner_temp.slice(0) = K_valid * inner * K_valid.t();
    avar_inner_temp.slice(0) =  K_valid * Omega_valid * K_valid.t();
    
    //Results for any candidates besides valid and full
    for(int i = 0; i < n_add_cand; i++){
      //Fit 2SLS for current candidate
      mat z2_candidate = z2.cols(find(candidates.col(i)));
      tsls_fit candidate_fit(x, y, join_rows(z1, z2_candidate));
      //Extract and store quantities needed for FMSC
      mat K_candidate = n_obs * candidate_fit.C;
      mat Omega_candidate = candidate_fit.Omega_center();
      K_temp(i + 1) = K_candidate;
      Omega_temp(i + 1) = Omega_candidate;
      estimates_temp.col(i + 1) = candidate_fit.est();
      //Use find to convert vector of 0-1 indicators into
      //vector of *positions* of the non-zeros
      candidate = find(join_cols(z1_indicator, candidates.col(i)));
      inner = Bias_mat.submat(candidate, candidate);
      sqbias_inner_temp.slice(i + 1) = K_candidate * inner * 
                                                     K_candidate.t();
      avar_inner_temp.slice(i + 1) =  K_candidate * Omega_candidate * 
                                                    K_candidate.t();
    }
    
    //Results for full estimator - outside loop since already fitted
    K_temp(n_add_cand + 1) = K_full;
    Omega_temp(n_add_cand + 1) = Omega_full;
    estimates_temp.col(n_add_cand + 1) = full.est();
    //Use find to convert vector of 0-1 indicators into
    //vector of *positions* of the non-zeros
    candidate = find(join_cols(z1_indicator, full_indicator));
    inner = Bias_mat.submat(candidate, candidate);
    sqbias_inner_temp.slice(n_add_cand + 1) = K_full * inner * K_full.t();
    avar_inner_temp.slice(n_add_cand + 1) = K_full * Omega_full * 
                                                     K_full.t();
    //Store results                                     
    K = K_temp;
    Omega = Omega_temp;
    estimates = estimates_temp;
    sqbias_inner = sqbias_inner_temp;
    avar_inner = avar_inner_temp;
}


//Bare-bones R interface to the results of the class constructor
//for the fmsc_chooseIV class
// [[Rcpp::export]]
List fmsc_2SLS_cpp(mat x, colvec y, mat z1, mat z2, 
                                    umat candidates){
  
  fmsc_chooseIV results(x, y, z1, z2, candidates);  
  List out = List::create(results.estimates,
                          results.z2_indicators,
                          results.sqbias_inner,
                          results.avar_inner,
                          results.tau,
                          results.Psi,
                          results.tau_outer_est,
                          results.K,
                          results.Omega);
  out.names() = CharacterVector::create("est",
                                        "cand",
                                        "sqbias",
                                        "avar",
                                        "tau",
                                        "Psi",
                                        "tau.outer",
                                        "K",
                                        "Omega");
  return(out);
}
