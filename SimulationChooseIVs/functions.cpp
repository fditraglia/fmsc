/*-------------------------------------------------------
# Filename:        functions.cpp
# Author:          Frank DiTraglia
# First Version:   2014-02-12
# Last Updated:    2014-06-27
#--------------------------------------------------------
# C++ functions for choosing IVs simulation experiment.
#-------------------------------------------------------*/

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat mvrnorm_cpp(int n, vec mu, mat Sigma){
/*-------------------------------------------------------
# Generate draws from a multivariate normal distribution
#--------------------------------------------------------
#  n        number of samples
#  mu       mean vector
#  Sigma    covariance matrix
#--------------------------------------------------------
# Details:
#           This is essentially a stripped-down version
#           of the mvrnorm function from the MASS library
#           in R. Through the magic of Rcpp we're 
#           transforming the *same* standard normal draws
#           as the R version. However, since Armadillo
#           follows a different convention from R in its
#           definition of the eign-decomposition, the 
#           output of this function will *not* be the
#           same as that of its R counterpart. Since we
#           access R's function for generating normal
#           draws, we can set the seed from R.
#-------------------------------------------------------*/
  RNGScope scope;
  int p = Sigma.n_cols;
  mat X = reshape(vec(rnorm(p * n)), p, n);
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Sigma);
  X = eigvec * diagmat(sqrt(eigval)) * X;
  X.each_col() += mu;
  return(X.t());
}

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



//A C++ version of R's cancor function for canonical correlation analysis
class cancor {
  public: 
    cancor(const mat&, const mat&);
    rowvec cor;
    mat xcoef, ycoef;
  private:
    colvec d;
    mat Qx, Rx, Qy, Ry, U, V;
};
//Class constructor
cancor::cancor(const mat& X, const mat& Y){
  qr_econ(Qx, Rx, X - repmat(mean(X), X.n_rows, 1));
  qr_econ(Qy, Ry, Y - repmat(mean(Y), Y.n_rows, 1));
  svd_econ(U, d, V, Qx.t() * Qy);
  xcoef = solve(trimatu(Rx), U);
  ycoef = solve(trimatu(Ry), V);
  cor = d.t();
}

//Class for calculating the canonical correlation information criterion (CCIC)
//of Hall & Peixe (2003). The argument Z_c is the instrument matrix, while 
//D_resid is a matrix containing the partial derivatives of the GMM residual
//function u(theta) with respect to the parameter vector theta. 
//For linear GMM, D_resid is the matrix of regressors, X.
class CCIC {
  public:
    CCIC(const mat&, const mat&);
    double BIC(){return(first_term + n_overid * log(double(n_obs)));}
    double AIC(){return(first_term + n_overid * 2.0);}
    double HQ(){return(first_term + n_overid * 2.01 * log(log(double(n_obs))));}
  private:
    int n_obs;
    int n_overid;
    double first_term;
    cancor cc_results;
    vec r;
};
//Class constructor
CCIC::CCIC(const mat& D_resid, const mat& Z): cc_results(D_resid, Z){
  n_obs = Z.n_rows;
  n_overid = Z.n_cols - D_resid.n_cols;
  r = cc_results.cor;
  first_term = double(n_obs) * sum(log(ones(r.n_elem) - pow(r, 2)));
}


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


//Class for calculating Andrews (1999) GMM moment selection
//criteria in linear GMM models 
class linearGMM_msc {
  public:
    linearGMM_msc(const mat&, const colvec&, const mat&);
    double Jstat(){return(J);}
    double pJtest(){return(R::pchisq(J, n_overid, 1, 0));}
    //the arguments of R::pchisq are (q, df, lower.tail, log.p)
    double GMM_AIC(){return(J - 2 * n_overid);}
    double GMM_BIC(){return(J - log(n_obs) * n_overid);}
    double GMM_HQ(){return(J - 2.01 * log(log(n_obs)) * n_overid);}
    colvec est_1step(){return(b_1step);}
    colvec est_2step(){return(b_2step);}
  private:
    tsls_fit first_step;
    colvec b_1step, b_2step, resid_2step, v;
    mat L, R, Q, Omega_tilde, L_tilde;
    double J;
    int n_obs, n_overid;
};
//Class constructor
//Uses a initialization list to run the first-step
//of the efficient two-step GMM estimation procedure
linearGMM_msc::linearGMM_msc(const mat& X, const colvec& y,
                             const mat& Z): first_step(X, y, Z){
  n_obs = X.n_rows;
  n_overid = Z.n_cols - X.n_cols;
  b_1step = first_step.est();
  L = trans(chol(first_step.Omega_center()));
  qr_econ(Q, R, solve(trimatl(L), Z.t() * X));
  b_2step = solve(trimatu(R), Q.t() * 
                    solve(trimatl(L), Z.t() * y));
  resid_2step = y - X * b_2step;
  //Centered robust var matrix using second step residuals
  Omega_tilde = Z.t() * ( diagmat(pow(resid_2step, 2)) / n_obs   
       - ( resid_2step * resid_2step.t() ) / (n_obs * n_obs) ) * Z;
  L_tilde = trans(chol(Omega_tilde));
  v = solve(trimatl(L_tilde), Z.t() * resid_2step);
  J = dot(v, v) / n_obs;
}



//Class to carry out moment selection for linear GMM models
//using the Andrews (1999) criteria and CCIC of Hall & Peixe (2003)
//This class is really just a wrapper to repeatedly use the 
//linearGMM_msc class to calculate moment selection criteria over a 
//collection of candidate moment sets.
class linearGMM_select{
  public:
    linearGMM_select(const mat&, const colvec&, 
                                     const mat&, const umat&);
    colvec est_AIC(){return(est_selected(AIC, estimates_1step));}
    colvec est_BIC(){return(est_selected(BIC, estimates_1step));}
    colvec est_HQ(){return(est_selected(HQ, estimates_1step));}
    colvec est_CCIC_AIC(){return(est_selected(AIC_CCIC, estimates_1step));}
    colvec est_CCIC_BIC(){return(est_selected(BIC_CCIC, estimates_1step));}
    colvec est_CCIC_HQ(){return(est_selected(HQ_CCIC, estimates_1step));}
    uvec moments_AIC(){return(moments_selected(AIC));}
    uvec moments_BIC(){return(moments_selected(BIC));}
    uvec moments_HQ(){return(moments_selected(HQ));}
    uvec moments_CCIC_AIC(){return(moments_selected(AIC_CCIC));}
    uvec moments_CCIC_BIC(){return(moments_selected(BIC_CCIC));}
    uvec moments_CCIC_HQ(){return(moments_selected(HQ_CCIC));}
    //Public data members
    colvec J, pJtest, AIC, BIC, HQ, AIC_CCIC, BIC_CCIC, HQ_CCIC;
    mat estimates_1step, estimates_2step;
  private:
    int n_candidates, n_params;
    umat moment_sets_copy;
    colvec est_selected(colvec criterion, mat est){
        uword which_min;
        criterion.min(which_min);
        return(est.col(which_min));
    }
    uvec moments_selected(colvec criterion){
        uword which_min;
        criterion.min(which_min);
        return(moment_sets_copy.col(which_min));
    }
};
//Class constructor
linearGMM_select::linearGMM_select(const mat& X,
      const colvec& y, const mat& Z_full, const umat& moment_sets){
  //Each column of the matrix moment_sets is a moment set, namely
  //a vector of zeros and ones that indicates which columns of Z 
  //to use in estimation  
  n_candidates = moment_sets.n_cols;
  n_params = X.n_cols;
  moment_sets_copy = moment_sets;
  J.zeros(n_candidates);
  pJtest.zeros(n_candidates);
  AIC.zeros(n_candidates);
  BIC.zeros(n_candidates);
  HQ.zeros(n_candidates);
  AIC_CCIC.zeros(n_candidates);
  BIC_CCIC.zeros(n_candidates);
  HQ_CCIC.zeros(n_candidates);
  estimates_1step.zeros(n_params, n_candidates);
  estimates_2step.zeros(n_params, n_candidates);
  
  for(int i = 0; i < n_candidates; i++){
    //Andrews (1999) Criteria
    uvec cand_ind = find(moment_sets.col(i)); //indices of non-zeros
    linearGMM_msc candidate(X, y, Z_full.cols(cand_ind));
    J(i) = candidate.Jstat();
    pJtest(i) = candidate.pJtest();
    AIC(i) = candidate.GMM_AIC();
    BIC(i) = candidate.GMM_BIC();
    HQ(i) = candidate.GMM_HQ();
    //Parameter Estimates
    estimates_1step.col(i) = candidate.est_1step();
    estimates_2step.col(i) = candidate.est_2step();
    //Hall & Peixe (2003) Criteria
    CCIC candidate_CCIC(X, Z_full.cols(cand_ind));
    AIC_CCIC(i) = candidate_CCIC.AIC();
    BIC_CCIC(i) = candidate_CCIC.BIC();
    HQ_CCIC(i) = candidate_CCIC.HQ();
  }
}


//This class only allows target parameters that are weighted averages
//of the elements of beta. A more general solution would be to pass
//a function pointer instead of a vector of weights.
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
    
    mat K_valid = -1.0 * n_obs * valid.C;
    mat Omega_valid = valid.Omega_robust(); //no centering needed
    mat K_full = -1.0 * n_obs * full.C;
    mat Omega_full = full.Omega_center();
    Psi =  join_rows(z2.t() * x * K_valid / n_obs, eye(n_z2, n_z2));
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
      mat K_candidate = -1.0 * n_obs * candidate_fit.C;
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
    //Store results in private data members                                       
    K = K_temp;
    Omega = Omega_temp;
    estimates = estimates_temp;
    sqbias_inner = sqbias_inner_temp;
    avar_inner = avar_inner_temp;
}


//Class for post-FMSC CI construction in the case of a single
//endogenous regressor (whose coefficient is the target
//parameter), a single "suspect" instrument and any number of
//valid instruments. In this case, tau is a scalar.
class fmsc_CI_simple{
  public:
    fmsc_CI_cimple(const colvec&, const colvec& const mat&, 
                   const colvec&, int);
  private:  
    fmsc_chooseIV fmsc;
    mat Psi, PsiM, Omega;
    double mu_full, mu_valid;
  
  
};
//Class constructor
fmsc_CI_simple::fmsc_CI_simple(const colvec& x, const colvec& y, 
              const mat& z_valid, const colvec& z_suspect, 
              int n_sims = 500): fmsc(x, y, z_valid, z_suspect){
  
  
  

}



class dgp {
  public:
    dgp(double, vec, double, mat, mat, int);
    colvec x, y, w;
    mat z;
  private: 
    int n_z;
    mat e_v_w;
};
//Class constructor
dgp::dgp(double b, vec p, double g, mat V, mat Q, int n){
//b = scalar coef for single endog regressor
//p = vector of first-stage coeffs for exog instruments
//g = scalar first-stage coeff for potentially endog instrument w
//V = variance matrix (3x3) for (epsilon, v, w)'
//Q = variance matrix for exog instruments
//n = sample size
  n_z = Q.n_cols;
  z = mvrnorm_cpp(n, zeros<vec>(n_z), Q);
  e_v_w = mvrnorm_cpp(n, zeros<vec>(3), V);
  w = e_v_w.col(2);
  x = z * p + g * w + e_v_w.col(1);
  y = b * x + e_v_w.col(0);
}



// [[Rcpp::export]]
NumericVector mse_compare_cpp(double b, double g, vec p, mat V, mat Q,
                              int n, int n_reps){
//Function to run n_reps of the simulation study and calculate the MSE
//of various estimators
  
  colvec valid(n_reps);
  colvec full(n_reps);
  colvec fmsc(n_reps);
  colvec fmsc_pos(n_reps);
  colvec j90(n_reps);
  colvec j95(n_reps);
  colvec aic(n_reps);
  colvec bic(n_reps);
  colvec hq(n_reps);
  colvec aic_combine(n_reps);
  colvec bic_combine(n_reps);
  colvec hq_combine(n_reps);
  
  //Don't bother storing CCIC for each replication
  //  (Only need temporary storage to calculate combined estimators)
  double ccic_aic;
  double ccic_bic;
  double ccic_hq;
  
  //Candidate moment sets to pass to linear_GMM_select
  //   (not needed for fmsc_chooseIV since this defaults
  //    to valid and full models only)
  umat valid_full(Q.n_rows + 1, 2, fill::ones);
  valid_full(Q.n_rows, 0) = 0;
  
  //Weights for FMSC
  //  (Coef. on single endog. regressor is the target param.)
  colvec w_one(1, 1, fill::ones);

  for(int i = 0; i < n_reps; i++){
    
    dgp sim(b, p, g, V, Q, n);
    fmsc_chooseIV fmsc_results(sim.x, sim.y, sim.z, sim.w);
    linearGMM_select gmm_msc_results(sim.x, sim.y, 
                          join_rows(sim.z, sim.w), valid_full);
    
    valid(i) = fmsc_results.mu_valid(w_one);
    full(i) = fmsc_results.mu_full(w_one);
    fmsc(i) = fmsc_results.mu_fmsc(w_one);
    fmsc_pos(i) = fmsc_results.mu_fmsc_pos(w_one);
    aic(i) = as_scalar(gmm_msc_results.est_AIC());
    bic(i) = as_scalar(gmm_msc_results.est_BIC());
    hq(i) = as_scalar(gmm_msc_results.est_HQ());
    
    //Downward J-test Estimators
    //  95%
    if(gmm_msc_results.pJtest(1) < 0.05){
      j95(i) = valid(i);
    }else{
      j95(i) = full(i);
    }
    //  90%
    if(gmm_msc_results.pJtest(1) < 0.1){
      j90(i) = valid(i);
    }else{
      j90(i) = full(i);
    }
    
    //Combination CCIC and GMM-MSC Estimators
    //  (Select valid unless *both* GMM-MSC and CCIC choose full)
    ccic_aic = as_scalar(gmm_msc_results.est_CCIC_AIC());
    ccic_bic = as_scalar(gmm_msc_results.est_CCIC_BIC());
    ccic_hq = as_scalar(gmm_msc_results.est_CCIC_HQ());
    //  CCIC-GMM-MSC-AIC
    if(ccic_aic == aic(i) == full(i)){
      aic_combine(i) = full(i);
    }else{
      aic_combine(i) = valid(i);
    }
    //  CCIC-GMM-MSC-BIC
    if(ccic_bic == bic(i) == full(i)){
      bic_combine(i) = full(i);
    }else{
      bic_combine(i) = valid(i);
    }
    //  CCIC-GMM-MSC-HQ
    if(ccic_hq == hq(i) == full(i)){
      hq_combine(i) = full(i);
    }else{
      hq_combine(i) = valid(i);
    }  
    
  }
  
    double const trim_frac = 0; //Change this if you want trimmed MSE
    
    double MSE_valid = MSE_trim(valid, b, trim_frac);
    double MSE_full = MSE_trim(full, b, trim_frac);
    double MSE_fmsc = MSE_trim(fmsc, b, trim_frac);
    double MSE_fmsc_pos = MSE_trim(fmsc_pos, b, trim_frac);
    double MSE_j90 = MSE_trim(j90, b, trim_frac);
    double MSE_j95 = MSE_trim(j95, b, trim_frac);
    double MSE_aic = MSE_trim(aic, b, trim_frac);
    double MSE_bic = MSE_trim(bic, b, trim_frac);
    double MSE_hq = MSE_trim(hq, b, trim_frac);
    double MSE_aic_combine = MSE_trim(aic_combine, b, trim_frac);
    double MSE_bic_combine = MSE_trim(bic_combine, b, trim_frac);
    double MSE_hq_combine = MSE_trim(hq_combine, b, trim_frac);
        
  //Create and return vector of results
    NumericVector out = NumericVector::create(MSE_valid, 
                                              MSE_full,
                                              MSE_fmsc,
                                              MSE_fmsc_pos,
                                              MSE_j90,
                                              MSE_j95,
                                              MSE_aic,
                                              MSE_bic, 
                                              MSE_hq,
                                              MSE_aic_combine,
                                              MSE_bic_combine,
                                              MSE_hq_combine);
    out.names() = CharacterVector::create("Valid", 
                                          "Full",
                                          "FMSC",
                                          "posFMSC",
                                          "J90",
                                          "J95",
                                          "AIC",
                                          "BIC", 
                                          "HQ",
                                          "combAIC",
                                          "combBIC",
                                          "combHQ");
  return(out);
}

// [[Rcpp::export]]
NumericVector mse_compare_default_cpp(double g, double r, int n,
                                      int n_reps){
//This is simply a wrapper to mse_compare_cpp that runs the simulation
//with default values for the "uninteresting parameters."
//The setup is described in Section 3.6 of the paper.
  double b = 0.5;
  colvec p = 1.0 / 3.0 * ones(3);
  mat Q = 1.0 / 3.0 * eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5 - g * r << 8.0 / 9.0 - pow(g, 2.0) << 0 << endr
    << r << 0 << 1 << endr;
    
  NumericVector out = mse_compare_cpp(b, g, p, V, Q, n, n_reps);
  return(out);  
}