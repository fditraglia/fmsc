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


//The member functions as given here only allow for target parameters
//that are a linear function of beta. I've created a branch of this
//project with more general code using function pointers, but it's 
//a little unwieldy for simple examples.
class fmsc_chooseIV {
  public:
    fmsc_chooseIV(const mat&, const colvec&, const mat&, const mat&, umat); 
    //Target parameter estimators for all candidate specs
    colvec mu(colvec weights){return(weights.t() * estimates);}
    //Valid model estimator of target parameter
    double mu_valid(colvec weights){
      colvec mu = mu_est(weights);
      return(mu(0));
    }
    //Full model estimator of target parameter
    double mu_full(colvec weights){
      colvec mu = mu_est(weights);
      return(mu(mu.n_elem - 1));
    }
    //Squared Asymptotic Bias estimates for target 
    //parameter under each candidate specification
    colvec abias_sq(colvec weights){
        colvec out(z2_indicators.n_cols);
        for(int i = 0; i < K.n_elem; i++){
          out(i) = as_scalar(weights.t() * sqbias_inner(i) * weights);
        }
        return(out);
    }
    //Asymptotic Variance estimates for target parameter
    //under each candidate specification
    colvec avar(colvec weights){
      colvec out(z2_indicators.n_cols);
      for(int i = 0; i < K.n_elem; i++){
        out(i) = as_scalar(weights.t() * avar_inner(i) * weights);
      }
      return(out);
    }
    //Calculate fmsc allowing negative squared bias estimate
    colvec fmsc(colvec weights){
      return(abias_sq(weights) + avar(weights));
    }
    //Calculate fmsc setting negative squared bias to zero
    //This is the "positive-part" FMSC
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
//  private:
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
           const mat& z2, umat candidates = zeros(1,1)): 
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
    //Store results in private data members                                       
    K = K_temp;
    Omega = Omega_temp;
    estimates = estimates_temp;
    sqbias_inner = sqbias_inner_temp;
    avar_inner = avar_inner_temp;
}




class dgp {
  public:
    dgp(double, vec, double, mat, mat, int);
    colvec x, y, z2;
    mat z1;
  private: 
    int n_z1;
    mat u_e_z2;
};
//Class constructor
dgp::dgp(double b, vec p, double g, mat V, mat Q, int n){
//b = scalar coef for single endog regressor
//p = vector of first-stage coeffs for exog instruments
//g = scalar first-stage coeff for potentially endog instrument w
//V = variance matrix for (u, epsilon, w)'
//Q = variance matrix for exog instruments
//n = sample size
  RNGScope scope;
  n_z1 = Q.n_cols;
  z1 = trans(chol(Q) * reshape(colvec(rnorm(n * n_z1)), n_z1, n));
  u_e_z2 = trans(chol(V) * reshape(colvec(rnorm(3 * n)), 3, n));
  z2 = u_e_z2.col(2);
  x = z1 * p + g * z2 + u_e_z2.col(1);
  y = b * x + u_e_z2.col(0);
}




//Testing code - Make some of the member functions available to R
// [[Rcpp::export]]
colvec tsls_est_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.est());
}

// [[Rcpp::export]]
colvec tsls_SE_textbook_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.SE_textbook());
}


// [[Rcpp::export]]
colvec tsls_SE_robust_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.SE_robust());
}

// [[Rcpp::export]]
colvec tsls_SE_center_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.SE_center());
}


// [[Rcpp::export]]
colvec test_dgp(double g, double r, int n){
  colvec p = ones(3) / 10;
  double b = 1;
  mat Q = eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5  - g * r << 1 << 0 << endr
    << r << 0 << 1 << endr;
  dgp sims(b, p, g, V, Q, n);
  tsls_fit valid(sims.x, sims.y, sims.z1);
  return(valid.est());
}


// [[Rcpp::export]]
List cancor_cpp(mat X, mat Y){
  cancor results(X, Y);
  return List::create(Named("cor") = results.cor,
                      Named("xcoef") = results.xcoef,
                      Named("ycoef") = results.ycoef);
}



// [[Rcpp::export]]
List dgp_cpp(double g, double r, int n = 500){
  //This is the simulation setup from the original 
  //version of the paper (Section 3.4)
  double b = 1;
  colvec p = 0.1 * ones(3);
  mat Q = eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5  - g * r << 1 << 0 << endr
    << r << 0 << 1 << endr;
  
    dgp sims(b, p, g, V, Q, n);
  
  return List::create(Named("x") = sims.x,
                      Named("y") = sims.y,
                      Named("z1") = sims.z1, 
                      Named("z2") = sims.z2);
}


// [[Rcpp::export]]
colvec CCIC_test(double g, double r, int n = 500){
  
  //This is the simulation setup from the original 
  //version of the paper (Section 3.4)
  double b = 1;
  colvec p = 0.1 * ones(3);
  mat Q = eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5  - g * r << 1 << 0 << endr
    << r << 0 << 1 << endr;
  
  dgp sims(b, p, g, V, Q, n);
  
  CCIC HallPeixe(sims.x, join_rows(sims.z1, sims.z2));
  colvec out(3);
  out(0) = HallPeixe.BIC();
  out(1) = HallPeixe.AIC();
  out(2) = HallPeixe.HQ();
  return(out);
}

// [[Rcpp::export]]
colvec Andrews_test(double g, double r, int n = 500){
  
  //This is the simulation setup from the original 
  //version of the paper (Section 3.4)
  double b = 1;
  colvec p = 0.1 * ones(3);
  mat Q = eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5  - g * r << 1 << 0 << endr
    << r << 0 << 1 << endr;
  
  dgp sims(b, p, g, V, Q, n);
  
  linearGMM_msc Andrews(sims.x, sims.y,
                        join_rows(sims.z1, sims.z2));
  colvec out(7);
  out(0) = as_scalar(Andrews.est_1step());
  out(1) = as_scalar(Andrews.est_2step());
  out(2) = Andrews.Jstat();
  out(3) = Andrews.pJtest();
  out(4) = Andrews.GMM_AIC();
  out(5) = Andrews.GMM_BIC();
  out(6) = Andrews.GMM_HQ();
  return(out);
}


// [[Rcpp::export]]
List GMMselect_test(double g, double r, int n = 500){
  
  //This is the simulation setup from the original 
  //version of the paper (Section 3.4)
  double b = 1;
  colvec p = 0.1 * ones(3);
  mat Q = eye(3, 3);
  mat V(3,3); 
  V << 1 << 0.5 - g * r << r << endr
    << 0.5  - g * r << 1 << 0 << endr
    << r << 0 << 1 << endr;
  
  //Only two candidate specifications
  umat valid_full(4,2);
  valid_full << 1 << 1 << endr
             << 1 << 1 << endr
             << 1 << 1 << endr
             << 0 << 1 << endr;
  
  dgp sims(b, p, g, V, Q, n);
  
  linearGMM_select results(sims.x, sims.y,
                           join_rows(sims.z1, sims.z2),
                           valid_full);
//  colvec testy = conv_to<colvec>::from(results.moments_BIC());          
  return List::create(
//                      Named("x") = sims.x,
//                      Named("y") = sims.y,
//                      Named("z1") = sims.z1,
//                      Named("z2") = sims.z2,
                      Named("J") = results.J,
                      Named("pJtest") = results.pJtest,
                      Named("AIC") = results.AIC,
                      Named("BIC") = results.BIC,
                      Named("HQ") = results.HQ,
//                      Named("AIC_CCIC") = results.AIC_CCIC,
//                      Named("BIC_CCIC") = results.BIC_CCIC,
//                      Named("HQ_CCIC") = results.HQ_CCIC,
                      Named("onestep") = results.estimates_1step,
//                      Named("twostep") = results.estimates_2step,
                      Named("estAIC") = results.est_AIC(),
                      Named("estBIC") = results.est_BIC(), 
                      Named("estHQ") = results.est_HQ(),
//                      Named("estCCIC_AIC") = results.est_CCIC_AIC(),
//                      Named("estCCIC_HQ") = results.est_CCIC_HQ(),
//                      Named("estCCIC_BIC") = results.est_CCIC_BIC(),
                      Named("momentsAIC") = results.moments_AIC(),
                      Named("momentsBIC") = results.moments_BIC(),
                      Named("momentsHQ") = results.moments_HQ());
//                      Named("momentsCCIC_AIC") = results.moments_CCIC_AIC(),
//                      Named("momentsCCIC_BIC") = results.moments_CCIC_BIC(),
//                      Named("momentsCCIC_HQ") = results.moments_CCIC_HQ());                       
}

// [[Rcpp::export]]
List fmsc_test(mat x, colvec y, mat z1, mat z2, 
                      umat candidates){
  
  fmsc_chooseIV test(x, y, z1, z2, candidates);
  
  return List::create(Named("tau") = test.tau,
                      Named("Psi") = test.Psi,
                      Named("tau.outer") = test.tau_outer_est,
                      Named("Bias.mat") = test.Bias_mat,
                      Named("Estimates") = test.estimates,
                      Named("K") = test.K,
                      Named("Omega") = test.Omega,
                      Named("sq.bias.inner") = test.sqbias_inner,
                      Named("avar.inner") = test.avar_inner);
}