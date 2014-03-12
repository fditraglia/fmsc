// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


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
//D_resid is a matrix containing the partial derivatives of the GMM residuals
//with respect to the parameter vector. 
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
    //the arguments of R::chisq are (q, df, lower.tail, log.p)
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
    colvec est_selected(mat est, colvec criterion){
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
  n_candidates = moment_sets.n_elem;
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
    linearGMM_msc candidate(X, y, Z_full.cols(moment_sets.col(i)));
    J(i) = candidate.Jstat();
    pJtest(i) = candidate.pJtest();
    AIC(i) = candidate.GMM_AIC();
    BIC(i) = candidate.GMM_BIC();
    HQ(i) = candidate.GMM_HQ();
    //Parameter Estimates
    estimates_1step.col(i) = candidate.est_1step();
    estimates_2step.col(i) = candidate.est_2step();
    //Hall & Peixe (2003) Criteria
    CCIC candidate_CCIC(X, Z_full.cols(moment_sets.col(i)));
    AIC_CCIC(i) = candidate_CCIC.AIC();
    BIC_CCIC(i) = candidate_CCIC.BIC();
    HQ_CCIC(i) = candidate_CCIC.HQ();
  }
}




class dgp {
  public:
    dgp(double, vec, double, double, mat, mat, int);
    colvec x, y, z2;
    mat z1;
  private: 
    int n_z1;
    mat u_e_z2;
};
//Class constructor
dgp::dgp(double b, vec p, double g, double r, mat V, 
                  mat Q, int n){
  RNGScope scope;
  n_z1 = Q.n_cols;
  z1 = trans(chol(Q) * reshape(colvec(rnorm(n * n_z1)), n_z1, n));
  u_e_z2 = trans(chol(V) * reshape(colvec(rnorm(3 * n)), 3, n));
  z2 = u_e_z2.col(2);
  x = z1 * p + g * z2 + u_e_z2.col(1);
  y = b * x + u_e_z2.col(0);
}


class fmsc {
  public:
    fmsc(const mat&, const colvec&, const mat&, const mat&, umat); 
    colvec est_full(){return(full.est();)}
    colvec est_valid(){return(valid.est());}
  private:
    tsls_fit valid, full;
    colvec tau;
    mat Psi, tau_outer_est, Bias_mat;
    int n_obs, n_z1, n_z2, n_z;
};
//Class constructor - need to use initialization list here
//This ensures that the tsls_fit constructor is called to 
//set up valid and full *before* we enter the body of the
//present constructor. 
fmsc::fmsc(const mat& x, const colvec& y, const mat& Z1, 
           const mat& Z2, umat candidates = zeros(1,1)): 
                valid(x, y, z1), full(x, y, join_rows(z1, z2)){
    n_z1 = z1.n_cols;
    n_z2 = z2.n_cols;
    n_z = n_z1 + n_z2;
    n_obs = y.n_elem;
    tau = z2.t() * valid.resid();
    Psi =  join_rows(-1 * z2.t() * valid.C , eye(n_z2, n_z2));
    tau_outer_est = tau * tau.t() - Psi * full.Omega_center() * Psi.t();
    Bias_mat = mat(n_z, n_z, fill::zeros);
    Bias_mat(span(n_z1, n_z - 1), span(n_z1, n_z - 1)) = tau_outer_est;
    if(all(vectorise(candidates) == 0)){
    }else{
    }
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
  dgp sims(b, p, g, r, V, Q, n);
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