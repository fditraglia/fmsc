// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

class tsls_fit {
  public:
    tsls_fit(const mat&, const colvec&, const mat&);
    colvec est(); 
    colvec resid();
  private:
    int n, k;
    double s_sq;
    colvec b, residuals;
    mat Z_copy, Qz, Rz, Xtilde, Qtilde, Rtilde, Rtilde_inv, D, C;
    mat Omega_textbook();
    mat Omega_robust();
    mat Omega_center();
    mat V_textbook();
    mat V_robust();
    mat V_center();
    colvec SE_textbook();
    colvec SE_robust();
    colvec SE_center();
};


//Class constructor
tsls_fit::tsls_fit(const mat& X, const colvec& y, const mat& Z){
  n = residuals.n_elem;
  k = Z.n_rows;
  Z_copy = Z; 
  qr_econ(Qz, Rz, Z);
  Xtilde = Qz.t() * X;
  qr_econ(Qtilde, Rtilde, Xtilde);
  b = solve(trimatu(Rtilde), Qtilde.t() * Qz.t() * y);
  residuals = y - X * b;
  D =  diagmat(pow(residuals, 2));
  s_sq = dot(residuals, residuals) / n;
  Rtilde_inv = solve(trimatu(Rtilde), 
            eye(Rtilde.n_rows, Rtilde.n_cols));
  //C = solve(trimatu(Rtilde), Qtilde.t() * Rtilde_inv.t());
}

mat tsls_fit::Omega_textbook(){
//Assumes heteroskedastic errors
  return(s_sq * Rz.t() * Rz / n);
}

mat tsls_fit::Omega_robust(){
//Heteroskedasticity-robust
  return(Z_copy.t() * D * Z_copy / n);
}

mat tsls_fit::Omega_center(){
//Heteroskedasticity-robust and centered
  mat e_outer = residuals * residuals.t();
  return(Z_copy.t() * (D / n -  e_outer / (n * n)) * Z_copy);
}

mat tsls_fit::V_textbook(){
//Covariance matrix estimator for sqrt(n) * (b_hat - b_true)
  return(n * s_sq * Rtilde_inv * Rtilde_inv.t());
}

mat tsls_fit::V_robust(){
//Covariance matrix estimator for sqrt(n) * (b_hat - b_true)
  return(n * n * C * Omega_robust() * C.t());
}

mat tsls_fit::V_center(){
//Covariance matrix estimator for sqrt(n) * (b_hat - b_true)
  return(n * n * C * Omega_center() * C.t());
}

colvec tsls_fit::SE_textbook(){
  return(sqrt(diagvec(V_textbook()) / n));
}

colvec tsls_fit::SE_robust(){
  return(sqrt(diagvec(V_robust()) / n));
}

colvec tsls_fit::SE_center(){
  return(sqrt(diagvec(V_center()) / n));
}

colvec tsls_fit::est(){
  return(b);
}

colvec tsls_fit::resid(){
  return(residuals);
}

// [[Rcpp::export]]
colvec tsls_est_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.est());
}

//// [[Rcpp::export]]
//colvec tsls_SE_cpp(mat X, colvec y, mat Z) {
//   tsls_fit results(X, y, Z);
//   return(results.est());
//}
