// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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
    mat Omega_center() {return(Z_copy.t() * (D / n -  (residuals * residuals.t()) / (n * n)) * Z_copy);};
    mat V_textbook() {return(n * s_sq * Rtilde_inv * Rtilde_inv.t());}
    mat V_robust() {return(n * n * C * Omega_robust() * C.t());}
    mat V_center() {return(n * n * C * Omega_center() * C.t());}
  private:
    int n, k;
    double s_sq;
    colvec b, residuals;
    mat Z_copy, Qz, Rz, Qtilde, Rtilde, Rtilde_inv, D, C;
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