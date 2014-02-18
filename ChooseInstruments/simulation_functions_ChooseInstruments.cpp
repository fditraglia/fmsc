// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

class tsls_fit {
  public:
    tsls_fit(const mat&, const colvec&, const mat&);
    colvec estimate(); //Return tsls estimate
  private:
    int n, k;
    colvec b, residuals;
    mat Z_copy, Qz, Rz, Xtilde, Qtilde, Rtilde, C;
    mat Omega();
    mat Omega_robust();
    mat Omega_center();
    //mat V();
    //mat V_robust();
    //mat V_center();
};


//Class constructor
tsls_fit::tsls_fit(const mat& X, const colvec& y, const mat& Z){
  qr_econ(Qz, Rz, Z);
  Xtilde = Qz.t() * X;
  qr_econ(Qtilde, Rtilde, Xtilde);
  b = solve(trimatu(Rtilde), Qtilde.t() * Qz.t() * y);
  residuals = y - b * X;
  n = residuals.n_elem;
  k = Z.n_rows;
  Z_copy = Z; 
  C = solve(trimatu(Rtilde), Qtilde * 
              solve(trimatl(Rtilde.t()), eye(k, k)));
}

mat tsls_fit::Omega(){
//Assumes heteroskedastic errors
  double s_sq = dot(residuals, residuals) / n;
  return(s_sq * Rz.t() * Rz / n);
}

mat tsls_fit::Omega_robust(){
//Heteroskedasticity-robust
  mat D = diagmat(pow(residuals, 2));
  return(Z_copy.t() * D * Z_copy / n);
}

mat tsls_fit::Omega_center(){
//Heteroskedasticity-robust and centered
  mat e_outer = residuals * residuals.t();
  mat D = diagmat(pow(residuals, 2));
  return(Z_copy.t() * (D / n -  e_outer / (n * n)) * Z_copy);
}

colvec tsls_fit::estimate(){
  return(b);
}

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
colvec tsls_cpp(mat X, colvec y, mat Z) {
   tsls_fit results(X, y, Z);
   return(results.estimate());
}
