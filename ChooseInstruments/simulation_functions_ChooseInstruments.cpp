// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

class LinearGMM {
  public:
    LinearGMM(const mat&, const colvec&, const mat&);
    colvec tsls_est(); //Return tsls estimate
    //colvec twostep_est(); //Return efficient two-step GMM estimate
    
  private:
    int n, k;
    colvec b_tsls, resid_tsls;
    mat Z_copy, Qz, Rz, Xtilde, Qtilde, Rtilde, C;
    mat Omega_tsls();
    mat Omega_tsls_robust();
    mat Omega_tsls_center();
};


//Class constructor
LinearGMM::LinearGMM(const mat& X, const colvec& y, const mat& Z){
  qr_econ(Qz, Rz, Z);
  Xtilde = Qz.t() * X;
  qr_econ(Qtilde, Rtilde, Xtilde);
  b_tsls = solve(trimatu(Rtilde), Qtilde.t() * Qz.t() * y);
  resid_tsls = y - b_tsls * X;
  n = resid_tsls.n_elem;
  k = Z.n_rows;
  Z_copy = Z; 
  C = solve(trimatu(Rtilde), Qtilde * 
              solve(trimatl(Rtilde.t()), eye(k, k)));
}

mat LinearGMM::Omega_tsls(){
//Assumes heteroskedastic errors
  double s_sq = dot(resid_tsls, resid_tsls) / n;
  return(s_sq * Rz.t() * Rz / n);
}

mat LinearGMM::Omega_tsls_robust(){
//Heteroskedasticity-robust
  mat D = diagmat(pow(resid_tsls, 2));
  return(Z_copy.t() * D * Z_copy / n);
}

mat LinearGMM::Omega_tsls_center(){
//Heteroskedasticity-robust and centered
  mat e_outer = resid_tsls * resid_tsls.t();
  mat D = diagmat(pow(resid_tsls, 2));
  return(Z_copy.t() * (D / n -  e_outer / (n * n)) * Z_copy);
}

colvec LinearGMM::tsls_est(){
  return(b_tsls);
}

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
colvec tsls_cpp(mat X, colvec y, mat Z) {
   LinearGMM results(X, y, Z);
   return(results.tsls_est());
}
