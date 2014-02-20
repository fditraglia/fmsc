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
    fmsc(colvec, colvec, mat, mat);  
  private:
    tsls_fit valid, full;
    colvec tau;
    mat Psi, Omega, tau_outer_est, Bias_mat;
    int n, q1, q2, q;
};
//Class constructor - need to use initialization list here
//This ensures that the tsls_fit constructor is called to 
//set up valid and full *before* we enter the body of the
//present constructor. 
fmsc::fmsc(colvec x, colvec y, mat z1, mat z2):
  valid(x, y, z1), full(x, y, join_rows(z1, z2)){
    q1 = z1.n_cols;
    q2 = z2.n_cols;
    n = y.n_elem;
    Omega = full.Omega_center();
    tau = z2.t() * valid.resid();
    Psi =  join_rows(-1 * z2.t() * valid.C , eye(q2, q2));
    tau_outer_est = tau * tau.t() - Psi * Omega * Psi.t();
    Bias_mat = mat(q, q, fill::zeros);
    Bias_mat(span(q1, q - 1), span(q1, q - 1)) = tau_outer_est;
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