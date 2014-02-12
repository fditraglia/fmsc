// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

class LinearGMM {
  public:
    LinearGMM(mat&, colvec&, mat&);
    colvec tsls_est(); //Return tsls estimate
    colvec twostep_est(); //Return efficient two-step GMM estimate
  
  private:
  
}


//Class constructor
LinearGMM::LinearGMM(const mat& X, const colvec& y, const mat& Z){
  
  
}

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
int timesTwo(int x) {
   return x * 2;
}
