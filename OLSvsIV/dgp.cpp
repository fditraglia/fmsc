/*------------------------------------------------------------
Filename:        dgp.cpp
Author:          Frank DiTraglia
First Version:   2013-22-11
This Version:    2013-22-11
--------------------------------------------------------------
Generates simulated data for the OLS versus IV FMSC example.
------------------------------------------------------------*/
  
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List dgp_cpp(double b, NumericVector PI_r, NumericMatrix Ve_r, 
              NumericMatrix Vz_r, int n){
    /*------------------------------------------------------------------
      Arguments:
       b       coefficient on x in the second stage
       PI      vector of coefficients on z in the second stage
       V.e     variance-covariance matrix of the errors epsilon and v
       V.z     variance-covariance matrix of the instruments z
       n       sample size
      
      Returns: list of matrices x, y and z containing simulated dataset
    ------------------------------------------------------------------*/ 
    
    //Number of instruments
    int n_z = Vz_r.ncol();
    
    //Initialize armadillo matrices corresponding to R input matrices
    arma::mat Ve(Ve_r.begin(), 2, 2, false); //Error terms 
    arma::mat Vz(Vz_r.begin(), n_z, n_z, false); //Instruments
    arma::colvec p(PI_r.begin(), PI_r.size(), false); //Can't call it PI!
    
    //Generate exogenous RVs: errors and instruments
    arma::mat errors = trans(arma::chol(Ve) * arma::randn(2, n));
    arma::colvec e = errors.col(0); //Remember: zero indexing!
    arma::colvec v = errors.col(1);
    arma::mat z = trans(arma::chol(Vz) * arma::randn(n_z, n));
    
    //Generate endogenous variables: x and y
    arma::colvec x = z * p + v;
    arma::colvec y = b * x + e;
    
    return List::create(Named("x") = x, Named("y") = y, Named("z") = z);
                        
}