/*------------------------------------------------------------
Filename:        dgp_alt.cpp
Author:          Frank DiTraglia
First Version:   2013-23-11
This Version:    2013-23-11
--------------------------------------------------------------
Generates simulated data for the OLS versus IV FMSC example.
Unlike dgp.cpp, this function generates random numbers using
the Rcpp "sugar" function rnorm. This will allow me to get
the same random numbers from R as from C++, set the seed, etc.
------------------------------------------------------------*/
  
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List dgp_alt_cpp(double b, NumericVector PI_r, NumericMatrix Ve_r, 
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
    
    RNGScope scope; // also done by sourceCpp()
    
    //Number of instruments
    int n_z = Vz_r.ncol();
    
    //Initialize armadillo matrices corresponding to R input matrices
    //that reuse the original memory
    arma::mat Ve(Ve_r.begin(), 2, 2, false); //Error terms 
    arma::mat Vz(Vz_r.begin(), n_z, n_z, false); //Instruments
    arma::colvec p(PI_r.begin(), PI_r.size(), false); //Can't call it PI!
    
    //Generate errors
    arma::colvec stdnorm_errors = rnorm(n * 2);
    arma::mat errors = trans(arma::chol(Ve) * reshape(stdnorm_errors, 2, n));
    
    //Generate instruments
    arma::colvec stdnorm_z = rnorm(n * n_z);
    arma::mat z = trans(arma::chol(Vz) * reshape(stdnorm_z, n_z, n));

    //Generate endogenous variables: x and y
    arma::colvec x = z * p + errors.col(1); //Remember: zero indexing!
    arma::colvec y = b * x + errors.col(0); 
    
    return List::create(Named("x") = x, Named("y") = y, Named("z") = z);
                        
}