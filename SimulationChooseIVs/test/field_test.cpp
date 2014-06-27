// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat field_test(mat A, mat B, mat C, int which){
  
  field<mat> my_field(3);
  my_field(0) = A;
  my_field(1) = B;
  my_field(2) = C;
  mat out = my_field(which);
  return(out);
  
}