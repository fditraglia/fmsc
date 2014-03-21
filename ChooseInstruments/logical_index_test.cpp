// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat logical_subset(mat M, vec col_indices){
  uvec candidate_ind = find(col_indices);
  mat out = M.cols(candidate_ind);
  return(out);
  
}