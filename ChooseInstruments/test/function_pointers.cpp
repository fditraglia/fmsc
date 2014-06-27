//This file contains a simple example of function pointers in C++
//I learned about this topic from http://www.newty.de/fpt/intro.html
//and my first example is similar to one given there.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


//These will be passed as arguments to another function
double Plus     (double a, double b) { return a + b; }
double Minus    (double a, double b) { return a - b; }
double Multiply (double a, double b) { return a * b; }
double Divide   (double a, double b) { return a / b; }

//This function adds one to the result of applying a user-specified
//function that returns a double to two doubles a and b
double ResultPlusOne (double a, double b, 
                   double (*pt2Function)(double, double)){
  
  double result = pt2Function(a, b);
  return(result + 1);
  
}


//Now we can switch between different operations just by changing
//the function that we pass to ResultPlusOne

// [[Rcpp::export]]
double foo (double a, double b){
  double answer = ResultPlusOne(a, b, &Plus);
  return answer;
}

// [[Rcpp::export]]
double bar (double a, double b){
  double answer = ResultPlusOne(a, b, &Minus);
  return answer;
}

//Now how about a more interesting example. This is closer to what
//I'll actually use in my FMSC code

//This is the basic code for an "indicator vector"
colvec indicator (colvec b_valid, int which_index){
  colvec D_mu = zeros<colvec>(b_valid.n_elem);
  D_mu(which_index) = 1;
  return D_mu;
}

//From which particular indicators can be generated
colvec first_element (colvec b_valid){
  return indicator(b_valid, 0);
}
colvec second_element (colvec b_valid){
  return indicator(b_valid, 1);
}

//Here's an example that actually depends on b_valid
colvec squares (colvec b_valid){
  return pow(b_valid, 2);
}

//This function actually does the computations
double quadratic_form (colvec b_valid, mat Omega, 
                       colvec (*pt2Function)(colvec)){
  colvec D_mu = pt2Function(b_valid);
  return as_scalar(D_mu.t() * Omega * D_mu);
}

//And these functions simply pass arguments to quadratic_form
// [[Rcpp::export]]
double baz1 (colvec b_valid, mat Omega){
  return quadratic_form(b_valid, Omega, &first_element);
}
// [[Rcpp::export]]
double baz2 (colvec b_valid, mat Omega){
    return quadratic_form(b_valid, Omega, &second_element);  
}

// [[Rcpp::export]]
double baz3 (colvec b_valid, mat Omega){
    return quadratic_form(b_valid, Omega, &squares);  
}

//One more example: make sure we can pass a function pointer *through*
//one function and *into* another.
double EmptyShell (double a, double b, 
                   double (*pt2Function)(double, double)){
  double answer = ResultPlusOne(a, b, *pt2Function);
  return(answer);
}

// [[Rcpp::export]]
double PassThrough (double a, double b){
  return(EmptyShell(a, b, Plus));
}
