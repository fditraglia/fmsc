//This file contains a simple example of function pointers in C++
//It is similar to an example from: 
//http://www.newty.de/fpt/intro.html

#include <Rcpp.h>
using namespace Rcpp;


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