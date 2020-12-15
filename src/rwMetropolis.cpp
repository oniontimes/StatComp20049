#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

double f(double x) {
  return exp(-abs(x));
}

//' @title Rcpp function
//' @name rwMetropolis
//' @description an Rcpp function for Exercise 9.4 (page 277, Statistical Computing with R)
//' @param sigma standard deviation
//' @param x0 original value
//' @param N length
//' @return generated random numbers \code{x}
//' @examples
//' \dontrun{
//' sd=2
//' x0=25
//' N = 2000
//' rwMetropolis(sigma,x0,N)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwMetropolis (double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0; 
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (f(y[0]) / f(x[i-1]))){
      x[i] = y[0];
    }
    else { 
      x[i] = x[i-1]; 
    }
  }
  return(x);
} 