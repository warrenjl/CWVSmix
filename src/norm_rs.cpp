// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================

#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double norm_rs(double a, double b)
{
   double  x;
   x = Rf_rnorm(0.0, 1.0);
   while( (x < a) || (x > b) ) x = norm_rand();
   return x;
}

