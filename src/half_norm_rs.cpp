// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================

#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double half_norm_rs(double a, double b)
{
   double   x;
   x = fabs(norm_rand());
   while( (x<a) || (x>b) ) x = fabs(norm_rand());
   return x;
}

