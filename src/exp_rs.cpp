// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================

#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double exp_rs(double a, double b)
{
  double  z, u, rate;

//  Rprintf("in exp_rs");
  rate = 1/a;
//1/a

   // Generate a proposal on (0, b-a)
   z = R::rexp(rate);
   while(z > (b-a)) z = R::rexp(rate);
   u = R::runif(0.0, 1.0);

   while( log(u) > (-0.5*z*z))
   {
      z = R::rexp(rate);
      while(z > (b-a)) z = R::rexp(rate);
      u = R::runif(0.0,1.0);
   }
   return(z+a);
}





