/*
* https://github.com/slwu89/MCMC/blob/master/tmvrnormGibbs.cpp
* rtnorm_gibbs returns a sample of size n from the specificed truncated Gaussian distribution.
* n: integer number of samples to take
* mu: mean of distribution
* sigma: standard deviation of distribution
* a: lower truncation bound
* b: upper truncation bound
*/

#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec rtnorm_gibbs(int n, 
                       double mu, 
                       double sigma, 
                       double a, 
                       double b){
  
  //sample from uniform distribution on unit interval
  NumericVector F = runif(n);
  
  //Phi(a) and Phi(b)
  double Fa = R::pnorm(a,mu,sigma,1,0);
  double Fb = R::pnorm(b,mu,sigma,1,0);
  
  NumericVector F_out(F.length());
  for(int i=0; i < F.length(); i++){
    double p_i = F[i] * (Fb - Fa) + Fa;
    F_out[i] = R::qnorm(p_i,0.0,1.0,1,0);
  }
  
  NumericVector out(F.length());
  for(int i=0; i < out.length(); i++){
    out[i] = mu + sigma * F_out[i];
  }
  
  return(out);
}
