#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_epsilon_update(int n,
                             int p,
                             int m,
                             arma::vec y,
                             arma::mat x,
                             arma::mat z, 
                             int likelihood_indicator,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon,
                             arma::vec beta_old,
                             arma::vec eta_full,
                             arma::mat risk_sum){

double a_sigma2_epsilon_update = 0.50*n + 
                                 a_sigma2_epsilon;
  
arma::mat ident(m, m); ident.eye();

arma::vec mu = x*beta_old + 
               risk_sum*eta_full;
  
double b_sigma2_epsilon_update = 0.50*dot((y - mu), (y - mu)) + 
                                 b_sigma2_epsilon;

double sigma2_epsilon = 1.00/R::rgamma(a_sigma2_epsilon_update,
                                       (1.00/b_sigma2_epsilon_update));

return(sigma2_epsilon);

}



