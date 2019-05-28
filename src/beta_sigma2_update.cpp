#include "RcppArmadillo.h"
#include "CWMix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double beta_sigma2_update(int q,
                          double beta_beta_sigma2,
                          arma::vec sigma2_eta){

double alpha_beta_sigma2_update = 3.00*q +
                                  2.00*beta_beta_sigma2;

double beta_beta_sigma2_update = sum(1.00/sigma2_eta) + 
                                 beta_beta_sigma2;

double beta_sigma2 = R::rgamma(alpha_beta_sigma2_update,
                               (1.00/beta_beta_sigma2_update));

return(beta_sigma2);

}





