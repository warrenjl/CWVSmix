#include "RcppArmadillo.h"
#include "CWMix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_eta_update(int p,
                         arma::mat z,
                         arma::vec eta,
                         arma::mat corr_inv,
                         double beta_sigma2_old){

int m = z.n_cols/p;
double alpha_sigma2_eta_update = 0.50*m + 
                                 3.00;

double beta_sigma2_eta_update = 0.50*dot(eta, ((corr_inv)*eta)) + 
                                beta_sigma2_old;

double sigma2_eta = 1.00/R::rgamma(alpha_sigma2_eta_update,
                                   (1.00/beta_sigma2_eta_update));

return(sigma2_eta);

}





