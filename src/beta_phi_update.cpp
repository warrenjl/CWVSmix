#include "RcppArmadillo.h"
#include "CWMix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double beta_phi_update(int q,
                       arma::vec phi,
                       double alpha_beta_phi,
                       double beta_beta_phi){

double alpha_beta_phi_update = q +
                               alpha_beta_phi;

double beta_beta_phi_update = sum(phi) + 
                              beta_beta_phi;

double beta_phi = R::rgamma(alpha_beta_phi_update,
                            (1.00/beta_beta_phi_update));

return(beta_phi);

}





