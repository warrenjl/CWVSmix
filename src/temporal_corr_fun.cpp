#include "RcppArmadillo.h"
#include "CWMix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List temporal_corr_fun(int p_z,
                             double phi){

double log_deter = 0.00; 
double sign = 0.00;     
arma::mat temporal_corr(p_z, p_z);
for(int j = 0; j < p_z; ++ j){
   for(int k = 0; k < p_z; ++ k){
      temporal_corr(j,k) = exp(-phi*abs(j - k));
      }
   }

arma::mat temporal_corr_inv = inv_sympd(temporal_corr);
log_det(log_deter, 
        sign, 
        temporal_corr);

return Rcpp::List::create(Rcpp::Named("temporal_corr_inv") = temporal_corr_inv,
                          Rcpp::Named("log_deter") = log_deter);

}

