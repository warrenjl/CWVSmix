#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_update(double phi_old,
                      double alpha_phi,
                      double beta_phi,
                      arma::vec w,
                      Rcpp::List temporal_corr_info,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans){
  
int m = w.size();

/*Second*/
Rcpp::List temporal_corr_info_old = temporal_corr_info;
arma::mat corr_inv_old = temporal_corr_info_old[0];
double log_deter_old = temporal_corr_info_old[1];
double phi_trans_old = log(phi_old);

double second = -0.50*log_deter_old - 
                0.50*dot(w, (corr_inv_old*w)) + 
                alpha_phi*phi_trans_old -
                beta_phi*exp(phi_trans_old);

/*First*/
double phi_trans = R::rnorm(phi_trans_old, 
                            sqrt(metrop_var_phi_trans));
double phi = exp(phi_trans);
temporal_corr_info = temporal_corr_fun(m, 
                                       phi);
arma::mat corr_inv = temporal_corr_info[0];
double log_deter = temporal_corr_info[1];

double first = -0.50*log_deter - 
               0.50*dot(w, (corr_inv*w)) + 
               alpha_phi*phi_trans -
               beta_phi*exp(phi_trans);

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  
  phi = phi_old;
  temporal_corr_info = temporal_corr_info_old;
  acc = 0;
  
  }
acctot_phi_trans = acctot_phi_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("phi") = phi,
                          Rcpp::Named("acctot_phi_trans") = acctot_phi_trans,
                          Rcpp::Named("temporal_corr_info") = temporal_corr_info);

}
                 
  
