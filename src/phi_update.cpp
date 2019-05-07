#include "RcppArmadillo.h"
#include "CWMix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_update(double phi_old,
                      double a_phi,
                      double b_phi,
                      arma::vec eta,
                      double sigma2_eta,
                      Rcpp::List temporal_corr_info,
                      double beta_phi_old,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans){
  
int m = eta.size();

/*Second*/
Rcpp::List temporal_corr_info_old = temporal_corr_info;
arma::mat corr_inv_old = temporal_corr_info_old[0];
double log_deter_old = temporal_corr_info_old[1];
double phi_trans_old = log((phi_old - a_phi)/(b_phi - phi_old));

double second = -0.50*log_deter_old - 
                0.50*(1/sigma2_eta)*dot(eta, (corr_inv_old*eta)) + 
                phi_trans_old -
                2*log(1 + exp(phi_trans_old)) -
                (beta_phi_old - 1.00)*log(1 + exp(phi_trans_old));

/*First*/
double phi_trans = R::rnorm(phi_trans_old, 
                            sqrt(metrop_var_phi_trans));
double phi = (a_phi + b_phi*exp(phi_trans))/(1 + exp(phi_trans));
temporal_corr_info = temporal_corr_fun(m, 
                                       phi);
arma::mat corr_inv = temporal_corr_info[0];
double log_deter = temporal_corr_info[1];

double first = -0.50*log_deter - 
               0.50*(1/sigma2_eta)*dot(eta, (corr_inv*eta)) + 
               phi_trans -
               2*log(1 + exp(phi_trans)) -
               (beta_phi_old - 1.00)*log(1 + exp(phi_trans));

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
                 
  
