#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List rho_update(double rho_old,
                      int p,
                      int m,
                      double alpha_rho,
                      double beta_rho,
                      arma::mat lambda_star,
                      Rcpp::List temporal_corr_info,
                      double metrop_var_rho_trans,
                      int acctot_rho_trans){
  
/*Second*/
Rcpp::List temporal_corr_info_old = temporal_corr_info;
arma::mat temporal_corr_inv_old = temporal_corr_info_old[0];
double log_deter_old = temporal_corr_info_old[1];
double rho_trans_old = log(rho_old);

arma::vec lambda_star_vec(p*m); lambda_star_vec.fill(0.00);
for(int j = 0; j < m; ++ j){
   lambda_star_vec.subvec((p*j), (p*(j + 1) - 1)) = lambda_star.col(j);
   }

double second = -0.50*p*log_deter_old + 
                -0.50*dot(lambda_star_vec, (kron(temporal_corr_inv_old, eye(p, p))*lambda_star_vec)) + 
                alpha_rho*rho_trans_old +
                -beta_rho*exp(rho_trans_old);

/*First*/
double rho_trans = R::rnorm(rho_trans_old, 
                            sqrt(metrop_var_rho_trans));
double rho = exp(rho_trans);
temporal_corr_info = temporal_corr_fun(m, 
                                       rho);
arma::mat temporal_corr_inv = temporal_corr_info[0];
double log_deter = temporal_corr_info[1];

double first = -0.50*p*log_deter + 
               -0.50*dot(lambda_star_vec, (kron(temporal_corr_inv, eye(p, p))*lambda_star_vec)) + 
               alpha_rho*rho_trans +
               -beta_rho*exp(rho_trans);

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  
  rho = rho_old;
  temporal_corr_info = temporal_corr_info_old;
  acc = 0;
  
  }
acctot_rho_trans = acctot_rho_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("rho") = rho,
                          Rcpp::Named("acctot_rho_trans") = acctot_rho_trans,
                          Rcpp::Named("temporal_corr_info") = temporal_corr_info);

}
                 
  
