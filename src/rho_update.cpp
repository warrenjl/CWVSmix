#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List rho_update(double rho_old,
                      int p,
                      int m,
                      arma::mat alpha_lambda,
                      arma::mat lambda_star,
                      double metrop_var_rho_trans,
                      int acctot_rho_trans){

double second = 0.00;
double first = 0.00;

/*Second*/
double rho_trans_old = log(rho_old/(1 - rho_old));

for(int j = 1; j < m; ++ j){
   for(int k = 0; k < p; ++ k){
      second = second +
               -lgamma(alpha_lambda(k, j) + rho_old*lambda_star(k, (j-1))) +
               (alpha_lambda(k, j) + rho_old*lambda_star(k, (j-1)) - 1)*log(lambda_star(k,j)) +
               -lambda_star(k,j);
      }
   }
second = second +
         rho_trans_old + 
         -2*log(1 + exp(rho_trans_old));
         
/*First*/
double rho_trans = R::rnorm(rho_trans_old, 
                            sqrt(metrop_var_rho_trans));
double rho = 1/(1 + exp(-rho_trans));

for(int j = 1; j < m; ++ j){
   for(int k = 0; k < p; ++ k){
      first = first +
              -lgamma(alpha_lambda(k, j) + rho*lambda_star(k, (j-1))) +
              (alpha_lambda(k, j) + rho*lambda_star(k, (j-1)) - 1)*log(lambda_star(k,j)) +
              -lambda_star(k,j);
      }
   }
first = first +
        rho_trans + 
        -2*log(1 + exp(rho_trans));
  
/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
     
  rho = rho_old;
  acc = 0;
     
  }
acctot_rho_trans = acctot_rho_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("rho") = rho,
                          Rcpp::Named("acctot_rho_trans") = acctot_rho_trans);

}
