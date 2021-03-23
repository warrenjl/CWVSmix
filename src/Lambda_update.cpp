#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List lambda_update(arma::mat lambda_star,
                         arma::mat lambda_old,
                         int ind,
                         int p,
                         int m,
                         arma::mat x,
                         arma::mat z,
                         arma::vec w,
                         arma::vec gamma,
                         arma::vec beta,
                         arma::vec eta_full,
                         arma::mat risk_sum,
                         arma::mat temporal_corr_inv_old,
                         double metrop_var_lambda_trans,
                         int acctot_lambda){
  
arma::mat lambda_star_old = lambda_star;
arma::mat lambda = lambda_old;
  
arma::mat corr_inv_old = kron(temporal_corr_inv_old, eye(p, p));
  
arma::vec lambda_star_vec_old(p*m); lambda_star_vec_old.fill(0.00);
for(int j = 0; j < m; ++ j){
   lambda_star_vec_old.subvec((p*j), (p*(j + 1) - 1)) = lambda_star_old.col(j);
   }
arma::vec lambda_star_vec = lambda_star_vec_old;
  
arma::mat risk_sum_old = risk_sum;
  
double second = 0.00;
double first = 0.00;
  
/*Second*/
second = -0.50*dot((gamma - x*beta - risk_sum_old*eta_full), w%(gamma - x*beta - risk_sum_old*eta_full)) +
         -0.50*dot(lambda_star_vec_old, (corr_inv_old*lambda_star_vec_old));
  
/*First*/
arma::vec lambda_star_max(p); lambda_star_max.fill(0.00);
  
for(int j = 0; j < p; ++ j){
      
   lambda_star(j, ind) = R::rnorm(lambda_star_old(j, ind),
                                  sqrt(metrop_var_lambda_trans));
   lambda_star_vec(p*ind + j) = lambda_star(j, ind);
   lambda_star_max(j) = (lambda_star(j, ind) > 0.00)*lambda_star(j, ind);
      
   }
  
if(sum(lambda_star_max) > 0.00){
  for(int j = 0; j < p; ++ j){
     lambda(j, ind) = lambda_star_max(j)/sum(lambda_star_max);
     }
  }
  
if(sum(lambda_star_max) == 0.00){
  lambda.col(ind).fill(1.00/p);
  }
  
risk_sum.col(ind) =  z.cols(p*ind, (p*(ind + 1) - 1))*lambda.col(ind);
  
first = -0.50*dot((gamma - x*beta - risk_sum*eta_full), w%(gamma - x*beta - risk_sum*eta_full)) + 
        -0.50*dot(lambda_star_vec, (corr_inv_old*lambda_star_vec));
    
/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if((ratio < R::runif(0.00, 1.00))){
      
  risk_sum = risk_sum_old;
  lambda_star.col(ind) = lambda_star_old.col(ind);
  lambda = lambda_old;
  acc = 0;
      
  }
acctot_lambda = acctot_lambda + 
                acc;
    
return Rcpp::List::create(Rcpp::Named("risk_sum") = risk_sum.col(ind),
                          Rcpp::Named("lambda_star") = lambda_star.col(ind),
                          Rcpp::Named("lambda") = lambda.col(ind),
                          Rcpp::Named("acctot_lambda") = acctot_lambda);
    
}       