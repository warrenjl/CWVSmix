#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List lambda_update(arma::mat lambda_star,
                         arma::mat lambda_old,
                         int ind,
                         int n,
                         int p,
                         int m,
                         arma::mat x,
                         arma::mat z,
                         arma::vec alpha_lambda,
                         arma::vec w,
                         arma::vec gamma,
                         arma::vec beta,
                         arma::vec eta_full,
                         arma::mat risk_sum,
                         arma::vec metrop_var_lambda,
                         arma::vec acctot_lambda){
  
arma::mat ident(m, m); ident.eye();
arma::mat lambda_star_old = lambda_star;
arma::mat lambda = lambda_old;

double second = 0.00;
double first = 0.00;

for(int j = 0; j < p; ++ j){
  
   arma::mat lambda_old = lambda;
   arma::mat risk_sum_old = risk_sum;
  
   /*Second*/
   if((ind > 0) & (ind < (m - 1))){
     second = -0.50*dot((gamma - x*beta - risk_sum_old*eta_full), w%(gamma - x*beta - risk_sum_old*eta_full)) +
              log(lambda_star_old(j, ind))*(alpha_lambda(j) + lambda_star_old(j, (ind - 1))) +
              -lambda_star_old(j, ind) +
              -lgamma(alpha_lambda(j) + lambda_star_old(j, (ind - 1))) +
              log(lambda_star_old(j, (ind + 1)))*(alpha_lambda(j) + lambda_star_old(j, ind)) +
              -lambda_star_old(j, (ind + 1)) +
              -lgamma(alpha_lambda(j) + lambda_star_old(j, ind));
     }
   
   if(ind == (m - 1)){
     second = -0.50*dot((gamma - x*beta - risk_sum_old*eta_full), w%(gamma - x*beta - risk_sum_old*eta_full)) +
              log(lambda_star_old(j, ind))*(alpha_lambda(j) + lambda_star_old(j, (ind - 1))) +
              -lambda_star_old(j, ind) +
              -lgamma(alpha_lambda(j) + lambda_star_old(j, (ind - 1)));
     }
   
   if(ind == 0){
     second = -0.50*dot((gamma - x*beta - risk_sum_old*eta_full), w%(gamma - x*beta - risk_sum_old*eta_full)) +
              log(lambda_star_old(j, ind))*alpha_lambda(j) +
              -lambda_star_old(j, ind) +
              log(lambda_star_old(j, (ind + 1)))*(alpha_lambda(j) + lambda_star_old(j, ind)) +
              -lambda_star_old(j, (ind + 1)) +
              -lgamma(alpha_lambda(j) + lambda_star_old(j, ind));
     }
  
   /*First*/
   lambda_star(j, ind) = exp(R::rnorm(log(lambda_star_old(j, ind)),
                                      sqrt(metrop_var_lambda(j))));

   for(int k = 0; k < p; ++ k){
      lambda(k, ind) = lambda_star(k, ind)/sum(lambda_star.col(ind));
      }
   risk_sum.col(ind) =  z.cols(p*ind, (p*(ind + 1) - 1))*lambda.col(ind);
   
   if((ind > 0) & (ind < (m - 1))){
     first = -0.50*dot((gamma - x*beta - risk_sum*eta_full), w%(gamma - x*beta - risk_sum*eta_full)) + 
             log(lambda_star(j, ind))*(alpha_lambda(j) + lambda_star(j, (ind - 1))) +
             -lambda_star(j, ind) +
             -lgamma(alpha_lambda(j) + lambda_star(j, (ind - 1))) +
             log(lambda_star(j, (ind + 1)))*(alpha_lambda(j) + lambda_star(j, ind)) +
             -lambda_star(j, (ind + 1)) +
             -lgamma(alpha_lambda(j) + lambda_star(j, ind));
     }
   
   if(ind == (m - 1)){
     first = -0.50*dot((gamma - x*beta - risk_sum*eta_full), w%(gamma - x*beta - risk_sum*eta_full)) +
             log(lambda_star(j, ind))*(alpha_lambda(j) + lambda_star(j, (ind - 1))) +
             -lambda_star(j, ind) +
             -lgamma(alpha_lambda(j) + lambda_star(j, (ind - 1)));
     }
   
   if(ind == 0){
     first = -0.50*dot((gamma - x*beta - risk_sum*eta_full), w%(gamma - x*beta - risk_sum*eta_full)) +
             log(lambda_star(j, ind))*alpha_lambda(j) +
             -lambda_star(j, ind) +
             log(lambda_star(j, (ind + 1)))*(alpha_lambda(j) + lambda_star(j, ind)) +
             -lambda_star(j, (ind + 1)) +
             -lgamma(alpha_lambda(j) + lambda_star(j, ind));
     }
               
   /*Decision*/
   double ratio = exp(first - second);   
   int acc = 1;
   if((ratio < R::runif(0.00, 1.00))){
     
     risk_sum = risk_sum_old;
     lambda_star(j, ind) = lambda_star_old(j, ind);
     lambda = lambda_old;
     acc = 0;
     
     }
   acctot_lambda(j) = acctot_lambda(j) + 
                      acc;
   
   }

return Rcpp::List::create(Rcpp::Named("risk_sum") = risk_sum.col(ind),
                          Rcpp::Named("lambda_star") = lambda_star.col(ind),
                          Rcpp::Named("lambda_vec") = lambda.col(ind),
                          Rcpp::Named("acctot_lambda") = acctot_lambda);

}
                 
  
