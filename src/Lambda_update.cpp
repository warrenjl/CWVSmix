#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List Lambda_update(arma::vec lambda_star,
                         arma::mat Lambda_old,
                         int ind,
                         int p,
                         int q,
                         int m,
                         arma::mat x,
                         arma::mat z,
                         arma::vec alpha_Lambda,
                         arma::vec w,
                         arma::vec gamma,
                         arma::vec beta,
                         arma::vec eta_full,
                         arma::vec metrop_var_Lambda,
                         arma::vec acctot_Lambda){
  
arma::mat ident(m, m); ident.eye();
arma::vec lambda_star_old = lambda_star;
arma::mat Lambda = Lambda_old;

for(int j = ind; j < p; ++ j){
  
   arma::mat Lambda_old = Lambda;
  
   /*Second*/
   double second = -0.50*dot((gamma - x*beta - z*((kron(ident, Lambda_old))*eta_full)), w%(gamma - x*beta - z*((kron(ident, Lambda_old))*eta_full))) +
                   alpha_Lambda(j)*log(lambda_star_old(j)) +
                   -lambda_star_old(j);
  
   /*First*/
   lambda_star(j) = exp(R::rnorm(log(lambda_star_old(j)),
                                 sqrt(metrop_var_Lambda(j))));

   for(int k = ind; k < p; ++ k){
      Lambda(k, ind) = lambda_star(k)/sum(lambda_star.subvec(ind, (p - 1)));
      }
   
   double first = -0.50*dot((gamma - x*beta - z*((kron(ident, Lambda))*eta_full)), w%(gamma - x*beta - z*((kron(ident, Lambda))*eta_full))) + 
                  alpha_Lambda(j)*log(lambda_star(j)) +
                  -lambda_star(j);
               
   /*Decision*/
   double ratio = exp(first - second);   
   int acc = 1;
   if((ratio < R::runif(0.00, 1.00))){
     
     lambda_star(j) = lambda_star_old(j);
     Lambda = Lambda_old;
     acc = 0;
     
     }
   acctot_Lambda(j) = acctot_Lambda(j) + 
                      acc;
   
   }

return Rcpp::List::create(Rcpp::Named("lambda_star") = lambda_star,
                          Rcpp::Named("Lambda_vec") = Lambda.col(ind),
                          Rcpp::Named("acctot_Lambda") = acctot_Lambda);

}
                 
  
