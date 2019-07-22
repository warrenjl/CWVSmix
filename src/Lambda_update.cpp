#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List Lambda_update(arma::vec stick,
                         arma::mat Lambda_old,
                         int ind,
                         int p,
                         int q,
                         int m,
                         arma::mat x,
                         arma::mat z,
                         arma::vec w,
                         arma::vec gamma,
                         arma::vec beta,
                         arma::vec eta_full,
                         arma::vec metrop_scale_Lambda,
                         arma::vec acctot_Lambda){
  
arma::mat ident(m, m); ident.eye();
arma::vec stick_old = stick;
arma::mat Lambda = Lambda_old;

for(int j = ind; j < (p - 1); ++ j){
  
   arma::mat Lambda_old = Lambda;
  
   /*Second*/
   double second = -0.50*dot((gamma - x*beta - z*((kron(ident, Lambda_old))*eta_full)), w%(gamma - x*beta - z*((kron(ident, Lambda_old))*eta_full)));
  
   /*First*/
   stick(j) = R::runif(stick_old(j) - metrop_scale_Lambda(j),
                       stick_old(j) + metrop_scale_Lambda(j));
     
   if(stick(j) < 0){
     stick(j) = abs(stick(j));
     }
   
   if(stick(j) > 1){
     stick(j) = 2 - stick(j);
     }

   for(int k = j; k < p; ++ k){
     
      if(k > ind){
        Lambda(k, ind) = stick(k)*prod(1 - stick.subvec(ind, (k - 1)));
        }
   
      if(k == ind){
        Lambda(k, ind) = stick(k);
        }
      
      }
   
   double first = -0.50*dot((gamma - x*beta - z*((kron(ident, Lambda))*eta_full)), w%(gamma - x*beta - z*((kron(ident, Lambda))*eta_full)));
               
   /*Decision*/
   double ratio = exp(first - second);   
   int acc = 1;
   if((ratio < R::runif(0.00, 1.00))){
     
     stick(j) = stick_old(j);
     Lambda = Lambda_old;
     acc = 0;
     
     }
   acctot_Lambda(j) = acctot_Lambda(j) + 
                      acc;
   
   }

return Rcpp::List::create(Rcpp::Named("stick") = stick,
                          Rcpp::Named("Lambda_vec") = Lambda.col(ind),
                          Rcpp::Named("acctot_Lambda") = acctot_Lambda);

}
                 
  
