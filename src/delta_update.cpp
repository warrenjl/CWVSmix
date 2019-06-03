#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List delta_update(arma::mat delta_old,
                        int p,
                        int q,
                        arma::vec y,
                        arma::mat x,
                        arma::mat z,
                        arma::vec w,
                        arma::vec gamma,
                        arma::vec beta,
                        arma::mat Lambda,
                        arma::mat w1_old,
                        arma::mat w2_old,
                        arma::vec A11_old,
                        arma::vec A22_old,
                        arma::vec A21_old){
  
int m = z.n_cols/p;
arma::mat ident(m, m); ident.eye();

arma::vec delta_diag(m*q); delta_diag.fill(0.00); 
arma::vec w1_full(m*q); w1_full.fill(0.00);
arma::vec w2_full(m*q); w2_full.fill(0.00);
arma::vec A11_diag(m*q); A11_diag.fill(0.00);
for(int j = 0; j < m; ++ j){
    
   delta_diag.subvec((j*q), (q*(j + 1) - 1)) = trans(delta_old.row(j));
   w1_full.subvec((j*q), (q*(j + 1) - 1)) = trans(w1_old.row(j));
   w2_full.subvec((j*q), (q*(j + 1) - 1)) = trans(w2_old.row(j));
   A11_diag.subvec((j*q), (q*(j + 1) - 1)) = A11_old;
    
   }
arma::vec eta_full = (delta_diag%A11_diag%w1_full);
arma::mat delta = delta_old;

arma::mat pi(m, q); pi.fill(0.00);
for(int j = 0; j < q; ++ j){
  
   arma::vec alpha = A21_old(j)*w1_old.col(j) +
                     A22_old(j)*w2_old.col(j);
  
   for(int k = 0; k < m; ++ k){
      pi(k,j) = R::pnorm(alpha(k),
                         0.00,
                         1.00,
                         true,
                         false);
      }
   
   }

arma::vec pieces(2); pieces.fill(0.00);
arma::vec log_pi(2); log_pi.fill(0.00);
arma::vec probs(2);

int counter = (m*q) - 1;
for(int j = (q - 1); j > -1; -- j){

   for(int k = (m - 1); k > -1; -- k){

      pieces.fill(0.00);
      
      //Middle
      if((j > 0) & (j < (q - 1))){
        
        log_pi(0) = log(1.00 - (sum(delta.col(j - 1)) > 0.00)*pi(k,j)) +
                    sum((1.00 - delta.col(j + 1))%log(1.00 - (sum(delta.col(j)) > 0.00)*pi.col(j + 1))) +
                    sum(delta.col(j + 1))*log(sum(delta.col(j)) > 0.00);
        
        log_pi(1) = log((sum(delta.col(j - 1)) > 0.00)*pi(k,j)) + 
                    sum((1.00 - delta.col(j + 1))%log(1.00 - pi.col(j + 1)));
        
        }
      
      //Start, q > 1
      if((j == 0) & (q > 1)){
        
        log_pi(0) = log(1.00 - pi(k,j)) +
                    sum((1.00 - delta.col(j + 1))%log(1.00 - (sum(delta.col(j)) > 0.00)*pi.col(j + 1))) +
                    sum(delta.col(j + 1))*log(sum(delta.col(j)) > 0.00);
        
        log_pi(1) = log(pi(k,j)) + 
                    sum((1.00 - delta.col(j + 1))%log(1.00 - pi.col(j + 1)));
        
        }
      
      //q=1
      if(q == 1){
        
        log_pi(0) = log(1.00 - pi(k,j));
        
        log_pi(1) = log(pi(k,j));
        
        }
      
      //End, q > 1
      if((j == (q - 1)) & (q > 1)){
        
        log_pi(0) = log(1.00 - (sum(delta.col(j - 1)) > 0.00)*pi(k,j));
        
        log_pi(1) = log((sum(delta.col(j - 1)) > 0.00)*pi(k,j));
        
        }
   
      for(int l = 0; l < 2; ++ l){
     
         delta_diag(counter) = l;
         eta_full = (delta_diag%A11_diag%w1_full);
         pieces(l) = -0.50*dot((gamma - x*beta - z*((kron(ident, Lambda))*eta_full)), w%(gamma - x*beta - z*((kron(ident, Lambda))*eta_full))) +
                     log_pi(l);
      
         }

      probs.fill(0.00);

      for(int l = 0; l < 2; ++l){
     
         probs(l) = 1.00/(sum(exp(pieces - pieces(l))));
  
         if(arma::is_finite(probs(l)) == 0.00){
           probs(l) = 0.00;  /*Computational Correction*/
           }
      
         }

      delta_diag(counter) = as<double>(Rcpp::rbinom(1,
                                                    1,
                                                    probs(1)));
      delta(k,j) = delta_diag(counter);
      -- counter;
      
      }
  
   }

eta_full = (delta_diag%A11_diag%w1_full);
  
return Rcpp::List::create(Rcpp::Named("delta") = delta,
                          Rcpp::Named("eta_full") = eta_full);

}



