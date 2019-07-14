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
arma::vec A11_diag(m*q); A11_diag.fill(0.00);
for(int j = 0; j < m; ++ j){
    
   delta_diag.subvec((j*q), (q*(j + 1) - 1)) = trans(delta_old.row(j));
   w1_full.subvec((j*q), (q*(j + 1) - 1)) = trans(w1_old.row(j));
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
double log_pi = 0.00;
arma::vec probs(2);

for(int j = 0; j < q; ++ j){

   for(int k = 0; k < m; ++ k){
     
      int counter = q*k + j;

      pieces.fill(0.00);
      for(int l = 0; l < 2; ++ l){
       
         delta(k,j) = l;
         delta_diag(counter) = l;
         eta_full = (delta_diag%A11_diag%w1_full);
         
         if(l == 0){
       
           //Middle
           if((j > 0) & (j < (q - 1))){
             
             arma::vec temp = trans(delta.submat(k, 0, k, (j - 1)));
             log_pi = log(1.00 - pi(k,j)*prod(1.00 - temp)); 
             
             for(int t = (j + 1); t < q; ++ t){
               
                arma::vec temp = trans(delta.submat(k, 0, k, (t - 1)));
                log_pi = log_pi +
                         log(pow((pi(k,t)*prod(1.00 - temp)), delta(k,t))*pow((1 - pi(k,t)*prod(1.00 - temp)), (1.00 - delta(k,t)))); 
               
                }
             
             }
      
           //Start, q > 1
           if((j == 0) & (q > 1)){
             
             log_pi = log(1.00 - pi(k,j));
             for(int t = (j + 1); t < q; ++ t){
               
                arma::vec temp = trans(delta.submat(k, 0, k, (t - 1)));
               
                log_pi = log_pi +
                         log(pow((pi(k,t)*prod(1.00 - temp)), delta(k,t))*pow((1.00 - pi(k,t)*prod(1.00 - temp)), (1.00 - delta(k,t)))); 
                
                }
             
             }
           
           //q=1
           if(q == 1){
             log_pi = log(1.00 - pi(k,j));
             }
      
           //End, q > 1
           if((j == (q - 1)) & (q > 1)){
             
             arma::vec temp = trans(delta.submat(k, 0, k, (j - 1)));
             log_pi = log(1.00 - prod(1.00 - temp)*pi(k,j));
             
             }
           
           }
         
         
         
         if(l == 1){
           
           //Middle
           if((j > 0) & (j < (q - 1))){
             
             arma::vec temp1 = trans(delta.submat(k, 0, k, (j - 1)));
             arma::vec temp2 = trans(delta.submat(k, (j + 1), k, (q - 1)));
             log_pi = log(prod(1.00 - temp1)*pi(k,j)*prod(1.00 - temp2));
             
             }
           
           //Start, q > 1
           if((j == 0) & (q > 1)){
             
             arma::vec temp = trans(delta.submat(k, (j + 1), k, (q - 1)));
             log_pi = log(pi(k,j)*prod(1.00 - temp));
             
             }
           
           //q=1
           if(q == 1){
             log_pi = log(pi(k,j));
             }
           
           //End, q > 1
           if((j == (q - 1)) & (q > 1)){
             
             arma::vec temp = trans(delta.submat(k, 0, k, (j - 1)));
             log_pi = log(prod(1.00 - temp)*pi(k,j));
             
             }
           
           }
   
         pieces(l) = -0.50*dot((gamma - x*beta - z*((kron(ident, Lambda))*eta_full)), w%(gamma - x*beta - z*((kron(ident, Lambda))*eta_full))) +
                     log_pi;
      
         }

      probs.fill(0.00);
      for(int l = 0; l < 2; ++ l){
     
         probs(l) = 1.00/(sum(exp(pieces - pieces(l))));
  
         if(arma::is_finite(probs(l)) == 0.00){
           probs(l) = 0.00;  /*Computational Correction*/
           }
      
         }
      
      delta_diag(counter) = as<double>(Rcpp::rbinom(1,
                                                    1,
                                                    probs(1)));
      delta(k,j) = delta_diag(counter);
      
      }
  
   }

eta_full = (delta_diag%A11_diag%w1_full);
  
return Rcpp::List::create(Rcpp::Named("delta") = delta,
                          Rcpp::Named("eta_full") = eta_full);

}



