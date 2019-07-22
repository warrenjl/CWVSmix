#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List delta_update(arma::mat delta_old,
                        int p,
                        int q,
                        int m,
                        arma::vec y,
                        arma::mat x,
                        arma::mat z,
                        arma::vec w,
                        arma::vec gamma,
                        arma::vec beta,
                        arma::mat Lambda,
                        arma::vec w1_old,
                        arma::vec w2_old,
                        double A11_old,
                        double A22_old,
                        double A21_old){
  
arma::mat ident(m, m); ident.eye();

arma::vec eta_full(m*q); eta_full.fill(0.00); 
for(int j = 0; j < m; ++ j){
   eta_full.subvec((j*q), (q*(j + 1) - 1)) = A11_old*w1_old(j)*trans(delta_old.row(j));
   }
arma::mat delta = delta_old;

arma::vec pi(m); pi.fill(0.00);

arma::vec alpha = A21_old*w1_old +
                  A22_old*w2_old;

for(int k = 0; k < m; ++ k){
   pi(k) = R::pnorm(alpha(k),
                    0.00,
                    1.00,
                    true,
                    false);
   }

arma::vec pieces(2); pieces.fill(0.00);
double log_pi = 0.00;
arma::vec probs(2);

for(int j = 0; j < q; ++ j){

   for(int k = 0; k < m; ++ k){
     
      pieces.fill(0.00);
      for(int l = 0; l < 2; ++ l){
       
         delta(k, j) = l;
         eta_full((k*q) + j) = A11_old*w1_old(k)*delta(k, j);
         
         if(l == 0){
       
           //Middle
           if((j > 0) & (j < (q - 1))){
             
             arma::vec temp = trans(delta.submat(k, 0, k, (j - 1)));
             log_pi = log(1.00 - pi(k)*prod(1.00 - temp)); 
             
             for(int t = (j + 1); t < q; ++ t){
               
                arma::vec temp = trans(delta.submat(k, 0, k, (t - 1)));
                log_pi = log_pi +
                         log(pow((pi(k)*prod(1.00 - temp)), delta(k, t))*pow((1 - pi(k)*prod(1.00 - temp)), (1.00 - delta(k, t)))); 
               
                }
             
             }
      
           //Start, q > 1
           if((j == 0) & (q > 1)){
             
             log_pi = log(1.00 - pi(k));
             for(int t = (j + 1); t < q; ++ t){
               
                arma::vec temp = trans(delta.submat(k, 0, k, (t - 1)));
               
                log_pi = log_pi +
                         log(pow((pi(k)*prod(1.00 - temp)), delta(k, t))*pow((1.00 - pi(k)*prod(1.00 - temp)), (1.00 - delta(k, t)))); 
                
                }
             
             }
           
           //q=1
           if(q == 1){
             log_pi = log(1.00 - pi(k));
             }
      
           //End, q > 1
           if((j == (q - 1)) & (q > 1)){
             
             arma::vec temp = trans(delta.submat(k, 0, k, (j - 1)));
             log_pi = log(1.00 - prod(1.00 - temp)*pi(k));
             
             }
           
           }
         
         
         
         if(l == 1){
           
           //Middle
           if((j > 0) & (j < (q - 1))){
             
             arma::vec temp1 = trans(delta.submat(k, 0, k, (j - 1)));
             arma::vec temp2 = trans(delta.submat(k, (j + 1), k, (q - 1)));
             log_pi = log(prod(1.00 - temp1)*pi(k)*prod(1.00 - temp2));
             
             }
           
           //Start, q > 1
           if((j == 0) & (q > 1)){
             
             arma::vec temp = trans(delta.submat(k, (j + 1), k, (q - 1)));
             log_pi = log(pi(k)*prod(1.00 - temp));
             
             }
           
           //q=1
           if(q == 1){
             log_pi = log(pi(k));
             }
           
           //End, q > 1
           if((j == (q - 1)) & (q > 1)){
             
             arma::vec temp = trans(delta.submat(k, 0, k, (j - 1)));
             log_pi = log(prod(1.00 - temp)*pi(k));
             
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
      
      delta(k,j) = as<double>(Rcpp::rbinom(1,
                                           1,
                                           probs(1)));
      eta_full((k*q) + j) = A11_old*w1_old(k)*delta(k,j);
       
      }
  
   }

return Rcpp::List::create(Rcpp::Named("delta") = delta,
                          Rcpp::Named("eta_full") = eta_full);

}



