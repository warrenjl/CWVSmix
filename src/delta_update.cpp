#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List delta_update(arma::vec delta_old,
                        int n,
                        int p,
                        int m,
                        arma::vec y,
                        arma::mat x,
                        arma::mat z,
                        arma::vec w,
                        arma::vec gamma,
                        arma::vec beta,
                        arma::vec w1_old,
                        arma::vec w2_old,
                        double A11_old,
                        double A22_old,
                        double A21_old,
                        arma::vec eta_full,
                        arma::mat risk_sum){
   
arma::mat ident(m, m); ident.eye();
   
arma::vec delta = delta_old;
   
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
   
for(int j = 0; j < m; ++ j){
         
   pieces.fill(0.00);
   for(int k = 0; k < 2; ++ k){
            
      delta(j) = k;
      eta_full(j) = A11_old*w1_old(j)*delta(j);
            
      if(k == 0){
        log_pi = log(1.00 - pi(j));
        }
            
      if(k == 1){
        log_pi = log(pi(j));
        }
            
      pieces(k) = -0.50*dot((gamma - x*beta - risk_sum*eta_full), w%(gamma - x*beta - risk_sum*eta_full)) +
                  log_pi;
            
      }
         
   probs.fill(0.00);
   for(int k = 0; k < 2; ++ k){
            
      probs(k) = 1.00/(sum(exp(pieces - pieces(k))));
            
      if(arma::is_finite(probs(k)) == 0.00){
        probs(k) = 0.00;  /*Computational Correction*/
        }
            
      }
         
   delta(j) = as<double>(Rcpp::rbinom(1,
                                      1,
                                      probs(1)));
         
   eta_full(j) = A11_old*w1_old(j)*delta(j);
      
   }
   
return Rcpp::List::create(Rcpp::Named("delta") = delta,
                          Rcpp::Named("eta_full") = eta_full);
   
}


