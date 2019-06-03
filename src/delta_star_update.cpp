#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat delta_star_update(arma::mat delta,
                            arma::mat w1_old,
                            arma::mat w2_old,
                            arma::vec A22_old,
                            arma::vec A21_old){
  
int q = delta.n_cols;
int m = delta.n_rows;
arma::mat delta_star(m, q); delta_star.fill(0.00);
for(int j = 0; j < q; ++ j){
  
   arma::vec alpha = A21_old(j)*w1_old.col(j) +
                     A22_old(j)*w2_old.col(j);
  
   for(int k = 0; k < m; ++ k){
   
      if(delta(k,j) == 1.00){
        while(delta_star(k,j) <= 0.00){
             delta_star(k,j) = R::rnorm(alpha(k),
                                        sqrt(1.00));
             }
        }
   
      //All Other
      if(j > 0){
        
        if((delta(k,j) == 0.00) & (sum(delta.col(j - 1)) > 0.00)){
          while(delta_star(k,j) >= 0.00){
               delta_star(k,j) = R::rnorm(alpha(k),
                                          sqrt(1.00));
               }
          }
      
        if((delta(k,j) == 0.00) & (sum(delta.col(j - 1)) == 0.00)){
          delta_star(k,j) = R::rnorm(alpha(k),
                                     sqrt(1.00));
          }
        
        }
      
      //Start
      if(j == 0){
        
        if(delta(k,j) == 0.00){
          while(delta_star(k,j) >= 0.00){
               delta_star(k,j) = R::rnorm(alpha(k),
                                          sqrt(1.00));
               }
           }
        
        }
   
      }
   
   }
   
return(delta_star);

}



