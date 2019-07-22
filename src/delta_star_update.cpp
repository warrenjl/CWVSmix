#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat delta_star_update(int q,
                            int m,
                            arma::mat delta,
                            arma::vec w1_old,
                            arma::vec w2_old,
                            double A22_old,
                            double A21_old){
  
arma::mat delta_star(m, q); delta_star.fill(0.00);

arma::vec alpha = A21_old*w1_old +
                  A22_old*w2_old;

for(int j = 0; j < q; ++ j){
  
   for(int k = 0; k < m; ++ k){
   
      if(delta(k, j) == 1.00){
        delta_star(k, j) = rnorm_trunc(alpha(k),
                                       1.00,
                                       0.00,
                                       datum::inf);
        }
   
      //All Other
      if(j > 0){
        
        if((delta(k, j) == 0.00) & (sum(delta.row(k).subvec(0, (j - 1))) > 0.00)){
          delta_star(k, j) = R::rnorm(alpha(k),
                                      sqrt(1.00));
          }
      
        if((delta(k, j) == 0.00) & (sum(delta.row(k).subvec(0, (j - 1))) == 0.00)){
          delta_star(k, j) = rnorm_trunc(alpha(k),
                                         1.00,
                                         -datum::inf,
                                         0.00);
          }
        
        }
      
      //Start
      if(j == 0){
        
        if(delta(k, j) == 0.00){
          delta_star(k, j) = rnorm_trunc(alpha(k),
                                         1.00,
                                         -datum::inf,
                                         0.00);
          }
        
        }
   
      }
   
   }
   
return(delta_star);

}



