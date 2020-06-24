#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec delta_star_update(int m,
                            arma::vec delta,
                            arma::vec w1_old,
                            arma::vec w2_old,
                            double A22_old,
                            double A21_old){
  
arma::vec delta_star(m); delta_star.fill(0.00);

arma::vec alpha = A21_old*w1_old +
                  A22_old*w2_old;

for(int j = 0; j < m; ++ j){
   
   if(delta(j) == 1.00){
     delta_star(j) = rnorm_trunc(alpha(j),
                                 1.00,
                                 0.00,
                                 datum::inf);
     }
   
    
   if(delta(j) == 0.00){
     delta_star(j) = rnorm_trunc(alpha(j),
                                 1.00,
                                 -datum::inf,
                                 0.00);
     }
        
   }
   
return(delta_star);

}



