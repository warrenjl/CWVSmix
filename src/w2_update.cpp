#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec w2_update(int p,
                    int q,
                    int m,
                    arma::mat z,
                    arma::mat delta_star,
                    arma::vec w1,
                    double A22_old,
                    double A21_old,
                    arma::mat corr_inv2){
  
arma::mat ident(m, m); ident.eye();
arma::vec delta_star_piece(m); delta_star_piece.fill(0.00); 
for(int j = 0; j < q; ++ j){
   delta_star_piece = delta_star_piece +
                      (delta_star.col(j) - A21_old*w1);
   }
  
arma::mat cov_w2 = inv_sympd(A22_old*A22_old*q*ident + 
                             corr_inv2);

arma::vec mean_w2 = cov_w2*(A22_old*delta_star_piece);

arma::mat ind_norms = arma::randn(1, m);
arma::vec w2 = mean_w2 + 
               trans(ind_norms*arma::chol(cov_w2));

return(w2);

}






