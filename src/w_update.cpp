#include "RcppArmadillo.h"
#include "CWMix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w_update(int p,
                    int q,
                    arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta_old,
                    arma::mat eta_old,
                    arma::mat Lambda_old){
  
int m = z.n_cols/p;  
arma::mat ident(m,m); ident.eye();

arma::vec eta_old_full(q*m); eta_old_full.fill(0.00);
for(int j = 0; j < m; ++ j){
   eta_old_full.subvec(j*q, (q*(j + 1) - 1)) = eta_old.col(j);
   } 

arma::vec mean_w = x*beta_old + 
                   z*((kron(ident, Lambda_old))*eta_old_full);

arma::vec w = rcpp_pgdraw(1.00,
                          mean_w);

arma::vec gamma = (y - 0.50)/w;

return Rcpp::List::create(Rcpp::Named("w") = w,
                          Rcpp::Named("gamma") = gamma);

}
































































