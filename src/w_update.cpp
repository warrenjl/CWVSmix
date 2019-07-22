#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w_update(int p,
                    int q,
                    int m,
                    arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta_old,
                    arma::mat Lambda_old,
                    arma::vec eta_full){
  
arma::mat ident(m, m); ident.eye();

arma::vec mean_w = x*beta_old + 
                   z*((kron(ident, Lambda_old))*eta_full);

arma::vec w = rcpp_pgdraw(1.00,
                          mean_w);

arma::vec gamma = (y - 0.50)/w;

return Rcpp::List::create(Rcpp::Named("w") = w,
                          Rcpp::Named("gamma") = gamma);

}
































































