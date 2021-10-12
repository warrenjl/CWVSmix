#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w_update(int n,
                    int p,
                    int m,
                    arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta_old,
                    arma::vec eta_full,
                    arma::mat risk_sum){
  
arma::mat ident(m, m); ident.eye();

arma::vec mean_w = x*beta_old + 
                   risk_sum*eta_full;

arma::vec input(1); input.fill(1.00);
arma::vec w = rcpp_pgdraw(input,
                          mean_w);

arma::vec gamma = (y - 0.50)/w;

return Rcpp::List::create(Rcpp::Named("w") = w,
                          Rcpp::Named("gamma") = gamma);

}


















































































































