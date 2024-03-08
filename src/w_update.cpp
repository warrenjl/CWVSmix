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
                    arma::vec off_set,
                    arma::vec tri_als,
                    int likelihood_indicator,
                    int r,
                    arma::vec beta_old,
                    arma::vec eta_full,
                    arma::mat risk_sum){
  
arma::mat ident(m, m); ident.eye();

arma::vec mean_w = off_set +
                   x*beta_old + 
                   risk_sum*eta_full;

arma::vec input0 = tri_als;
arma::vec input2 = (r + y);

arma::vec w(n); w.fill(0.00);
arma::vec gamma(n); gamma.fill(0.00);

if(likelihood_indicator == 0){
  
  w = rcpp_pgdraw(input0,
                  mean_w);

  gamma = (y - 0.50*tri_als)/w;

  }

if(likelihood_indicator == 2){
  
  w = rcpp_pgdraw(input2,
                  mean_w);
  gamma = 0.50*(y - r)/w;
  
  } 

return Rcpp::List::create(Rcpp::Named("w") = w,
                          Rcpp::Named("gamma") = gamma);

}


















































































































