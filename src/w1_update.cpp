#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w1_update(int n,
                     int p,
                     int m,
                     arma::mat x,
                     arma::mat z,
                     arma::vec w,
                     arma::vec gamma,
                     arma::vec beta,
                     arma::vec delta,
                     arma::vec delta_star,
                     arma::vec w2_old,
                     double A11_old,
                     double A22_old,
                     double A21_old,
                     arma::mat risk_sum,
                     arma::mat corr_inv1){

arma::mat ident(m, m); ident.eye();
arma::vec delta_star_piece = (delta_star - A22_old*w2_old);

arma::mat delta_diag(m, m); delta_diag.fill(0.00); 
for(int j = 0; j < m; ++ j){
  delta_diag(j, j) = delta(j);
  }

arma::mat cov_piece = A11_old*(risk_sum*delta_diag);
arma::mat cov_piece_trans = trans(cov_piece);

arma::mat w_mat(n, m);
for(int j = 0; j < m; ++j){
   w_mat.col(j) = w;
   }

arma::mat cov_w1 = inv_sympd(cov_piece_trans*(w_mat%cov_piece) + 
                             A21_old*A21_old*ident + 
                             corr_inv1);

arma::vec mean_w1 = cov_w1*(cov_piece_trans*(w%(gamma - x*beta)) + 
                            A21_old*delta_star_piece);

arma::mat ind_norms = arma::randn(1, 
                                  m);
arma::vec w1 = mean_w1 +                   
               trans(ind_norms*arma::chol(cov_w1));

arma::vec eta_full(m); eta_full.fill(0.00); 
for(int j = 0; j < m; ++ j){
   eta_full(j) = A11_old*w1(j)*delta(j);
   }

return Rcpp::List::create(Rcpp::Named("w1") = w1,
                          Rcpp::Named("eta_full") = eta_full);

}

