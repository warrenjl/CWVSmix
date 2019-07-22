#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w1_update(int n,
                     int p,
                     int q,
                     int m,
                     arma::mat x,
                     arma::mat z,
                     arma::vec w,
                     arma::vec gamma,
                     arma::vec beta,
                     arma::mat Lambda,
                     arma::mat delta,
                     arma::mat delta_star,
                     arma::vec w2_old,
                     double A11_old,
                     double A22_old,
                     double A21_old,
                     arma::mat corr_inv1){

arma::mat ident(m, m); ident.eye();
arma::vec delta_star_piece(m); delta_star_piece.fill(0.00); 
for(int j = 0; j < q; ++ j){
   delta_star_piece = delta_star_piece +
                      (delta_star.col(j) - A22_old*w2_old);
   }

arma::mat delta_diag(m*q, m); delta_diag.fill(0.00); 
for(int j = 0; j < m; ++ j){
  delta_diag.col(j).subvec((j*q), (q*(j + 1) - 1)) = trans(delta.row(j));
  }

arma::mat cov_piece = A11_old*z*kron(ident, Lambda)*delta_diag;
arma::mat cov_piece_trans = trans(cov_piece);

arma::mat w_mat(n, m);
for(int j = 0; j < m; ++j){
   w_mat.col(j) = w;
   }

arma::mat cov_w1 = inv_sympd(cov_piece_trans*(w_mat%cov_piece) + 
                             A21_old*A21_old*q*ident + 
                             corr_inv1);

arma::vec mean_w1 = cov_w1*(cov_piece_trans*(w%(gamma - x*beta)) + 
                            A21_old*delta_star_piece);

arma::mat ind_norms = arma::randn(1, m);
arma::vec w1 = mean_w1 +                   
               trans(ind_norms*arma::chol(cov_w1));

arma::vec eta_full(m*q); eta_full.fill(0.00); 
for(int j = 0; j < m; ++ j){
   eta_full.subvec((j*q), (q*(j + 1) - 1)) = A11_old*w1(j)*trans(delta.row(j));
   }

return Rcpp::List::create(Rcpp::Named("w1") = w1,
                          Rcpp::Named("eta_full") = eta_full);

}

