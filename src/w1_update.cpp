#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w1_update(int p,
                     int q,
                     arma::mat x,
                     arma::mat z,
                     arma::vec w,
                     arma::vec gamma,
                     arma::vec beta,
                     arma::mat Lambda,
                     arma::mat delta,
                     arma::mat delta_star,
                     arma::mat w2_old,
                     arma::vec A11_old,
                     arma::vec A21_old,
                     arma::vec A22_old,
                     Rcpp::List temporal_corr_info1){

int m = z.n_cols/p;
arma::mat ident(m, m); ident.eye();
int n = w.size();

arma::vec delta_diag(m*q); delta_diag.fill(0.00); 
arma::vec w2_full(m*q); w2_full.fill(0.00);
arma::vec A11_diag(m*q); A11_diag.fill(0.00);
arma::vec A22_diag(m*q); A22_diag.fill(0.00);
arma::vec A21_diag(m*q); A21_diag.fill(0.00);
for(int j = 0; j < m; ++ j){
  
   delta_diag.subvec((j*q), (q*(j + 1) - 1)) = delta.row(j);
   w2_full.subvec((j*q), (q*(j + 1) - 1)) = w2_old.row(j);
   A11_diag.subvec((j*q), (q*(j + 1) - 1)) = A11_old;
   A22_diag.subvec((j*q), (q*(j + 1) - 1)) = A22_old;
   A21_diag.subvec((j*q), (q*(j + 1) - 1)) = A21_old;
   
   }

arma::mat Sigma0_inv((m*q), (m*q)); Sigma0_inv.fill(0.00);
for(int j = 0; j < q; ++ j){
  
   Rcpp::List temporal_corr_info1_temp = temporal_corr_info1[j];
   Sigma0_inv.submat((j*m), (j*m), (m*(j + 1) - 1), (m*(j + 1) - 1)) = Rcpp::as<arma::mat>(temporal_corr_info1_temp[0]);
  
   }
arma::vec sort_set(m*q); sort_set.fill(0.00);
for(int j = 0; j < m; ++ j){
   sort_set.subvec((j*q), (q*(j + 1) - 1)) = regspace(j, m, ((q*m) - 1));
   }

arma::uvec sort_set_final = conv_to<arma::uvec>::from(sort_set);
arma::mat Sigma1_inv = Sigma0_inv.submat(sort_set_final, sort_set_final);

arma::mat cov_piece = z*kron(ident, Lambda)*diagmat(delta_diag%A11_diag);
arma::mat cov_piece_trans = trans(cov_piece);

arma::mat w_mat(n, (m*q));
for(int j = 0; j < (m*q); ++j){
   w_mat.col(j) = w;
   }

arma::mat cov_w1 = inv_sympd(cov_piece_trans*(w_mat%cov_piece) + 
                             diagmat(A21_diag%A21_diag) + 
                             Sigma1_inv);

arma::vec mean_w1 = cov_w1*(cov_piece_trans*(w%(gamma - x*beta)) + 
                            diagmat(A21_diag)*(delta_star - A22_diag%w2_full));

arma::mat ind_norms = arma::randn(1, (m*q));
arma::vec w1_full = mean_w1 +                   
                    trans(ind_norms*arma::chol(cov_w1));

arma::vec eta_full = (delta_diag%A11_diag%w1_full);

arma::mat w1(m, q); w1.fill(0.00);
for(int j = 0; j < q; ++ j){
   arma::vec subset = regspace(j, q, ((m*q) - 1));
   arma::uvec subset_final = conv_to<arma::uvec>::from(subset);
   w1.col(j) = w1_full.elem(subset_final);
   }

return Rcpp::List::create(Rcpp::Named("w1") = w1,
                          Rcpp::Named("eta_full") = eta_full);

}

