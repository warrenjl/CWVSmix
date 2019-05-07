#include "RcppArmadillo.h"
#include "CWMix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat eta_update(int p,
                     int q,
                     arma::mat x, 
                     arma::mat z,
                     arma::vec w,
                     arma::vec gamma,
                     arma::vec beta,
                     arma::vec sigma2_eta_old,
                     Rcpp::List temporal_corr_info,
                     arma::mat Lambda_old){

int m = z.n_cols/p;  
int n = w.size();
arma::mat w_mat(n, (m*q));
for(int j = 0; j < (m*q); ++ j){
   w_mat.col(j) = w;
   }

arma::mat ident(m,m); ident.eye();

arma::mat cov_piece = z*kron(ident, Lambda_old);
arma::mat cov_piece_trans = trans(cov_piece);

arma::mat Sigma0_inv((m*q), (m*q)); Sigma0_inv.fill(0.00);
for(int j = 0; j < q; ++ j){
  
   arma::mat temp(m,m); temp.fill(1.00/sigma2_eta_old(j));
   Rcpp::List temporal_corr_info_temp = temporal_corr_info[j];
   Sigma0_inv.submat(j*m, j*m, (m*(j + 1) - 1), (m*(j + 1) - 1)) = temp%Rcpp::as<arma::mat>(temporal_corr_info_temp[0]);
   
   }
arma::vec sort_set(m*q); sort_set.fill(0.00);
for(int j = 0; j < m; ++ j){
  sort_set.subvec(j*q, (q*(j + 1) - 1)) = regspace(j, m, (q*m - 1));
  }

arma::uvec sort_set_final = conv_to<arma::uvec>::from(sort_set);
arma::mat Sigma1_inv = Sigma0_inv.submat(sort_set_final, sort_set_final);

arma::mat cov_eta = inv_sympd(cov_piece_trans*(w_mat%cov_piece) + 
                              Sigma1_inv);

arma::vec mean_eta = cov_eta*(cov_piece_trans*(w%(gamma - x*beta)));

arma::mat ind_norms = arma::randn(1, (m*q));
arma::vec eta_full = mean_eta + 
                     trans(ind_norms*arma::chol(cov_eta));

arma::mat eta(q,m); eta.fill(0.00);
for(int j = 0; j < m; ++ j){
   eta.col(j) = eta_full.subvec(j*q, (q*(j + 1) - 1));
   } 

return(eta);

}



