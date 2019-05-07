#include "RcppArmadillo.h"
#include "CWMix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec beta_update(int p,
                      int q,
                      arma::mat x, 
                      arma::mat z,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec gamma,
                      arma::mat eta_old,
                      arma::mat Lambda_old){

int m = z.n_cols/p;  
arma::mat ident(m, m); ident.eye();
int p_x = x.n_cols;
int n = w.size();
arma::mat w_mat(n, p_x);
for(int j = 0; j < p_x; ++j){
   w_mat.col(j) = w;
   }

arma::vec eta_old_full(q*m); eta_old_full.fill(0.00);
for(int j = 0; j < m; ++j){
   eta_old_full.subvec(j*q, (q*(j + 1) - 1)) = eta_old.col(j);
   } 

arma::mat x_trans = trans(x);

arma::mat cov_beta = inv_sympd(x_trans*(w_mat%x) + 
                               (1.00/sigma2_beta)*eye(p_x, p_x));

arma::vec mean_beta = cov_beta*(x_trans*(w%(gamma - z*((kron(ident, Lambda_old))*eta_old_full))));

arma::mat ind_norms = arma::randn(1, p_x);
arma::vec beta = mean_beta + 
                 trans(ind_norms*arma::chol(cov_beta));

return(beta);

}



