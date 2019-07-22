#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec beta_update(int n,
                      int p,
                      int q,
                      int m,
                      int p_x,
                      arma::mat x, 
                      arma::mat z,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec gamma,
                      arma::mat Lambda_old,
                      arma::vec eta_full){

arma::mat ident(m, m); ident.eye();
arma::mat w_mat(n, p_x);
for(int j = 0; j < p_x; ++ j){
   w_mat.col(j) = w;
   }

arma::mat x_trans = trans(x);

arma::mat cov_beta = inv_sympd(x_trans*(w_mat%x) + 
                               (1.00/sigma2_beta)*eye(p_x, p_x));

arma::vec mean_beta = cov_beta*(x_trans*(w%(gamma - z*((kron(ident, Lambda_old))*eta_full))));

arma::mat ind_norms = arma::randn(1, 
                                  p_x);

arma::vec beta = mean_beta + 
                 trans(ind_norms*arma::chol(cov_beta));

return(beta);

}



