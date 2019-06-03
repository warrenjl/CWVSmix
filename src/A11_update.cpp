#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List A11_update(arma::vec A11_old,
                      int p,
                      int q,
                      arma::mat x,
                      arma::mat z,
                      double sigma2_A,
                      arma::vec w,
                      arma::vec gamma,
                      arma::vec beta,
                      arma::mat Lambda,
                      arma::mat delta,
                      arma::mat w1,
                      arma::vec metrop_var_A11_trans,
                      arma::vec acctot_A11_trans){
  
int m = z.n_cols/p;  
arma::mat ident(m, m); ident.eye();
arma::vec A11 = A11_old;

arma::vec delta_diag(m*q); delta_diag.fill(0.00); 
arma::vec w1_full(m*q); w1_full.fill(0.00);
for(int j = 0; j < m; ++ j){
  
   delta_diag.subvec((j*q), (q*(j + 1) - 1)) = trans(delta.row(j));
   w1_full.subvec((j*q), (q*(j + 1) - 1)) = trans(w1.row(j));
  
   }

arma::vec eta_full(m*q); eta_full.fill(0.00);
for(int j = 0; j < q; ++ j){
  
   /*Second*/
   double A11_trans_old = log(A11_old(j));
  
   arma::vec A11_diag(m*q); A11_diag.fill(0.00);
   for(int k = 0; k < m; ++ k){
      A11_diag.subvec((k*q), (q*(k + 1) - 1)) = A11_old;
      }
   arma::vec eta_full_old = (delta_diag%A11_diag%w1_full);
  
   arma::vec mean_piece_old = gamma - 
                              x*beta - 
                              z*((kron(ident, Lambda))*eta_full_old);

   double second = -0.50*dot(mean_piece_old, w%mean_piece_old) - 
                   0.50*(1.00/sigma2_A)*(A11_trans_old*A11_trans_old);

   /*First*/
   double A11_trans = R::rnorm(A11_trans_old, 
                               sqrt(metrop_var_A11_trans(j)));
   A11(j) = exp(A11_trans);
   for(int k = 0; k < m; ++ k){
      A11_diag.subvec((k*q), (q*(k + 1) - 1)) = A11;
      }
   eta_full = (delta_diag%A11_diag%w1_full);
   
   arma::vec mean_piece = gamma - 
                          x*beta - 
                          z*((kron(ident, Lambda))*eta_full);

   double first = -0.50*dot(mean_piece, w%mean_piece) - 
                  0.50*(1.00/sigma2_A)*(A11_trans*A11_trans);

   /*Decision*/
   double ratio = exp(first - second);   
   int acc = 1;
   if(ratio < R::runif(0.00, 1.00)){
     A11(j) = A11_old(j);
     eta_full = eta_full_old;
     acc = 0;
     }
   acctot_A11_trans(j) = acctot_A11_trans(j) + 
                         acc;
   
   }

return Rcpp::List::create(Rcpp::Named("A11") = A11,
                          Rcpp::Named("eta_full") = eta_full,
                          Rcpp::Named("acctot_A11_trans") = acctot_A11_trans);

}



