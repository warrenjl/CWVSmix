#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List Lambda_update(arma::mat Lambda_old,
                         int ind,
                         int p,
                         int q,
                         arma::mat x,
                         arma::mat z,
                         double alpha_Lambda,
                         arma::vec w,
                         arma::vec gamma,
                         arma::vec beta,
                         arma::vec eta_full,
                         double metrop_scale_Lambda,
                         int acctot_Lambda){
  
int m = z.n_cols/p;  
arma::mat ident(m, m); ident.eye();
double epsilon = 0.01;

/*First*/
arma::vec Lambda_vec_temp(p - ind); Lambda_vec_temp.fill(0.00);
for(int j = ind; j < p; ++ j){
   Lambda_vec_temp(j - ind) = R::rgamma((metrop_scale_Lambda*Lambda_old(j, ind) + epsilon),
                                        1.00);
   }
arma::mat Lambda = Lambda_old;
for(int j = ind; j < p; ++ j){
   Lambda(j, ind) = Lambda_vec_temp(j - ind)/sum(Lambda_vec_temp);
   }

arma::vec mh_piece(p - ind); mh_piece.fill(0.00);
for(int j = 0; j < (p - ind); ++ j){
   mh_piece(j) = log(tgamma(metrop_scale_Lambda*Lambda((j + ind), ind) + epsilon));
   }

double first = -0.50*dot((gamma - x*beta - z*((kron(ident, Lambda))*eta_full)), w%(gamma - x*beta - z*((kron(ident, Lambda))*eta_full))) + 
               (alpha_Lambda - 1.00)*sum(log(Lambda.col(ind).subvec(ind, (p-1)))) + 
               log(tgamma(sum(metrop_scale_Lambda*Lambda.col(ind).subvec(ind, (p-1)) + epsilon))) - 
               sum(mh_piece) +
               sum((metrop_scale_Lambda*Lambda.col(ind).subvec(ind, (p-1)) + epsilon - 1.00)%log(Lambda_old.col(ind).subvec(ind, (p-1))));

/*Second*/
arma::vec mh_piece_old(p - ind); mh_piece_old.fill(0.00);
for(int j = 0; j < (p - ind); ++ j){
  mh_piece_old(j) = log(tgamma(metrop_scale_Lambda*Lambda_old((j + ind), ind) + epsilon));
  }

double second = -0.50*dot((gamma - x*beta - z*((kron(ident, Lambda_old))*eta_full)), w%(gamma - x*beta - z*((kron(ident, Lambda_old))*eta_full))) + 
                (alpha_Lambda - 1.00)*sum(log(Lambda_old.col(ind).subvec(ind, (p-1)))) + 
                log(tgamma(sum(metrop_scale_Lambda*Lambda_old.col(ind).subvec(ind, (p-1)) + epsilon))) - 
                sum(mh_piece_old) +
                sum((metrop_scale_Lambda*Lambda_old.col(ind).subvec(ind, (p-1)) + epsilon - 1.00)%log(Lambda.col(ind).subvec(ind, (p-1))));

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if((ratio < R::runif(0.00, 1.00))){
  Lambda = Lambda_old;
  acc = 0;
  }
acctot_Lambda = acctot_Lambda + 
                acc;

return Rcpp::List::create(Rcpp::Named("Lambda_vec") = Lambda.col(ind),
                          Rcpp::Named("acctot_Lambda") = acctot_Lambda);

}
                 
  
