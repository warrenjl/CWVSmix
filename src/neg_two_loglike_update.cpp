#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update(int p,
                              int q,
                              arma::vec y,
                              arma::mat x,
                              arma::mat z, 
                              arma::vec beta,
                              arma::mat Lambda,
                              arma::vec eta_full){

int n = y.size();
int m = z.n_cols/p;  
arma::mat ident(m, m); ident.eye();
arma::vec dens(n); dens.fill(0.00);

arma::vec logit_probs = x*beta + 
                        z*((kron(ident, Lambda))*eta_full);

arma::vec probs = exp(logit_probs)/(1.00 + exp(logit_probs));

for(int j = 0; j < n; ++ j){
   dens(j) = R::dbinom(y(j),
                       1,
                       probs(j),
                       TRUE);
   }

double neg_two_loglike = -2.00*sum(dens);

return neg_two_loglike;

}

























































