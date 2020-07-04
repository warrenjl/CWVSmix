#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List A11_update(double A11_old,
                      int n,
                      int p,
                      int m,
                      arma::mat x,
                      arma::mat z,
                      double sigma2_A,
                      arma::vec w,
                      arma::vec gamma,
                      arma::vec beta,
                      arma::vec delta,
                      arma::vec w1,
                      arma::mat risk_sum,
                      double metrop_var_A11_trans,
                      int acctot_A11_trans){
  
arma::mat ident(m, m); ident.eye();
arma::vec eta_sub(m); eta_sub.fill(0.00); 
for(int j = 0; j < m; ++ j){
   eta_sub(j) = w1(j)*delta(j);
   }

/*Second*/
double A11_trans_old = log(A11_old);
arma::vec eta_full_old = A11_old*eta_sub;
  
arma::vec mean_piece_old = gamma - 
                           x*beta - 
                           risk_sum*eta_full_old;

double second = -0.50*dot(mean_piece_old, w%mean_piece_old) + 
                -0.50*(1.00/sigma2_A)*(A11_trans_old*A11_trans_old);

/*First*/
double A11_trans = R::rnorm(A11_trans_old, 
                            sqrt(metrop_var_A11_trans));
double A11 = exp(A11_trans);
arma::vec eta_full = A11*eta_sub;
   
arma::vec mean_piece = gamma - 
                       x*beta - 
                       risk_sum*eta_full;

double first = -0.50*dot(mean_piece, w%mean_piece) + 
               -0.50*(1.00/sigma2_A)*(A11_trans*A11_trans);

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
     
  A11 = A11_old;
  eta_full = eta_full_old;
  acc = 0;
     
  }
acctot_A11_trans = acctot_A11_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("A11") = A11,
                          Rcpp::Named("eta_full") = eta_full,
                          Rcpp::Named("acctot_A11_trans") = acctot_A11_trans);

}

