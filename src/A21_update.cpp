#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double A21_update(int q,
                  double sigma2_A,
                  arma::mat delta_star,
                  arma::vec w1,
                  arma::vec w2,
                  double A22){

double A21_var = 1.00/(q*sum(w1%w1) + (1.00/sigma2_A));

double A21_mean_piece = 0.00;
for(int j = 0; j < q; ++ j){
   A21_mean_piece = A21_mean_piece + 
                    sum(w1%(delta_star.col(j) - A22*w2));
   }
double A21_mean = A21_var*A21_mean_piece;

double A21 = R::rnorm(A21_mean, 
                      sqrt(A21_var));

return A21;

}



