#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List CWVSmix(int mcmc_samples,
                   int p,
                   int q,
                   arma::vec y,
                   arma::mat x,
                   arma::mat z,
                   arma::mat metrop_var_Lambda,
                   double metrop_var_A11_trans,
                   double metrop_var_A22_trans,
                   double metrop_var_phi1_trans,
                   double metrop_var_phi2_trans,
                   Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> alpha_Lambda_prior = R_NilValue,
                   Rcpp::Nullable<double> sigma2_A_prior = R_NilValue,
                   Rcpp::Nullable<double> alpha_phi1_prior = R_NilValue,
                   Rcpp::Nullable<double> beta_phi1_prior = R_NilValue,
                   Rcpp::Nullable<double> alpha_phi2_prior = R_NilValue,
                   Rcpp::Nullable<double> beta_phi2_prior = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> Lambda_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> delta_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> w1_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> w2_init = R_NilValue,
                   Rcpp::Nullable<double> A11_init = R_NilValue,
                   Rcpp::Nullable<double> A22_init = R_NilValue,
                   Rcpp::Nullable<double> A21_init = R_NilValue,
                   Rcpp::Nullable<double> phi1_init = R_NilValue,
                   Rcpp::Nullable<double> phi2_init = R_NilValue){

//Defining Parameters and Quantities of Interest
int n = y.size();
int m = z.n_cols/p;
double max_time = (m - 1);
int p_x = x.n_cols;

arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
Rcpp::List Lambda(mcmc_samples);
Rcpp::List delta(mcmc_samples);
for(int j = 0; j < mcmc_samples; ++ j){
  
   arma::mat Lambda_temp(p, q); Lambda_temp.fill(0.00);
   Lambda[j] = Lambda_temp;
   
   arma::mat delta_temp(m, q); delta_temp.fill(0.00);
   delta[j] = delta_temp;
  
   }
arma::mat w1(m, mcmc_samples); w1.fill(0.00);
arma::mat w2(m, mcmc_samples); w2.fill(0.00);
arma::vec A11(mcmc_samples); A11.fill(0.00);
arma::vec A22(mcmc_samples); A22.fill(0.00);
arma::vec A21(mcmc_samples); A21.fill(0.00);
arma::vec phi1(mcmc_samples); phi1.fill(0.00);
arma::vec phi2(mcmc_samples); phi2.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Prior Information
double sigma2_beta = 10000.00;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }

arma::mat alpha_Lambda(p, q); alpha_Lambda.fill(0.00);
for(int j = 0; j < q; ++ j){
   for(int k = j; k < p; ++ k){
      alpha_Lambda(k, j) = 1.00;
      }
   }
if(alpha_Lambda_prior.isNotNull()){
  alpha_Lambda = Rcpp::as<arma::mat>(alpha_Lambda_prior);
  }

double sigma2_A = 1.00;
if(sigma2_A_prior.isNotNull()){
  sigma2_A = Rcpp::as<double>(sigma2_A_prior);
  }

double alpha_phi1 = 1.00;  
if(alpha_phi1_prior.isNotNull()){
  alpha_phi1 = Rcpp::as<double>(alpha_phi1_prior);
  }

double beta_phi1 = 1.00;
if(beta_phi1_prior.isNotNull()){
  beta_phi1 = Rcpp::as<double>(beta_phi1_prior);
  }

double alpha_phi2 = 1.00;  
if(alpha_phi2_prior.isNotNull()){
  alpha_phi2 = Rcpp::as<double>(alpha_phi2_prior);
  }

double beta_phi2 = 1.00;
if(beta_phi2_prior.isNotNull()){
  beta_phi2 = Rcpp::as<double>(beta_phi2_prior);
  }

//Initial Values
beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

arma::mat Lambda_temp(p, q); Lambda_temp.fill(0.00);
for(int j = 0; j < q; ++ j){
   for(int k = j; k < p; ++ k){
      Lambda_temp(k, j) = alpha_Lambda(k, j)/sum(alpha_Lambda.col(j).subvec(j, (p - 1)));
      }
   }
if(Lambda_init.isNotNull()){
  Lambda_temp = Rcpp::as<arma::mat>(Lambda_init);
  }
Lambda[0] = Lambda_temp;
arma::mat Lambda_star(p, q); Lambda_star.fill(0.00);
for(int j = 0; j < q; ++ j){
   for(int k = j; k < p; ++ k){
      Lambda_star(k, j) = Lambda_temp(k, j);
      }
   }

arma::mat delta_temp(m, q); delta_temp.fill(0.00); delta_temp.col(0).fill(1.00);
if(delta_init.isNotNull()){
  delta_temp = Rcpp::as<arma::mat>(delta_init);
  }
delta[0] = delta_temp;

w1.col(0).fill(0.00);
if(w1_init.isNotNull()){
  w1.col(0) = Rcpp::as<arma::vec>(w1_init);
  }

w2.col(0).fill(0.00);
if(w2_init.isNotNull()){
  w2.col(0) = Rcpp::as<arma::vec>(w2_init);
  }

A11(0) = 1.00;
if(A11_init.isNotNull()){
  A11(0) = Rcpp::as<double>(A11_init);
  }

A22(0) = 1.00;
if(A22_init.isNotNull()){
  A22(0) = Rcpp::as<double>(A22_init);
  }

A21(0) = 0.00;
if(A21_init.isNotNull()){
  A21(0) = Rcpp::as<double>(A21_init);
  }

phi1(0) = -log(0.05)/max_time;  //Effective range equal to largest temporal distance in dataset (strong temporal correlation)
if(phi1_init.isNotNull()){
  phi1(0) = Rcpp::as<double>(phi1_init);
  }

phi2(0) = -log(0.05)/max_time;  //Effective range equal to largest temporal distance in dataset (strong temporal correlation)
if(phi2_init.isNotNull()){
  phi2(0) = Rcpp::as<double>(phi2_init);
  }

Rcpp::List temporal_corr_info1 = temporal_corr_fun(m, 
                                                   phi1(0));
  
Rcpp::List temporal_corr_info2 = temporal_corr_fun(m, 
                                                   phi2(0));

delta_temp = Rcpp::as<arma::mat>(delta[0]);
arma::vec eta_full(m*q); eta_full.fill(0.00); 
for(int j = 0; j < m; ++ j){
   eta_full.subvec((j*q), (q*(j + 1) - 1)) = A11(0)*w1(j, 0)*trans(delta_temp.row(j));
   }

neg_two_loglike(0) = neg_two_loglike_update(n,
                                            p,
                                            q,
                                            m,
                                            y,
                                            x,
                                            z, 
                                            beta.col(0),
                                            Lambda[0],
                                            eta_full);

//Metropolis Settings
arma::mat acctot_Lambda(p, q); acctot_Lambda.fill(0);
arma::vec acctot_Lambda_vec(q*p + q - q*(q + 1)/2); acctot_Lambda_vec.fill(0);
int acctot_A11_trans = 0;
int acctot_A22_trans = 0; 
int acctot_phi1_trans = 0;
int acctot_phi2_trans = 0; 

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++ j){
  
   //w Update
   Rcpp::List w_output = w_update(p,
                                  q,
                                  m,
                                  y,
                                  x,
                                  z,
                                  beta.col(j-1),
                                  Lambda[j-1],
                                  eta_full);
  
   arma::vec w = w_output[0];
   arma::vec gamma = w_output[1];
  
   //beta Update
   beta.col(j) = beta_update(n,
                             p,
                             q,
                             m,
                             p_x,
                             x, 
                             z,
                             sigma2_beta,
                             w,
                             gamma,
                             Lambda[j-1],
                             eta_full);
   
   //Lambda Update
   arma::mat Lambda_temp = Lambda[j-1];
   int counter = 0;
   for(int k = 0; k < q; ++ k){
     
      Rcpp::List Lambda_output = Lambda_update(Lambda_star.col(k),
                                               Lambda_temp,
                                               k,
                                               p,
                                               q,
                                               m,
                                               x,
                                               z,
                                               alpha_Lambda.col(k),
                                               w,
                                               gamma,
                                               beta.col(j),
                                               eta_full,
                                               metrop_var_Lambda.col(k),
                                               acctot_Lambda.col(k));
     
      Lambda_star.col(k) = Rcpp::as<arma::vec>(Lambda_output[0]);
      Lambda_temp.col(k) = Rcpp::as<arma::vec>(Lambda_output[1]);
      acctot_Lambda.col(k) = Rcpp::as<arma::vec>(Lambda_output[2]);
      
      acctot_Lambda_vec.subvec((0 + counter), (p - 1 - k + counter)) = acctot_Lambda.col(k).subvec(k, (p - 1));
      counter = counter + 
                (p - k - 1);
      
     
      }
   Lambda[j] = Lambda_temp;
   
   //delta Update
   Rcpp::List delta_output = delta_update(delta[j-1],
                                          p,
                                          q,
                                          m,
                                          y,
                                          x,
                                          z,
                                          w,
                                          gamma,
                                          beta.col(j),
                                          Lambda[j],
                                          w1.col(j-1),
                                          w2.col(j-1),
                                          A11(j-1),
                                          A22(j-1),
                                          A21(j-1));
        
   delta[j] = Rcpp::as<arma::mat>(delta_output[0]);
   eta_full = Rcpp::as<arma::vec>(delta_output[1]);
   
   //delta_star Update
   arma::mat delta_star = delta_star_update(q,
                                            m,
                                            delta[j],
                                            w1.col(j-1),
                                            w2.col(j-1),
                                            A22(j-1),
                                            A21(j-1));
   
   //w1 Update
   Rcpp::List w1_output = w1_update(n,
                                    p,
                                    q,
                                    m,
                                    x,
                                    z,
                                    w,
                                    gamma,
                                    beta.col(j),
                                    Lambda[j],
                                    delta[j],
                                    delta_star,
                                    w2.col(j-1),
                                    A11(j-1),
                                    A22(j-1),
                                    A21(j-1),
                                    temporal_corr_info1[0]);
   
   w1.col(j) = Rcpp::as<arma::vec>(w1_output[0]);
   eta_full = Rcpp::as<arma::vec>(w1_output[1]);
   
   //w2 Update
   w2.col(j) = w2_update(p,
                         q,
                         m,
                         z,
                         delta_star,
                         w1.col(j),
                         A22(j-1),
                         A21(j-1),
                         temporal_corr_info2[0]);
   
   //A11 Update
   Rcpp::List A11_output = A11_update(A11(j-1),
                                      p,
                                      q,
                                      m,
                                      x,
                                      z,
                                      sigma2_A,
                                      w,
                                      gamma,
                                      beta.col(j),
                                      Lambda[j],
                                      delta[j],
                                      w1.col(j),
                                      metrop_var_A11_trans,
                                      acctot_A11_trans);
   
   A11(j) = Rcpp::as<double>(A11_output[0]);
   eta_full = Rcpp::as<arma::vec>(A11_output[1]);
   acctot_A11_trans = Rcpp::as<double>(A11_output[2]);
   
   //A22 Update
   Rcpp::List A22_output = A22_update(A22(j-1),
                                      q,
                                      m,
                                      sigma2_A,
                                      delta_star,
                                      w1.col(j),
                                      w2.col(j),
                                      A21(j-1),
                                      metrop_var_A22_trans,
                                      acctot_A22_trans);
   
   A22(j) = Rcpp::as<double>(A22_output[0]);
   acctot_A22_trans = Rcpp::as<double>(A22_output[1]);
   
   //A21 Update
   A21(j) = A21_update(q,
                       sigma2_A,
                       delta_star,
                       w1.col(j),
                       w2.col(j),
                       A22(j));
   
   //phi1 Update
   Rcpp::List phi1_output = phi_update(phi1(j-1),
                                       m,
                                       alpha_phi1,
                                       beta_phi1,
                                       w1.col(j),
                                       temporal_corr_info1,
                                       metrop_var_phi1_trans,
                                       acctot_phi1_trans);
     
   phi1(j) = Rcpp::as<double>(phi1_output[0]);
   acctot_phi1_trans = Rcpp::as<int>(phi1_output[1]);
   temporal_corr_info1 = phi1_output[2];
   
   //phi2 Update
   Rcpp::List phi2_output = phi_update(phi2(j-1),
                                       m,
                                       alpha_phi2,
                                       beta_phi2,
                                       w2.col(j),
                                       temporal_corr_info2,
                                       metrop_var_phi2_trans,
                                       acctot_phi2_trans);
     
   phi2(j) = Rcpp::as<double>(phi2_output[0]);
   acctot_phi2_trans = Rcpp::as<int>(phi2_output[1]);
   temporal_corr_info2 = phi2_output[2];
  
   //neg_two_loglike Update
   neg_two_loglike(j) = neg_two_loglike_update(n,
                                               p,
                                               q,
                                               m,
                                               y,
                                               x,
                                               z, 
                                               beta.col(j),
                                               Lambda[j],
                                               eta_full);
   
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
  
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
     
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     
     if((q == 1) & (p == 2)){
       
       double accrate_Lambda_min = round(100*(min(acctot_Lambda_vec)/(double)j));
       Rcpp::Rcout << "Lambda Acceptance: " << accrate_Lambda_min << "%" << std::endl;
       
       double accrate_A11_trans = round(100*(min(acctot_A11_trans)/(double)j));
       Rcpp::Rcout << "A11 Acceptance: " << accrate_A11_trans << "%" << std::endl;
       
       double accrate_A22_trans = round(100*(min(acctot_A22_trans)/(double)j));
       Rcpp::Rcout << "A22 Acceptance: " << accrate_A22_trans << "%" << std::endl;
       
       double accrate_phi1_trans = round(100*(min(acctot_phi1_trans)/(double)j));
       Rcpp::Rcout << "phi1 Acceptance: " << accrate_phi1_trans << "%" << std::endl;
       
       double accrate_phi2_trans = round(100*(min(acctot_phi2_trans)/(double)j));
       Rcpp::Rcout << "phi2 Acceptance: " << accrate_phi2_trans << "%" << std::endl;
       
       Rcpp::Rcout << "***********************" << std::endl;
       
     }
    
     if((q > 1) | (p > 2)){
       
       double accrate_Lambda_min = round(100*(min(acctot_Lambda_vec)/(double)j));
       Rcpp::Rcout << "Lambda Acceptance (min): " << accrate_Lambda_min << "%" << std::endl;
       
       double accrate_Lambda_max = round(100*(max(acctot_Lambda_vec)/(double)j));
       Rcpp::Rcout << "Lambda Acceptance (max): " << accrate_Lambda_max << "%" << std::endl;
       
       double accrate_A11_trans = round(100*(min(acctot_A11_trans)/(double)j));
       Rcpp::Rcout << "A11 Acceptance: " << accrate_A11_trans << "%" << std::endl;
       
       double accrate_A22_trans = round(100*(min(acctot_A22_trans)/(double)j));
       Rcpp::Rcout << "A22 Acceptance: " << accrate_A22_trans << "%" << std::endl;
       
       double accrate_phi1_trans = round(100*(min(acctot_phi1_trans)/(double)j));
       Rcpp::Rcout << "phi1 Acceptance: " << accrate_phi1_trans << "%" << std::endl;
       
       double accrate_phi2_trans = round(100*(min(acctot_phi2_trans)/(double)j));
       Rcpp::Rcout << "phi2 Acceptance: " << accrate_phi2_trans << "%" << std::endl;
      
       Rcpp::Rcout << "*****************************" << std::endl;
      
       }
    
     }
  
   }
                                  
return Rcpp::List::create(Rcpp::Named("beta") = beta,
                          Rcpp::Named("Lambda") = Lambda,
                          Rcpp::Named("delta") = delta,
                          Rcpp::Named("w1") = w1,
                          Rcpp::Named("w2") = w2,
                          Rcpp::Named("A11") = A11,
                          Rcpp::Named("A22") = A22,
                          Rcpp::Named("A21") = A21,
                          Rcpp::Named("phi1") = phi1,
                          Rcpp::Named("phi2") = phi2,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_Lambda") = acctot_Lambda,
                          Rcpp::Named("acctot_A11_trans") = acctot_A11_trans,
                          Rcpp::Named("acctot_A22_trans") = acctot_A22_trans,
                          Rcpp::Named("acctot_phi1_trans") = acctot_phi1_trans,
                          Rcpp::Named("acctot_phi2_trans") = acctot_phi2_trans);

}
