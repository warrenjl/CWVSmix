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
                   arma::vec mh_scale_Lambda,
                   arma::vec metrop_var_A11_trans,
                   arma::vec metrop_var_A22_trans,
                   arma::vec metrop_var_phi1_trans,
                   arma::vec metrop_var_phi2_trans,
                   Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                   Rcpp::Nullable<double> alpha_Lambda_prior = R_NilValue,
                   Rcpp::Nullable<double> sigma2_A_prior = R_NilValue,
                   Rcpp::Nullable<double> alpha_phi1_prior = R_NilValue,
                   Rcpp::Nullable<double> beta_phi1_prior = R_NilValue,
                   Rcpp::Nullable<double> alpha_phi2_prior = R_NilValue,
                   Rcpp::Nullable<double> beta_phi2_prior = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> Lambda_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> delta_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> w1_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> w2_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> A11_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> A22_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> A21_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> phi1_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> phi2_init = R_NilValue){

//Defining Parameters and Quantities of Interest
int p_x = x.n_cols;
int m = z.n_cols/p;
double max_time = (m - 1);

arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
Rcpp::List Lambda(mcmc_samples);
Rcpp::List delta(mcmc_samples);
Rcpp::List w1(mcmc_samples);
Rcpp::List w2(mcmc_samples);
Rcpp::List eta(mcmc_samples);
for(int j = 0; j < mcmc_samples; ++ j){
  
   arma::mat Lambda_temp(p, q); Lambda_temp.fill(0.00);
   Lambda[j] = Lambda_temp;
   
   arma::mat delta_temp(m, q); delta_temp.fill(0.00);
   delta[j] = delta_temp;
   
   arma::mat w1_temp(m, q); w1_temp.fill(0.00);
   w1[j] = w1_temp;
   
   arma::mat w2_temp(m, q); w2_temp.fill(0.00);
   w2[j] = w2_temp;
   
   arma::mat eta_temp(m, q); eta_temp.fill(0.00);
   eta[j] = eta_temp;
  
   }
arma::mat A11(q, mcmc_samples); A11.fill(0.00);
arma::mat A22(q, mcmc_samples); A22.fill(0.00);
arma::mat A21(q, mcmc_samples); A21.fill(0.00);
arma::mat phi1(q, mcmc_samples); phi1.fill(0.00);
arma::mat phi2(q, mcmc_samples); phi2.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Prior Information
double sigma2_beta = 10000.00;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }

double alpha_Lambda = 1.00;
if(alpha_Lambda_prior.isNotNull()){
  alpha_Lambda = Rcpp::as<double>(alpha_Lambda_prior);
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
     Lambda_temp(k, j) = (1.00/(p - j));
     }
   }
if(Lambda_init.isNotNull()){
  Lambda_temp = Rcpp::as<arma::mat>(Lambda_init);
  }
Lambda[0] = Lambda_temp;

arma::mat delta_temp(m, q); delta_temp.fill(0.00); delta_temp.col(0).fill(1.00);
if(delta_init.isNotNull()){
  delta_temp = Rcpp::as<arma::mat>(delta_init);
  }
delta[0] = delta_temp;

arma::mat w1_temp(m, q); w1_temp.fill(0.00);
if(w1_init.isNotNull()){
  w1_temp = Rcpp::as<arma::mat>(w1_init);
  }
w1[0] = w1_temp;

arma::mat w2_temp(m, q); w2_temp.fill(0.00);
if(w2_init.isNotNull()){
  w2_temp = Rcpp::as<arma::mat>(w2_init);
  }
w2[0] = w2_temp;

A11.col(0).fill(1.00);
if(A11_init.isNotNull()){
  A11.col(0) = Rcpp::as<arma::vec>(A11_init);
  }

A22.col(0).fill(1.00);
if(A22_init.isNotNull()){
  A22.col(0) = Rcpp::as<arma::vec>(A22_init);
  }

A21.col(0).fill(0.00);
if(A21_init.isNotNull()){
  A21.col(0) = Rcpp::as<arma::vec>(A21_init);
  }

phi1.col(0).fill(-log(0.05)/max_time);  //Effective range equal to largest temporal distance in dataset (strong temporal correlation)
if(phi1_init.isNotNull()){
  phi1.col(0).fill(Rcpp::as<double>(phi1_init));
  }

phi2.col(0).fill(-log(0.05)/max_time);  //Effective range equal to largest temporal distance in dataset (strong temporal correlation)
if(phi2_init.isNotNull()){
  phi2.col(0).fill(Rcpp::as<double>(phi2_init));
  }

Rcpp::List temporal_corr_info1(q); 
Rcpp::List temporal_corr_info2(q); 
for(int j = 0; j < q; ++ j){
  
   temporal_corr_info1[j] = temporal_corr_fun(m, 
                                              phi1(j, 0));
  
   temporal_corr_info2[j] = temporal_corr_fun(m, 
                                              phi2(j, 0));
  
   }

delta_temp = Rcpp::as<arma::mat>(delta[0]);
w1_temp = Rcpp::as<arma::mat>(w1[0]);
arma::vec delta_diag(m*q); delta_diag.fill(0.00); 
arma::vec w1_full(m*q); w1_full.fill(0.00);
arma::vec A11_diag(m*q); A11_diag.fill(0.00);
for(int j = 0; j < m; ++ j){
  
   delta_diag.subvec((j*q), (q*(j + 1) - 1)) = trans(delta_temp.row(j));
   w1_full.subvec((j*q), (q*(j + 1) - 1)) = trans(w1_temp.row(j));
   A11_diag.subvec((j*q), (q*(j + 1) - 1)) = A11.col(0);
  
   }
arma::vec eta_full = (delta_diag%A11_diag%w1_full);

neg_two_loglike(0) = neg_two_loglike_update(p,
                                            q,
                                            y,
                                            x,
                                            z, 
                                            beta.col(0),
                                            Lambda[0],
                                            eta_full);

arma::mat eta_temp(m, q); eta_temp.fill(0.00);
for(int j = 0; j < q; ++ j){
  
   arma::vec subset = regspace(j, q, ((m*q) - 1));
   arma::uvec subset_final = conv_to<arma::uvec>::from(subset);
   eta_temp.col(j) = eta_full.elem(subset_final);
  
 }
eta[0] = eta_temp;

//Metropolis Settings
arma::vec acctot_Lambda(q); acctot_Lambda.fill(0);
arma::vec acctot_A11_trans(q); acctot_A11_trans.fill(0);
arma::vec acctot_A22_trans(q); acctot_A22_trans.fill(0);
arma::vec acctot_phi1_trans(q); acctot_phi1_trans.fill(0);
arma::vec acctot_phi2_trans(q); acctot_phi2_trans.fill(0);

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++ j){
  
   //w Update
   Rcpp::List w_output = w_update(p,
                                  q,
                                  y,
                                  x,
                                  z,
                                  beta.col(j-1),
                                  Lambda[j-1],
                                  eta_full);
  
   arma::vec w = w_output[0];
   arma::vec gamma = w_output[1];
  
   //beta Update
   beta.col(j) = beta_update(p,
                             q,
                             x, 
                             z,
                             sigma2_beta,
                             w,
                             gamma,
                             Lambda[j-1],
                             eta_full);
   
   //Lambda Update
   arma::mat Lambda_temp = Lambda[j-1];
   for(int k = 0; k < q; ++ k){
     
      Rcpp::List Lambda_output = Lambda_update(Lambda_temp,
                                               k,
                                               p,
                                               q,
                                               x,
                                               z,
                                               alpha_Lambda,
                                               w,
                                               gamma,
                                               beta.col(j),
                                               eta_full,
                                               mh_scale_Lambda(k),
                                               acctot_Lambda(k));
     
      Lambda_temp.col(k) = Rcpp::as<arma::vec>(Lambda_output[0]);
      acctot_Lambda(k) = Rcpp::as<int>(Lambda_output[1]);
     
      }
   Lambda[j] = Lambda_temp;
   
   //delta Update
   Rcpp::List delta_output = delta_update(delta[j-1],
                                          p,
                                          q,
                                          y,
                                          x,
                                          z,
                                          w,
                                          gamma,
                                          beta.col(j),
                                          Lambda[j],
                                          w1[j-1],
                                          w2[j-1],
                                          A11.col(j-1),
                                          A22.col(j-1),
                                          A21.col(j-1));
        
   delta[j] = Rcpp::as<arma::mat>(delta_output[0]);
   eta_full = Rcpp::as<arma::vec>(delta_output[1]);
   
   //delta_star Update
   arma::mat delta_star = delta_star_update(delta[j],
                                            w1[j-1],
                                            w2[j-1],
                                            A22.col(j-1),
                                            A21.col(j-1));
   
   //w1 Update
   Rcpp::List w1_output = w1_update(p,
                                    q,
                                    x,
                                    z,
                                    w,
                                    gamma,
                                    beta.col(j),
                                    Lambda[j],
                                    delta[j],
                                    delta_star,
                                    w2[j-1],
                                    A11.col(j-1),
                                    A22.col(j-1),
                                    A21.col(j-1),
                                    temporal_corr_info1);
   
   w1[j] = w1_output[0];
   eta_full = Rcpp::as<arma::vec>(w1_output[1]);
   
   //w2 Update
   arma::mat w2_temp = w2[j-1];
   arma::mat w1_temp = Rcpp::as<arma::mat>(w1[j]);
   for(int k = 0; k < q; ++ k){
      
      Rcpp::List temporal_corr_info2_temp = temporal_corr_info2[k];
      w2_temp.col(k) = w2_update(p,
                                 z,
                                 delta_star.col(k),
                                 w1_temp.col(k),
                                 A22(k, (j-1)),
                                 A21(k, (j-1)),
                                 temporal_corr_info2_temp[0]);
      
      }
   w2[j] = w2_temp;
   
   //A11 Update
   Rcpp::List A11_output = A11_update(A11.col(j-1),
                                      p,
                                      q,
                                      x,
                                      z,
                                      sigma2_A,
                                      w,
                                      gamma,
                                      beta.col(j),
                                      Lambda[j],
                                      delta[j],
                                      w1[j],
                                      metrop_var_A11_trans,
                                      acctot_A11_trans);
   
   A11.col(j) = Rcpp::as<arma::vec>(A11_output[0]);
   eta_full = Rcpp::as<arma::vec>(A11_output[1]);
   acctot_A11_trans = Rcpp::as<arma::vec>(A11_output[2]);
   
   //A22 Update
   w1_temp = Rcpp::as<arma::mat>(w1[j]);
   w2_temp = Rcpp::as<arma::mat>(w2[j]);
   for(int k = 0; k < q; ++ k){
   
      Rcpp::List A22_output = A22_update(A22(k, (j-1)),
                                         sigma2_A,
                                         delta_star.col(k),
                                         w1_temp.col(k),
                                         w2_temp.col(k),
                                         A21(k, (j-1)),
                                         metrop_var_A22_trans(k),
                                         acctot_A22_trans(k));
   
      A22(k,j) = Rcpp::as<double>(A22_output[0]);
      acctot_A22_trans(k) = Rcpp::as<double>(A22_output[1]);
      
      }
   
   //A21 Update
   w1_temp = Rcpp::as<arma::mat>(w1[j]);
   w2_temp = Rcpp::as<arma::mat>(w2[j]);
   for(int k = 0; k < q; ++ k){
     
      A21(k,j) = A21_update(sigma2_A,
                            delta_star.col(k),
                            w1_temp.col(k),
                            w2_temp.col(k),
                            A22(k,j));
      
      }
   
   //phi1 Update
   w1_temp = Rcpp::as<arma::mat>(w1[j]);
   for(int k = 0; k < q; ++ k){
     
      Rcpp::List phi1_output = phi_update(phi1(k, (j-1)),
                                          alpha_phi1,
                                          beta_phi1,
                                          w1_temp.col(k),
                                          temporal_corr_info1[k],
                                          metrop_var_phi1_trans(k),
                                          acctot_phi1_trans(k));
     
      phi1(k,j) = Rcpp::as<double>(phi1_output[0]);
      acctot_phi1_trans(k) = Rcpp::as<int>(phi1_output[1]);
      temporal_corr_info1[k] = phi1_output[2];
     
      }
   
   //phi2 Update
   w2_temp = Rcpp::as<arma::mat>(w2[j]);
   for(int k = 0; k < q; ++ k){
     
      arma::mat w2_temp = w2[j];
      Rcpp::List phi2_output = phi_update(phi2(k, (j-1)),
                                          alpha_phi2,
                                          beta_phi2,
                                          w2_temp.col(k),
                                          temporal_corr_info2[k],
                                          metrop_var_phi2_trans(k),
                                          acctot_phi2_trans(k));
     
      phi2(k,j) = Rcpp::as<double>(phi2_output[0]);
      acctot_phi2_trans(k) = Rcpp::as<int>(phi2_output[1]);
      temporal_corr_info2[k] = phi2_output[2];
     
      }
  
   //neg_two_loglike Update
   neg_two_loglike(j) = neg_two_loglike_update(p,
                                               q,
                                               y,
                                               x,
                                               z, 
                                               beta.col(j),
                                               Lambda[j],
                                               eta_full);
   
   //eta Update
   arma::mat eta_temp(m, q); eta_temp.fill(0.00);
   for(int j = 0; j < q; ++ j){
     
     arma::vec subset = regspace(j, q, ((m*q) - 1));
     arma::uvec subset_final = conv_to<arma::uvec>::from(subset);
     eta_temp.col(j) = eta_full.elem(subset_final);
     
     }
   eta[j] = eta_temp;
   
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
  
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
     
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
    
     if(q == 1){
       
       double accrate_Lambda = round(100*(min(acctot_Lambda)/(double)j));
       Rcpp::Rcout << "Lambda Acceptance: " << accrate_Lambda << "%" << std::endl;
       
       double accrate_A11_trans = round(100*(min(acctot_A11_trans)/(double)j));
       Rcpp::Rcout << "A11 Acceptance: " << accrate_A11_trans << "%" << std::endl;
       
       double accrate_A22_trans = round(100*(min(acctot_A22_trans)/(double)j));
       Rcpp::Rcout << "A22 Acceptance: " << accrate_A22_trans << "%" << std::endl;
       
       double accrate_phi1_trans = round(100*(min(acctot_phi1_trans)/(double)j));
       Rcpp::Rcout << "phi1 Acceptance: " << accrate_phi1_trans << "%" << std::endl;
       
       double accrate_phi2_trans = round(100*(min(acctot_phi2_trans)/(double)j));
       Rcpp::Rcout << "phi2 Acceptance: " << accrate_phi2_trans << "%" << std::endl;
      
       Rcpp::Rcout << "**********************" << std::endl;
      
       }
    
     if(q > 1){
       
       double accrate_Lambda_min = round(100*(min(acctot_Lambda)/(double)j));
       Rcpp::Rcout << "Lambda Acceptance (min): " << accrate_Lambda_min << "%" << std::endl;
       
       double accrate_Lambda_max = round(100*(max(acctot_Lambda)/(double)j));
       Rcpp::Rcout << "Lambda Acceptance (max): " << accrate_Lambda_max << "%" << std::endl;
       
       double accrate_A11_trans_min = round(100*(min(acctot_A11_trans)/(double)j));
       Rcpp::Rcout << "A11 Acceptance (min): " << accrate_A11_trans_min << "%" << std::endl;
       
       double accrate_A11_trans_max = round(100*(max(acctot_A11_trans)/(double)j));
       Rcpp::Rcout << "A11 Acceptance (max): " << accrate_A11_trans_max << "%" << std::endl;
       
       double accrate_A22_trans_min = round(100*(min(acctot_A22_trans)/(double)j));
       Rcpp::Rcout << "A22 Acceptance (min): " << accrate_A22_trans_min << "%" << std::endl;
       
       double accrate_A22_trans_max = round(100*(max(acctot_A22_trans)/(double)j));
       Rcpp::Rcout << "A22 Acceptance (max): " << accrate_A22_trans_max << "%" << std::endl;

       double accrate_phi1_trans_min = round(100*(min(acctot_phi1_trans)/(double)j));
       Rcpp::Rcout << "phi1 Acceptance (min): " << accrate_phi1_trans_min << "%" << std::endl;
    
       double accrate_phi1_trans_max = round(100*(max(acctot_phi1_trans)/(double)j));
       Rcpp::Rcout << "phi1 Acceptance (max): " << accrate_phi1_trans_max << "%" << std::endl;
       
       double accrate_phi2_trans_min = round(100*(min(acctot_phi2_trans)/(double)j));
       Rcpp::Rcout << "phi2 Acceptance (min): " << accrate_phi2_trans_min << "%" << std::endl;
       
       double accrate_phi2_trans_max = round(100*(max(acctot_phi2_trans)/(double)j));
       Rcpp::Rcout << "phi2 Acceptance (max): " << accrate_phi2_trans_max << "%" << std::endl;
      
       Rcpp::Rcout << "****************************" << std::endl;
      
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
                          Rcpp::Named("acctot_phi2_trans") = acctot_phi2_trans,
                          Rcpp::Named("eta") = eta);

}
