#include "RcppArmadillo.h"
#include "CWMix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List CWMix(int mcmc_samples,
                 int p,
                 int q,
                 arma::vec y,
                 arma::mat x,
                 arma::mat z,
                 arma::vec metrop_var_phi_trans,
                 arma::vec metrop_scale_Lambda,
                 Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                 Rcpp::Nullable<double> alpha_beta_sigma2_prior = R_NilValue,
                 Rcpp::Nullable<double> beta_beta_sigma2_prior = R_NilValue,
                 Rcpp::Nullable<double> alpha_beta_phi_prior = R_NilValue,
                 Rcpp::Nullable<double> beta_beta_phi_prior = R_NilValue,
                 Rcpp::Nullable<double> alpha_Lambda_prior = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericMatrix> eta_init = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericVector> sigma2_eta_init = R_NilValue,
                 Rcpp::Nullable<double> beta_sigma2_init = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericVector> phi_init = R_NilValue,
                 Rcpp::Nullable<double> beta_phi_init = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericMatrix> Lambda_init = R_NilValue){

//Defining Parameters and Quantities of Interest
int p_x = x.n_cols;
int m = z.n_cols/p;
double max_time = (m - 1);
arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
Rcpp::List eta(mcmc_samples);
for(int j = 0; j < mcmc_samples; ++j){
   arma::mat eta_temp(q,m); eta_temp.fill(0.00);
   eta[j] = eta_temp;
   }
arma::mat sigma2_eta(q, mcmc_samples); sigma2_eta.fill(0.00);
arma::vec beta_sigma2(mcmc_samples); beta_sigma2.fill(0.00);
arma::mat phi(q, mcmc_samples); phi.fill(0.00);
arma::vec beta_phi(mcmc_samples); beta_phi.fill(0.00);
Rcpp::List Lambda(mcmc_samples);
for(int j = 0; j < mcmc_samples; ++j){
   arma::mat Lambda_temp(p,q); Lambda_temp.fill(0.00);
   Lambda[j] = Lambda_temp;
   }
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Prior Information
double sigma2_beta = 10000.00;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }

double alpha_beta_sigma2 = 1.00;
if(alpha_beta_sigma2_prior.isNotNull()){
  alpha_beta_sigma2 = Rcpp::as<double>(alpha_beta_sigma2_prior);
  }
  
double beta_beta_sigma2 = 1.00;
if(beta_beta_sigma2_prior.isNotNull()){
  beta_beta_sigma2 = Rcpp::as<double>(beta_beta_sigma2_prior);
  }

double alpha_beta_phi = 1.00;
if(alpha_beta_phi_prior.isNotNull()){
  alpha_beta_phi = Rcpp::as<double>(alpha_beta_phi_prior);
  }

double beta_beta_phi = 1.00;
if(beta_beta_phi_prior.isNotNull()){
  beta_beta_phi = Rcpp::as<double>(beta_beta_phi_prior);
  }

double alpha_Lambda = 1.00;
if(alpha_Lambda_prior.isNotNull()){
  alpha_Lambda = Rcpp::as<double>(alpha_Lambda_prior);
  }

//Initial Values
beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

arma::mat eta_temp(q, m); eta_temp.fill(0.00);
if(eta_init.isNotNull()){
  eta_temp = Rcpp::as<arma::mat>(eta_init);
  }
eta[0] = eta_temp;

sigma2_eta.col(0).fill(1.00);
if(sigma2_eta_init.isNotNull()){
  sigma2_eta.col(0) = Rcpp::as<arma::vec>(sigma2_eta_init);
  }

beta_sigma2(0) = 1.00;
if(beta_sigma2_init.isNotNull()){
  beta_sigma2(0) = Rcpp::as<double>(beta_sigma2_init);
  }

phi.col(0).fill(-log(0.05)/max_time);  //Effective range equal to largest temporal distance in dataset (strong temporal correlation)
if(phi_init.isNotNull()){
  phi.col(0) = Rcpp::as<arma::vec>(phi_init);
  }

beta_phi(0) = 1.00;
if(beta_phi_init.isNotNull()){
  beta_phi(0) = Rcpp::as<double>(beta_phi_init);
  }

arma::mat Lambda_temp(p, q); Lambda_temp.fill(0.00);
for(int j = 0; j < q; ++ j){
   for(int k = j; k < p; ++ k){
      Lambda_temp(k,j) = (1.00/(p - j));
      }
   }
if(Lambda_init.isNotNull()){
  Lambda_temp = Rcpp::as<arma::mat>(Lambda_init);
  }
Lambda[0] = Lambda_temp;

Rcpp::List temporal_corr_info(q); 
for(int j = 0; j < q; ++j){
   temporal_corr_info[j] = temporal_corr_fun(m, 
                                             phi(j,0));
   }

neg_two_loglike(0) = neg_two_loglike_update(p,
                                            q,
                                            y,
                                            x,
                                            z, 
                                            beta.col(0),
                                            eta[0],
                                            Lambda[0]);

//Metropolis Settings
arma::vec acctot_phi_trans(q); acctot_phi_trans.fill(0);
arma::vec acctot_Lambda(q); acctot_Lambda.fill(0);

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++ j){
  
  //w Update
  Rcpp::List w_output = w_update(p,
                                 q,
                                 y,
                                 x,
                                 z,
                                 beta.col(j-1),
                                 eta[j-1],
                                 Lambda[j-1]);
  
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
                            eta[j-1],
                            Lambda[j-1]);
   
  //eta Update
  eta[j] = eta_update(p,
                      q,
                      x, 
                      z,
                      w,
                      gamma,
                      beta.col(j),
                      sigma2_eta.col(j-1),
                      temporal_corr_info,
                      Lambda[j-1]);
  
  //sigma2_eta Updates
  arma::mat eta_temp = eta[j];
  arma::mat eta_temp_trans = trans(eta_temp);
  for(int k = 0; k < q; ++ k){
    
     Rcpp::List temporal_corr_info_temp = temporal_corr_info[k];
     sigma2_eta(k,j) = sigma2_eta_update(p,
                                         z,
                                         eta_temp_trans.col(k),
                                         temporal_corr_info_temp[0],
                                         beta_sigma2(j-1));
     
     }
   
  //beta_sigma2 Update
  beta_sigma2(j) = beta_sigma2_update(q,
                                      sigma2_eta.col(j),
                                      alpha_beta_sigma2,
                                      beta_beta_sigma2);
  
  //phi Updates
  for(int k = 0; k < q; ++ k){
     
     Rcpp::List phi_output = phi_update(phi(k, (j-1)),
                                        eta_temp_trans.col(k),
                                        sigma2_eta(k,j),
                                        temporal_corr_info[k],
                                        beta_phi(j-1),
                                        metrop_var_phi_trans(k),
                                        acctot_phi_trans(k));
  
     phi(k,j) = Rcpp::as<double>(phi_output[0]);
     acctot_phi_trans(k) = Rcpp::as<int>(phi_output[1]);
     temporal_corr_info[k] = phi_output[2];
     
     }
  
  //beta_phi Update
  beta_phi(j) = beta_phi_update(q,
                                phi.col(j),
                                alpha_beta_phi,
                                beta_beta_phi);
  
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
                                              eta[j],
                                              metrop_scale_Lambda(k),
                                              acctot_Lambda(k));
    
     Lambda_temp.col(k) = Rcpp::as<arma::vec>(Lambda_output[0]);
     acctot_Lambda(k) = Rcpp::as<int>(Lambda_output[1]);
    
     }
  Lambda[j] = Lambda_temp;
  
  //neg_two_loglike Update
  neg_two_loglike(j) = neg_two_loglike_update(p,
                                              q,
                                              y,
                                              x,
                                              z, 
                                              beta.col(j),
                                              eta[j],
                                              Lambda[j]);
  
  //Progress
  if((j + 1) % 10 == 0){ 
    Rcpp::checkUserInterrupt();
    }
  
  if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
    double completion = round(100*((j + 1)/(double)mcmc_samples));
    Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
    
    if(q == 1){
      
      double accrate_phi_trans = round(100*(min(acctot_phi_trans)/(double)j));
      Rcpp::Rcout << "phi Acceptance: " << accrate_phi_trans << "%" << std::endl;
      
      double accrate_Lambda = round(100*(min(acctot_Lambda)/(double)j));
      Rcpp::Rcout << "Lambda Acceptance: " << accrate_Lambda << "%" << std::endl;
      
      Rcpp::Rcout << "**********************" << std::endl;
      
      }
    
    if(q > 1){

      double accrate_phi_trans_min = round(100*(min(acctot_phi_trans)/(double)j));
      Rcpp::Rcout << "phi Acceptance (min): " << accrate_phi_trans_min << "%" << std::endl;
    
      double accrate_phi_trans_max = round(100*(max(acctot_phi_trans)/(double)j));
      Rcpp::Rcout << "phi Acceptance (max): " << accrate_phi_trans_max << "%" << std::endl;
    
      double accrate_Lambda_min = round(100*(min(acctot_Lambda)/(double)j));
      Rcpp::Rcout << "Lambda Acceptance (min): " << accrate_Lambda_min << "%" << std::endl;
    
      double accrate_Lambda_max = round(100*(max(acctot_Lambda)/(double)j));
      Rcpp::Rcout << "Lambda Acceptance (max): " << accrate_Lambda_max << "%" << std::endl;
      
      Rcpp::Rcout << "****************************" << std::endl;
      
      }
    
    }
  
  }
                                  
return Rcpp::List::create(Rcpp::Named("beta") = beta,
                          Rcpp::Named("eta") = eta,
                          Rcpp::Named("sigma2_eta") = sigma2_eta,
                          Rcpp::Named("beta_sigma2") = beta_sigma2,
                          Rcpp::Named("phi") = phi,
                          Rcpp::Named("beta_phi") = beta_phi,
                          Rcpp::Named("Lambda") = Lambda,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_phi_trans") = acctot_phi_trans,
                          Rcpp::Named("acctot_Lambda") = acctot_Lambda);

}
