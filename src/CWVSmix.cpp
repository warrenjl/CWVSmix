#include "RcppArmadillo.h"
#include "CWVSmix.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List CWVSmix(int mcmc_samples,
                   int p,
                   arma::vec y,
                   arma::mat x,
                   arma::mat z,
                   arma::vec metrop_var_lambda_trans,
                   double metrop_var_rho_trans,
                   double metrop_var_A11_trans,
                   double metrop_var_A22_trans,
                   double metrop_var_phi1_trans,
                   double metrop_var_phi2_trans,
                   int interaction_indicator,
                   int likelihood_indicator,
                   Rcpp::Nullable<Rcpp::NumericVector> offset = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> trials = R_NilValue,
                   Rcpp::Nullable<double> a_r_prior = R_NilValue,
                   Rcpp::Nullable<double> b_r_prior = R_NilValue,
                   Rcpp::Nullable<double> a_sigma2_epsilon_prior = R_NilValue,
                   Rcpp::Nullable<double> b_sigma2_epsilon_prior = R_NilValue,
                   Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                   Rcpp::Nullable<double> alpha_rho_prior = R_NilValue,
                   Rcpp::Nullable<double> beta_rho_prior = R_NilValue,
                   Rcpp::Nullable<double> sigma2_A_prior = R_NilValue,
                   Rcpp::Nullable<double> alpha_phi1_prior = R_NilValue,
                   Rcpp::Nullable<double> beta_phi1_prior = R_NilValue,
                   Rcpp::Nullable<double> alpha_phi2_prior = R_NilValue,
                   Rcpp::Nullable<double> beta_phi2_prior = R_NilValue,
                   Rcpp::Nullable<double> r_init = R_NilValue,
                   Rcpp::Nullable<double> sigma2_epsilon_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                   Rcpp::Nullable<double> rho_init = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericVector> delta_init = R_NilValue,
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
int q = 0;
if(interaction_indicator == 1){
  q = (-1 + sqrt(1 + 8*p))/2;  //Solving for the Number of Main Effect Pollutants
  } 
   
arma::vec r(mcmc_samples); r.fill(0.00);
arma::vec sigma2_epsilon(mcmc_samples); sigma2_epsilon.fill(0.00);
arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
Rcpp::List lambda(mcmc_samples);
for(int j = 0; j < mcmc_samples; ++ j){
      
   arma::mat lambda_temp(p, m); lambda_temp.fill(0.00);
   lambda[j] = lambda_temp;
      
   }
arma::vec rho(mcmc_samples); rho.fill(0.00);
arma::mat delta(m, mcmc_samples);
arma::mat w1(m, mcmc_samples); w1.fill(0.00);
arma::mat w2(m, mcmc_samples); w2.fill(0.00);
arma::vec A11(mcmc_samples); A11.fill(0.00);
arma::vec A22(mcmc_samples); A22.fill(0.00);
arma::vec A21(mcmc_samples); A21.fill(0.00);
arma::vec phi1(mcmc_samples); phi1.fill(0.00);
arma::vec phi2(mcmc_samples); phi2.fill(0.00);
arma::mat alpha(m, mcmc_samples); alpha.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);
   
arma::vec off_set(n); off_set.fill(0.00);
if(offset.isNotNull()){
  off_set = Rcpp::as<arma::vec>(offset);
  }

arma::vec tri_als(n); tri_als.fill(1);
if(trials.isNotNull()){
  tri_als = Rcpp::as<arma::vec>(trials);
  }

//Prior Information
int a_r = 1;
if(a_r_prior.isNotNull()){
  a_r = Rcpp::as<int>(a_r_prior);
  }

int b_r = 100;
if(b_r_prior.isNotNull()){
  b_r = Rcpp::as<int>(b_r_prior);
  }

double a_sigma2_epsilon = 0.01;
if(a_sigma2_epsilon_prior.isNotNull()){
  a_sigma2_epsilon = Rcpp::as<double>(a_sigma2_epsilon_prior);
  }
   
double b_sigma2_epsilon = 0.01;
if(b_sigma2_epsilon_prior.isNotNull()){
  b_sigma2_epsilon = Rcpp::as<double>(b_sigma2_epsilon_prior);
  }
   
double sigma2_beta = 10000.00;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }
   
double alpha_rho = 1.00;  
if(alpha_rho_prior.isNotNull()){
  alpha_rho = Rcpp::as<double>(alpha_rho_prior);
  }
   
double beta_rho = 1.00;
if(beta_rho_prior.isNotNull()){
  beta_rho = Rcpp::as<double>(beta_rho_prior);
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
r(0) = b_r;
if(r_init.isNotNull()){
  r(0) = Rcpp::as<int>(r_init);
  }

sigma2_epsilon(0) = 1.00;
if(sigma2_epsilon_init.isNotNull()){
  sigma2_epsilon(0) = Rcpp::as<double>(sigma2_epsilon_init);
  }
   
beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }
   
arma::mat lambda_star = arma::randn(p, m);
arma::mat lambda_star_max(p, m); lambda_star_max.fill(0.00);
   
if(interaction_indicator == 0){
  for(int j = 0; j < p; ++ j){
     for(int k = 0; k < m; ++ k){
        lambda_star_max(j, k) = (lambda_star(j, k) > 0.00)*lambda_star(j, k);
        }   
     }
  }
   
if(interaction_indicator == 1){
      
  for(int j = 0; j < q; ++ j){
     for(int k = 0; k < m; ++ k){
        lambda_star_max(j, k) = (lambda_star(j, k) > 0.00)*lambda_star(j, k);
        }   
     }
      
  int counter = q;
  for(int j = 0; j < (q - 1); ++ j){
     for(int k = (j + 1); k < q; ++ k){
            
        for(int l = 0; l < m; ++ l){
           lambda_star_max(counter, l) = (lambda_star(j, l) > 0.00)*(lambda_star(k, l) > 0.00)*(lambda_star(counter, l) > 0.00)*lambda_star(counter, l);
           }
            
        counter = counter + 
                  1;
            
        }
     }
      
  }
   
arma::mat lambda_temp(p, m); lambda_temp.fill(0.00);
for(int j = 0; j < m; ++ j){
      
   if(sum(lambda_star_max.col(j)) > 0.00){
     for(int k = 0; k < p; ++ k){
        lambda_temp(k, j) = lambda_star_max(k, j)/sum(lambda_star_max.col(j));
        }
     }
      
   if(sum(lambda_star_max.col(j)) == 0.00){
     lambda_temp.col(j).fill(1.00/p);
     }
      
  }
lambda[0] = lambda_temp;
   
rho(0) = -log(0.05)/max_time;  //Effective range equal to largest temporal distance in dataset (strong temporal correlation)
if(rho_init.isNotNull()){
  rho(0) = Rcpp::as<double>(rho_init);
  }
   
delta.col(0).fill(1.00);
if(delta_init.isNotNull()){
  delta.col(0) = Rcpp::as<arma::vec>(delta_init);
  }
   
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
   
Rcpp::List temporal_corr_info = temporal_corr_fun(m, 
                                                  rho(0));
  
Rcpp::List temporal_corr_info1 = temporal_corr_fun(m, 
                                                   phi1(0));
   
Rcpp::List temporal_corr_info2 = temporal_corr_fun(m, 
                                                   phi2(0));
   
arma::vec eta_full = A11(0)*(w1.col(0)%delta.col(0));
   
arma::mat risk_sum(n, m); risk_sum.fill(0.00);
for(int j = 0; j < m; ++ j){
   risk_sum.col(j) = z.cols((p*j), (p*(j + 1) - 1))*Rcpp::as<arma::mat>(lambda[0]).col(j);
   }
   
neg_two_loglike(0) = neg_two_loglike_update(n,
                                            p,
                                            m,
                                            y,
                                            x,
                                            z,
                                            off_set,
                                            tri_als,
                                            likelihood_indicator,
                                            r(0),
                                            sigma2_epsilon(0),
                                            beta.col(0),
                                            eta_full,
                                            risk_sum);
   
//Metropolis Settings
arma::vec acctot_lambda(m); acctot_lambda.fill(0);
int acctot_rho_trans = 0;
int acctot_A11_trans = 0;
int acctot_A22_trans = 0; 
int acctot_phi1_trans = 0;
int acctot_phi2_trans = 0; 
   
//Main Sampling Loop
arma::vec w(y.size()); w.fill(0.00);
arma::vec gamma = y;
if(likelihood_indicator == 2){
  
  Rcpp::List w_output = w_update(n,
                                 p,
                                 m,
                                 y,
                                 x,
                                 z,
                                 off_set,
                                 tri_als,
                                 likelihood_indicator,
                                 r(0),
                                 beta.col(0),
                                 eta_full,
                                 risk_sum);
  w = Rcpp::as<arma::vec>(w_output[0]);
  gamma = Rcpp::as<arma::vec>(w_output[1]);
  
  }

for(int j = 1; j < mcmc_samples; ++ j){
      
   if(likelihood_indicator == 1){
         
     //sigma2_epsilon Update
     sigma2_epsilon(j) = sigma2_epsilon_update(n,
                                               p,
                                               m,
                                               y,
                                               x,
                                               z, 
                                               likelihood_indicator,
                                               a_sigma2_epsilon,
                                               b_sigma2_epsilon,
                                               beta.col(j-1),
                                               eta_full,
                                               risk_sum);
     w.fill(1.00/sigma2_epsilon(j));
         
     }
      
   if(likelihood_indicator == 0){
         
     //w Update
     Rcpp::List w_output = w_update(n,
                                    p,
                                    m,
                                    y,
                                    x,
                                    z,
                                    off_set,
                                    tri_als,
                                    likelihood_indicator,
                                    r(j-1),
                                    beta.col(j-1),
                                    eta_full,
                                    risk_sum);
         
     w = Rcpp::as<arma::vec>(w_output[0]);
     gamma = Rcpp::as<arma::vec>(w_output[1]);
         
     }
      
   //beta Update
   beta.col(j) = beta_update(n,
                             p,
                             m,
                             p_x,
                             x, 
                             z,
                             off_set,
                             sigma2_beta,
                             w,
                             gamma,
                             eta_full,
                             risk_sum);
      
   //lambda Update
   arma::mat lambda_temp = lambda[j-1];
   for(int k = 0; k < m; ++ k){
         
      Rcpp::List lambda_output = lambda_update(lambda_star,
                                               lambda_temp,
                                               k,
                                               p,
                                               q,
                                               m,
                                               x,
                                               z,
                                               interaction_indicator,
                                               off_set,
                                               w,
                                               gamma,
                                               beta.col(j),
                                               eta_full,
                                               risk_sum,
                                               temporal_corr_info[0],
                                               metrop_var_lambda_trans(k),
                                               acctot_lambda(k));
         
      risk_sum.col(k) = Rcpp::as<arma::vec>(lambda_output[0]);
      lambda_star.col(k) = Rcpp::as<arma::vec>(lambda_output[1]);
      lambda_temp.col(k) = Rcpp::as<arma::vec>(lambda_output[2]);
      acctot_lambda(k) = Rcpp::as<int>(lambda_output[3]);
         
      }
   lambda[j] = lambda_temp;
      
   //rho Update
   Rcpp::List rho_output = rho_update(rho(j-1),
                                      p,
                                      m,
                                      alpha_rho,
                                      beta_rho,
                                      lambda_star,
                                      temporal_corr_info,
                                      metrop_var_rho_trans,
                                      acctot_rho_trans);
      
   rho(j) = Rcpp::as<double>(rho_output[0]);
   acctot_rho_trans = Rcpp::as<int>(rho_output[1]);
   temporal_corr_info = rho_output[2];
      
   //delta Update
   Rcpp::List delta_output = delta_update(delta.col(j-1),
                                          n,
                                          p,
                                          m,
                                          y,
                                          x,
                                          z,
                                          off_set,
                                          w,
                                          gamma,
                                          beta.col(j),
                                          w1.col(j-1),
                                          w2.col(j-1),
                                          A11(j-1),
                                          A22(j-1),
                                          A21(j-1),
                                          eta_full,
                                          risk_sum);
      
   delta.col(j) = Rcpp::as<arma::vec>(delta_output[0]);
   eta_full = Rcpp::as<arma::vec>(delta_output[1]);
      
   //delta_star Update
   arma::vec delta_star = delta_star_update(m,
                                            delta.col(j),
                                            w1.col(j-1),
                                            w2.col(j-1),
                                            A22(j-1),
                                            A21(j-1));
      
   //w1 Update
   Rcpp::List w1_output = w1_update(n,
                                    p,
                                    m,
                                    x,
                                    z,
                                    off_set,
                                    w,
                                    gamma,
                                    beta.col(j),
                                    delta.col(j),
                                    delta_star,
                                    w2.col(j-1),
                                    A11(j-1),
                                    A22(j-1),
                                    A21(j-1),
                                    risk_sum,
                                    temporal_corr_info1[0]);
      
   w1.col(j) = Rcpp::as<arma::vec>(w1_output[0]);
   eta_full = Rcpp::as<arma::vec>(w1_output[1]);
      
   //w2 Update
   w2.col(j) = w2_update(p,
                         m,
                         z,
                         delta_star,
                         w1.col(j),
                         A22(j-1),
                         A21(j-1),
                         temporal_corr_info2[0]);
      
   //A11 Update
   Rcpp::List A11_output = A11_update(A11(j-1),
                                      n,
                                      p,
                                      m,
                                      x,
                                      z,
                                      off_set,
                                      sigma2_A,
                                      w,
                                      gamma,
                                      beta.col(j),
                                      delta.col(j),
                                      w1.col(j),
                                      risk_sum,
                                      metrop_var_A11_trans,
                                      acctot_A11_trans);
      
   A11(j) = Rcpp::as<double>(A11_output[0]);
   eta_full = Rcpp::as<arma::vec>(A11_output[1]);
   acctot_A11_trans = Rcpp::as<double>(A11_output[2]);
      
   //A22 Update
   Rcpp::List A22_output = A22_update(A22(j-1),
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
   A21(j) = A21_update(sigma2_A,
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
      
   //alpha Update
   alpha.col(j) = delta.col(j)%(A11(j)*w1.col(j));
   
   if(likelihood_indicator == 2){
     
     //r Update
     r(j) = r_update(n,
                     p,
                     m,
                     y,
                     x,
                     z,
                     off_set,
                     a_r,
                     b_r,
                     beta.col(j),
                     eta_full,
                     risk_sum);
     
     //w Update
     Rcpp::List w_output = w_update(n,
                                    p,
                                    m,
                                    y,
                                    x,
                                    z,
                                    off_set,
                                    tri_als,
                                    likelihood_indicator,
                                    r(j),
                                    beta.col(j),
                                    eta_full,
                                    risk_sum);
     w = Rcpp::as<arma::vec>(w_output[0]);
     gamma = Rcpp::as<arma::vec>(w_output[1]);
     
     }
      
   //neg_two_loglike Update
   neg_two_loglike(j) = neg_two_loglike_update(n,
                                               p,
                                               m,
                                               y,
                                               x,
                                               z, 
                                               off_set,
                                               tri_als,
                                               likelihood_indicator,
                                               r(j),
                                               sigma2_epsilon(j),
                                               beta.col(j),
                                               eta_full,
                                               risk_sum);
      
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
      
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
         
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
         
     double accrate_lambda_min = round(100*(min(acctot_lambda)/(double)j));
     Rcpp::Rcout << "lambda Acceptance (min): " << accrate_lambda_min << "%" << std::endl;
         
     double accrate_lambda_max = round(100*(max(acctot_lambda)/(double)j));
     Rcpp::Rcout << "lambda Acceptance (max): " << accrate_lambda_max << "%" << std::endl;
         
     double accrate_rho_trans = round(100*(min(acctot_rho_trans)/(double)j));
     Rcpp::Rcout << "rho Acceptance: " << accrate_rho_trans << "%" << std::endl;
         
     double accrate_A11_trans = round(100*(min(acctot_A11_trans)/(double)j));
     Rcpp::Rcout << "A11 Acceptance: " << accrate_A11_trans << "%" << std::endl;
         
     double accrate_A22_trans = round(100*(min(acctot_A22_trans)/(double)j));
     Rcpp::Rcout << "A22 Acceptance: " << accrate_A22_trans << "%" << std::endl;
         
     double accrate_phi1_trans = round(100*(min(acctot_phi1_trans)/(double)j));
     Rcpp::Rcout << "phi1 Acceptance: " << accrate_phi1_trans << "%" << std::endl;
         
     double accrate_phi2_trans = round(100*(min(acctot_phi2_trans)/(double)j));
     Rcpp::Rcout << "phi2 Acceptance: " << accrate_phi2_trans << "%" << std::endl;
         
     Rcpp::Rcout << "****************************" << std::endl;
         
     }
      
   }
   
return Rcpp::List::create(Rcpp::Named("sigma2_epsilon") = sigma2_epsilon,
                          Rcpp::Named("beta") = beta,
                          Rcpp::Named("lambda") = lambda,
                          Rcpp::Named("phi_lambda") = rho,  //Name Change to Reflect Model Details
                          Rcpp::Named("gamma") = delta,  //Name Change to Reflect Model Details
                          Rcpp::Named("delta1") = w1,  //Name Change to Reflect Model Details
                          Rcpp::Named("delta2") = w2,  //Name Change to Reflect Model Details
                          Rcpp::Named("A11") = A11,
                          Rcpp::Named("A22") = A22,
                          Rcpp::Named("A21") = A21,
                          Rcpp::Named("phi1") = phi1,
                          Rcpp::Named("phi2") = phi2,
                          Rcpp::Named("alpha") = alpha,
                          Rcpp::Named("r") = r,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_lambda") = acctot_lambda,
                          Rcpp::Named("acctot_A11_trans") = acctot_A11_trans,
                          Rcpp::Named("acctot_A22_trans") = acctot_A22_trans,
                          Rcpp::Named("acctot_phi1_trans") = acctot_phi1_trans,
                          Rcpp::Named("acctot_phi2_trans") = acctot_phi2_trans);
   
}