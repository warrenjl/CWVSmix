#ifndef __CWMix__
#define __CWMix__

arma::vec rcpp_pgdraw(double b, 
                      arma::vec c);

Rcpp::List temporal_corr_fun(int p_z,
                             double phi);

Rcpp::List w_update(int p,
                    int q,
                    arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta_old,
                    arma::mat eta_old,
                    arma::mat Lambda_old);

arma::vec beta_update(int p,
                      int q,
                      arma::mat x, 
                      arma::mat z,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec gamma,
                      arma::mat eta_old,
                      arma::mat Lambda_old);

arma::mat eta_update(int p,
                     int q,
                     arma::mat x, 
                     arma::mat z,
                     arma::vec w,
                     arma::vec gamma,
                     arma::vec beta,
                     arma::vec sigma2_eta_old,
                     Rcpp::List temporal_corr_info,
                     arma::mat Lambda_old);

double sigma2_eta_update(int p,
                         arma::mat z,
                         arma::vec eta,
                         arma::mat corr_inv,
                         double beta_sigma2_old);

double beta_sigma2_update(int q,
                          arma::vec sigma2_eta,
                          double alpha_beta_sigma2,
                          double beta_beta_sigma2);

Rcpp::List phi_update(double phi_old,
                      arma::vec eta,
                      double sigma2_eta,
                      Rcpp::List temporal_corr_info,
                      double beta_phi_old,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans);

double beta_phi_update(int q,
                       arma::vec phi,
                       double alpha_beta_phi,
                       double beta_beta_phi);

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
                         arma::mat eta,
                         double metrop_scale_Lambda,
                         int acctot_Lambda);

double neg_two_loglike_update(int q,
                              int p,
                              arma::vec y,
                              arma::mat x,
                              arma::mat z, 
                              arma::vec beta,
                              arma::mat eta,
                              arma::mat Lambda);

Rcpp::List CWMix(int mcmc_samples,
                 int p,
                 int q,
                 arma::vec y,
                 arma::mat x,
                 arma::mat z,
                 arma::vec metrop_var_phi_trans,
                 arma::vec metrop_scale_Lambda,
                 Rcpp::Nullable<double> sigma2_beta_prior,
                 Rcpp::Nullable<double> alpha_beta_sigma2_prior,
                 Rcpp::Nullable<double> beta_beta_sigma2_prior,
                 Rcpp::Nullable<double> alpha_beta_phi_prior,
                 Rcpp::Nullable<double> beta_beta_phi_prior,
                 Rcpp::Nullable<double> alpha_Lambda_prior,
                 Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                 Rcpp::Nullable<Rcpp::NumericMatrix> eta_init,
                 Rcpp::Nullable<Rcpp::NumericVector> sigma2_eta_init,
                 Rcpp::Nullable<double> beta_sigma2_init,
                 Rcpp::Nullable<Rcpp::NumericVector> phi_init,
                 Rcpp::Nullable<double> beta_phi_init,
                 Rcpp::Nullable<Rcpp::NumericMatrix> Lambda_init); 

#endif // __CWMix__
