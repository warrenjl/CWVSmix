#ifndef __CWVSmix__
#define __CWVSmix__

double rnorm_trunc(double mu, 
                   double sigma, 
                   double lower, 
                   double upper);

double norm_rs(double a, 
               double b);

double half_norm_rs(double a, 
                    double b);

double unif_rs(double a, 
               double b);

double exp_rs(double a, 
              double b);

arma::vec rcpp_pgdraw(arma::vec b, 
                      arma::vec c);

Rcpp::NumericVector sampleRcpp(Rcpp::NumericVector x,
                               int size,
                               bool replace,
                               Rcpp::NumericVector prob = Rcpp::NumericVector::create());

Rcpp::List temporal_corr_fun(int p_z,
                             double phi);

double sigma2_epsilon_update(int n,
                             int p,
                             int m,
                             arma::vec y,
                             arma::mat x,
                             arma::mat z, 
                             int likelihood_indicator,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon,
                             arma::vec beta_old,
                             arma::vec eta_full,
                             arma::mat risk_sum);

Rcpp::List w_update(int n,
                    int p,
                    int m,
                    arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec off_set,
                    arma::vec tri_als,
                    int likelihood_indicator,
                    int r,
                    arma::vec beta_old,
                    arma::vec eta_full,
                    arma::mat risk_sum);

arma::vec beta_update(int n,
                      int p,
                      int m,
                      int p_x,
                      arma::mat x, 
                      arma::mat z,
                      arma::vec off_set,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec gamma,
                      arma::vec eta_full,
                      arma::mat risk_sum);

Rcpp::List lambda_update(arma::mat lambda_star,
                         arma::mat lambda_old,
                         int ind,
                         int p,
                         int q,
                         int m,
                         arma::mat x,
                         arma::mat z,
                         int interaction_indicator,
                         arma::vec off_set,
                         arma::vec w,
                         arma::vec gamma,
                         arma::vec beta,
                         arma::vec eta_full,
                         arma::mat risk_sum,
                         arma::mat temporal_corr_inv_old,
                         double metrop_var_lambda_trans,
                         int acctot_lambda);

Rcpp::List rho_update(double rho_old,
                      int p,
                      int m,
                      double alpha_rho,
                      double beta_rho,
                      arma::mat lambda_star,
                      Rcpp::List temporal_corr_info,
                      double metrop_var_rho_trans,
                      int acctot_rho_trans);

Rcpp::List delta_update(arma::vec delta_old,
                        int n,
                        int p,
                        int m,
                        arma::vec y,
                        arma::mat x,
                        arma::mat z,
                        arma::vec off_set,
                        arma::vec w,
                        arma::vec gamma,
                        arma::vec beta,
                        arma::vec w1_old,
                        arma::vec w2_old,
                        double A11_old,
                        double A22_old,
                        double A21_old,
                        arma::vec eta_full,
                        arma::mat risk_sum);

arma::vec delta_star_update(int m,
                            arma::vec delta,
                            arma::vec w1_old,
                            arma::vec w2_old,
                            double A22_old,
                            double A21_old);

Rcpp::List w1_update(int n,
                     int p,
                     int m,
                     arma::mat x,
                     arma::mat z,
                     arma::vec off_set,
                     arma::vec w,
                     arma::vec gamma,
                     arma::vec beta,
                     arma::vec delta,
                     arma::vec delta_star,
                     arma::vec w2_old,
                     double A11_old,
                     double A22_old,
                     double A21_old,
                     arma::mat risk_sum,
                     arma::mat corr_inv1);

arma::vec w2_update(int p,
                    int m,
                    arma::mat z,
                    arma::vec delta_star,
                    arma::vec w1,
                    double A22_old,
                    double A21_old,
                    arma::mat corr_inv2);

Rcpp::List A11_update(double A11_old,
                      int n,
                      int p,
                      int m,
                      arma::mat x,
                      arma::mat z,
                      arma::vec off_set,
                      double sigma2_A,
                      arma::vec w,
                      arma::vec gamma,
                      arma::vec beta,
                      arma::vec delta,
                      arma::vec w1,
                      arma::mat risk_sum,
                      double metrop_var_A11_trans,
                      int acctot_A11_trans);

Rcpp::List A22_update(double A22_old,
                      int m,
                      double sigma2_A,
                      arma::vec delta_star,
                      arma::vec w1,
                      arma::vec w2,
                      double A21_old,
                      double metrop_var_A22_trans,
                      int acctot_A22_trans);

double A21_update(double sigma2_A,
                  arma::vec delta_star,
                  arma::vec w1,
                  arma::vec w2,
                  double A22);

Rcpp::List phi_update(double phi_old,
                      int m,
                      double alpha_phi,
                      double beta_phi,
                      arma::vec w,
                      Rcpp::List temporal_corr_info,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans);

int r_update(int n,
             int p,
             int m,
             arma::vec y,
             arma::mat x,
             arma::mat z,
             arma::vec off_set,
             int a_r,
             int b_r,
             arma::vec beta,
             arma::vec eta_full,
             arma::mat risk_sum);

double neg_two_loglike_update(int n,
                              int p,
                              int m,
                              arma::vec y,
                              arma::mat x,
                              arma::mat z, 
                              arma::vec off_set,
                              arma::vec tri_als,
                              int likelihood_indicator,
                              int r,
                              double sigma2_epsilon,
                              arma::vec beta,
                              arma::vec eta_full,
                              arma::mat risk_sum);

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
                   Rcpp::Nullable<Rcpp::NumericVector> offset,
                   Rcpp::Nullable<Rcpp::NumericVector> tri_als,
                   Rcpp::Nullable<double> a_r_prior,
                   Rcpp::Nullable<double> b_r_prior,
                   Rcpp::Nullable<double> a_sigma2_epsilon_prior,
                   Rcpp::Nullable<double> b_sigma2_epsilon_prior,
                   Rcpp::Nullable<double> sigma2_beta_prior,
                   Rcpp::Nullable<double> alpha_rho_prior,
                   Rcpp::Nullable<double> beta_rho_prior,
                   Rcpp::Nullable<double> sigma2_A_prior,
                   Rcpp::Nullable<double> alpha_phi1_prior,
                   Rcpp::Nullable<double> beta_phi1_prior,
                   Rcpp::Nullable<double> alpha_phi2_prior,
                   Rcpp::Nullable<double> beta_phi2_prior,
                   Rcpp::Nullable<double> r_init,
                   Rcpp::Nullable<double> sigma2_epsilon_init,
                   Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                   Rcpp::Nullable<double> rho_init,
                   Rcpp::Nullable<Rcpp::NumericVector> delta_init,
                   Rcpp::Nullable<Rcpp::NumericVector> w1_init,
                   Rcpp::Nullable<Rcpp::NumericVector> w2_init,
                   Rcpp::Nullable<double> A11_init,
                   Rcpp::Nullable<double> A22_init,
                   Rcpp::Nullable<double> A21_init,
                   Rcpp::Nullable<double> phi1_init,
                   Rcpp::Nullable<double> phi2_init); 

#endif // __CWVSmix__