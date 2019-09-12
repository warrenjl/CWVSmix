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

arma::vec rcpp_pgdraw(double b, 
                      arma::vec c);

Rcpp::List temporal_corr_fun(int p_z,
                             double phi);

Rcpp::List w_update(int p,
                    int q,
                    int m,
                    arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec beta_old,
                    arma::mat Lambda_old,
                    arma::vec eta_full);

arma::vec beta_update(int n,
                      int p,
                      int q,
                      int m,
                      int p_x,
                      arma::mat x, 
                      arma::mat z,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec gamma,
                      arma::mat Lambda_old,
                      arma::vec eta_full);

Rcpp::List Lambda_update(arma::vec lambda_star,
                         arma::mat Lambda_old,
                         int ind,
                         int p,
                         int q,
                         int m,
                         arma::mat x,
                         arma::mat z,
                         arma::vec alpha_Lambda,
                         arma::vec w,
                         arma::vec gamma,
                         arma::vec beta,
                         arma::vec eta_full,
                         arma::vec metrop_var_Lambda,
                         arma::vec acctot_Lambda);

Rcpp::List delta_update(int stable_ind,
                        arma::mat delta_old,
                        int p,
                        int q,
                        int m,
                        arma::vec y,
                        arma::mat x,
                        arma::mat z,
                        arma::vec w,
                        arma::vec gamma,
                        arma::vec beta,
                        arma::mat Lambda,
                        arma::vec w1_old,
                        arma::vec w2_old,
                        double A11_old,
                        double A22_old,
                        double A21_old);

arma::mat delta_star_update(int q,
                            int m,
                            arma::mat delta,
                            arma::vec w1_old,
                            arma::vec w2_old,
                            double A22_old,
                            double A21_old);

Rcpp::List w1_update(int n,
                     int p,
                     int q,
                     int m,
                     arma::mat x,
                     arma::mat z,
                     arma::vec w,
                     arma::vec gamma,
                     arma::vec beta,
                     arma::mat Lambda,
                     arma::mat delta,
                     arma::mat delta_star,
                     arma::vec w2_old,
                     double A11_old,
                     double A22_old,
                     double A21_old,
                     arma::mat corr_inv1);

arma::vec w2_update(int p,
                    int q,
                    int m,
                    arma::mat z,
                    arma::mat delta_star,
                    arma::vec w1,
                    double A22_old,
                    double A21_old,
                    arma::mat corr_inv2);

Rcpp::List A11_update(double A11_old,
                      int p,
                      int q,
                      int m,
                      arma::mat x,
                      arma::mat z,
                      double sigma2_A,
                      arma::vec w,
                      arma::vec gamma,
                      arma::vec beta,
                      arma::mat Lambda,
                      arma::mat delta,
                      arma::vec w1,
                      double metrop_var_A11_trans,
                      int acctot_A11_trans);

Rcpp::List A22_update(double A22_old,
                      int q,
                      int m,
                      double sigma2_A,
                      arma::mat delta_star,
                      arma::vec w1,
                      arma::vec w2,
                      double A21_old,
                      double metrop_var_A22_trans,
                      int acctot_A22_trans);
  
double A21_update(int q,
                  double sigma2_A,
                  arma::mat delta_star,
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

double neg_two_loglike_update(int n,
                              int p,
                              int q,
                              int m,
                              arma::vec y,
                              arma::mat x,
                              arma::mat z, 
                              arma::vec beta,
                              arma::mat Lambda,
                              arma::vec eta_full);

Rcpp::List CWVSmix(int mcmc_samples,
                   int stable,
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
                   Rcpp::Nullable<double> sigma2_beta_prior,
                   Rcpp::Nullable<Rcpp::NumericMatrix> alpha_Lambda_prior,
                   Rcpp::Nullable<double> sigma2_A_prior,
                   Rcpp::Nullable<double> alpha_phi1_prior,
                   Rcpp::Nullable<double> beta_phi1_prior,
                   Rcpp::Nullable<double> alpha_phi2_prior,
                   Rcpp::Nullable<double> beta_phi2_prior,
                   Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                   Rcpp::Nullable<Rcpp::NumericMatrix> Lambda_init,
                   Rcpp::Nullable<Rcpp::NumericMatrix> delta_init,
                   Rcpp::Nullable<Rcpp::NumericVector> w1_init,
                   Rcpp::Nullable<Rcpp::NumericVector> w2_init,
                   Rcpp::Nullable<double> A11_init,
                   Rcpp::Nullable<double> A22_init,
                   Rcpp::Nullable<double> A21_init,
                   Rcpp::Nullable<double> phi1_init,
                   Rcpp::Nullable<double> phi2_init); 

#endif // __CWVSmix__
