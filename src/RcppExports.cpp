// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// A11_update
Rcpp::List A11_update(double A11_old, int n, int p, int m, arma::mat x, arma::mat z, arma::vec off_set, double sigma2_A, arma::vec w, arma::vec gamma, arma::vec beta, arma::vec delta, arma::vec w1, arma::mat risk_sum, double metrop_var_A11_trans, int acctot_A11_trans);
RcppExport SEXP _CWVSmix_A11_update(SEXP A11_oldSEXP, SEXP nSEXP, SEXP pSEXP, SEXP mSEXP, SEXP xSEXP, SEXP zSEXP, SEXP off_setSEXP, SEXP sigma2_ASEXP, SEXP wSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP deltaSEXP, SEXP w1SEXP, SEXP risk_sumSEXP, SEXP metrop_var_A11_transSEXP, SEXP acctot_A11_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type A11_old(A11_oldSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type off_set(off_setSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_A(sigma2_ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type risk_sum(risk_sumSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_A11_trans(metrop_var_A11_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_A11_trans(acctot_A11_transSEXP);
    rcpp_result_gen = Rcpp::wrap(A11_update(A11_old, n, p, m, x, z, off_set, sigma2_A, w, gamma, beta, delta, w1, risk_sum, metrop_var_A11_trans, acctot_A11_trans));
    return rcpp_result_gen;
END_RCPP
}
// A21_update
double A21_update(double sigma2_A, arma::vec delta_star, arma::vec w1, arma::vec w2, double A22);
RcppExport SEXP _CWVSmix_A21_update(SEXP sigma2_ASEXP, SEXP delta_starSEXP, SEXP w1SEXP, SEXP w2SEXP, SEXP A22SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type sigma2_A(sigma2_ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta_star(delta_starSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w2(w2SEXP);
    Rcpp::traits::input_parameter< double >::type A22(A22SEXP);
    rcpp_result_gen = Rcpp::wrap(A21_update(sigma2_A, delta_star, w1, w2, A22));
    return rcpp_result_gen;
END_RCPP
}
// A22_update
Rcpp::List A22_update(double A22_old, int m, double sigma2_A, arma::vec delta_star, arma::vec w1, arma::vec w2, double A21_old, double metrop_var_A22_trans, int acctot_A22_trans);
RcppExport SEXP _CWVSmix_A22_update(SEXP A22_oldSEXP, SEXP mSEXP, SEXP sigma2_ASEXP, SEXP delta_starSEXP, SEXP w1SEXP, SEXP w2SEXP, SEXP A21_oldSEXP, SEXP metrop_var_A22_transSEXP, SEXP acctot_A22_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_A(sigma2_ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta_star(delta_starSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w2(w2SEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_A22_trans(metrop_var_A22_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_A22_trans(acctot_A22_transSEXP);
    rcpp_result_gen = Rcpp::wrap(A22_update(A22_old, m, sigma2_A, delta_star, w1, w2, A21_old, metrop_var_A22_trans, acctot_A22_trans));
    return rcpp_result_gen;
END_RCPP
}
// CWVSmix
Rcpp::List CWVSmix(int mcmc_samples, int p, arma::vec y, arma::mat x, arma::mat z, arma::vec metrop_var_lambda_trans, double metrop_var_rho_trans, double metrop_var_A11_trans, double metrop_var_A22_trans, double metrop_var_phi1_trans, double metrop_var_phi2_trans, int interaction_indicator, int likelihood_indicator, Rcpp::Nullable<Rcpp::NumericVector> offset, Rcpp::Nullable<Rcpp::NumericVector> trials, Rcpp::Nullable<double> a_r_prior, Rcpp::Nullable<double> b_r_prior, Rcpp::Nullable<double> a_sigma2_epsilon_prior, Rcpp::Nullable<double> b_sigma2_epsilon_prior, Rcpp::Nullable<double> sigma2_beta_prior, Rcpp::Nullable<double> alpha_rho_prior, Rcpp::Nullable<double> beta_rho_prior, Rcpp::Nullable<double> sigma2_A_prior, Rcpp::Nullable<double> alpha_phi1_prior, Rcpp::Nullable<double> beta_phi1_prior, Rcpp::Nullable<double> alpha_phi2_prior, Rcpp::Nullable<double> beta_phi2_prior, Rcpp::Nullable<double> r_init, Rcpp::Nullable<double> sigma2_epsilon_init, Rcpp::Nullable<Rcpp::NumericVector> beta_init, Rcpp::Nullable<double> rho_init, Rcpp::Nullable<Rcpp::NumericVector> delta_init, Rcpp::Nullable<Rcpp::NumericVector> w1_init, Rcpp::Nullable<Rcpp::NumericVector> w2_init, Rcpp::Nullable<double> A11_init, Rcpp::Nullable<double> A22_init, Rcpp::Nullable<double> A21_init, Rcpp::Nullable<double> phi1_init, Rcpp::Nullable<double> phi2_init);
RcppExport SEXP _CWVSmix_CWVSmix(SEXP mcmc_samplesSEXP, SEXP pSEXP, SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP metrop_var_lambda_transSEXP, SEXP metrop_var_rho_transSEXP, SEXP metrop_var_A11_transSEXP, SEXP metrop_var_A22_transSEXP, SEXP metrop_var_phi1_transSEXP, SEXP metrop_var_phi2_transSEXP, SEXP interaction_indicatorSEXP, SEXP likelihood_indicatorSEXP, SEXP offsetSEXP, SEXP trialsSEXP, SEXP a_r_priorSEXP, SEXP b_r_priorSEXP, SEXP a_sigma2_epsilon_priorSEXP, SEXP b_sigma2_epsilon_priorSEXP, SEXP sigma2_beta_priorSEXP, SEXP alpha_rho_priorSEXP, SEXP beta_rho_priorSEXP, SEXP sigma2_A_priorSEXP, SEXP alpha_phi1_priorSEXP, SEXP beta_phi1_priorSEXP, SEXP alpha_phi2_priorSEXP, SEXP beta_phi2_priorSEXP, SEXP r_initSEXP, SEXP sigma2_epsilon_initSEXP, SEXP beta_initSEXP, SEXP rho_initSEXP, SEXP delta_initSEXP, SEXP w1_initSEXP, SEXP w2_initSEXP, SEXP A11_initSEXP, SEXP A22_initSEXP, SEXP A21_initSEXP, SEXP phi1_initSEXP, SEXP phi2_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type mcmc_samples(mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type metrop_var_lambda_trans(metrop_var_lambda_transSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_rho_trans(metrop_var_rho_transSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_A11_trans(metrop_var_A11_transSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_A22_trans(metrop_var_A22_transSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_phi1_trans(metrop_var_phi1_transSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_phi2_trans(metrop_var_phi2_transSEXP);
    Rcpp::traits::input_parameter< int >::type interaction_indicator(interaction_indicatorSEXP);
    Rcpp::traits::input_parameter< int >::type likelihood_indicator(likelihood_indicatorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type trials(trialsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type a_r_prior(a_r_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type b_r_prior(b_r_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type a_sigma2_epsilon_prior(a_sigma2_epsilon_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type b_sigma2_epsilon_prior(b_sigma2_epsilon_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_beta_prior(sigma2_beta_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type alpha_rho_prior(alpha_rho_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type beta_rho_prior(beta_rho_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_A_prior(sigma2_A_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type alpha_phi1_prior(alpha_phi1_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type beta_phi1_prior(beta_phi1_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type alpha_phi2_prior(alpha_phi2_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type beta_phi2_prior(beta_phi2_priorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type r_init(r_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type sigma2_epsilon_init(sigma2_epsilon_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type rho_init(rho_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type delta_init(delta_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type w1_init(w1_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type w2_init(w2_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type A11_init(A11_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type A22_init(A22_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type A21_init(A21_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type phi1_init(phi1_initSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type phi2_init(phi2_initSEXP);
    rcpp_result_gen = Rcpp::wrap(CWVSmix(mcmc_samples, p, y, x, z, metrop_var_lambda_trans, metrop_var_rho_trans, metrop_var_A11_trans, metrop_var_A22_trans, metrop_var_phi1_trans, metrop_var_phi2_trans, interaction_indicator, likelihood_indicator, offset, trials, a_r_prior, b_r_prior, a_sigma2_epsilon_prior, b_sigma2_epsilon_prior, sigma2_beta_prior, alpha_rho_prior, beta_rho_prior, sigma2_A_prior, alpha_phi1_prior, beta_phi1_prior, alpha_phi2_prior, beta_phi2_prior, r_init, sigma2_epsilon_init, beta_init, rho_init, delta_init, w1_init, w2_init, A11_init, A22_init, A21_init, phi1_init, phi2_init));
    return rcpp_result_gen;
END_RCPP
}
// beta_update
arma::vec beta_update(int n, int p, int m, int p_x, arma::mat x, arma::mat z, arma::vec off_set, double sigma2_beta, arma::vec w, arma::vec gamma, arma::vec eta_full, arma::mat risk_sum);
RcppExport SEXP _CWVSmix_beta_update(SEXP nSEXP, SEXP pSEXP, SEXP mSEXP, SEXP p_xSEXP, SEXP xSEXP, SEXP zSEXP, SEXP off_setSEXP, SEXP sigma2_betaSEXP, SEXP wSEXP, SEXP gammaSEXP, SEXP eta_fullSEXP, SEXP risk_sumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type p_x(p_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type off_set(off_setSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_beta(sigma2_betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_full(eta_fullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type risk_sum(risk_sumSEXP);
    rcpp_result_gen = Rcpp::wrap(beta_update(n, p, m, p_x, x, z, off_set, sigma2_beta, w, gamma, eta_full, risk_sum));
    return rcpp_result_gen;
END_RCPP
}
// delta_star_update
arma::vec delta_star_update(int m, arma::vec delta, arma::vec w1_old, arma::vec w2_old, double A22_old, double A21_old);
RcppExport SEXP _CWVSmix_delta_star_update(SEXP mSEXP, SEXP deltaSEXP, SEXP w1_oldSEXP, SEXP w2_oldSEXP, SEXP A22_oldSEXP, SEXP A21_oldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w1_old(w1_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w2_old(w2_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    rcpp_result_gen = Rcpp::wrap(delta_star_update(m, delta, w1_old, w2_old, A22_old, A21_old));
    return rcpp_result_gen;
END_RCPP
}
// delta_update
Rcpp::List delta_update(arma::vec delta_old, int n, int p, int m, arma::vec y, arma::mat x, arma::mat z, arma::vec off_set, arma::vec w, arma::vec gamma, arma::vec beta, arma::vec w1_old, arma::vec w2_old, double A11_old, double A22_old, double A21_old, arma::vec eta_full, arma::mat risk_sum);
RcppExport SEXP _CWVSmix_delta_update(SEXP delta_oldSEXP, SEXP nSEXP, SEXP pSEXP, SEXP mSEXP, SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP off_setSEXP, SEXP wSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP w1_oldSEXP, SEXP w2_oldSEXP, SEXP A11_oldSEXP, SEXP A22_oldSEXP, SEXP A21_oldSEXP, SEXP eta_fullSEXP, SEXP risk_sumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type delta_old(delta_oldSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type off_set(off_setSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w1_old(w1_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w2_old(w2_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A11_old(A11_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_full(eta_fullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type risk_sum(risk_sumSEXP);
    rcpp_result_gen = Rcpp::wrap(delta_update(delta_old, n, p, m, y, x, z, off_set, w, gamma, beta, w1_old, w2_old, A11_old, A22_old, A21_old, eta_full, risk_sum));
    return rcpp_result_gen;
END_RCPP
}
// exp_rs
double exp_rs(double a, double b);
RcppExport SEXP _CWVSmix_exp_rs(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_rs(a, b));
    return rcpp_result_gen;
END_RCPP
}
// half_norm_rs
double half_norm_rs(double a, double b);
RcppExport SEXP _CWVSmix_half_norm_rs(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(half_norm_rs(a, b));
    return rcpp_result_gen;
END_RCPP
}
// lambda_update
Rcpp::List lambda_update(arma::mat lambda_star, arma::mat lambda_old, int ind, int p, int q, int m, arma::mat x, arma::mat z, int interaction_indicator, arma::vec off_set, arma::vec w, arma::vec gamma, arma::vec beta, arma::vec eta_full, arma::mat risk_sum, arma::mat temporal_corr_inv_old, double metrop_var_lambda_trans, int acctot_lambda);
RcppExport SEXP _CWVSmix_lambda_update(SEXP lambda_starSEXP, SEXP lambda_oldSEXP, SEXP indSEXP, SEXP pSEXP, SEXP qSEXP, SEXP mSEXP, SEXP xSEXP, SEXP zSEXP, SEXP interaction_indicatorSEXP, SEXP off_setSEXP, SEXP wSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP eta_fullSEXP, SEXP risk_sumSEXP, SEXP temporal_corr_inv_oldSEXP, SEXP metrop_var_lambda_transSEXP, SEXP acctot_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type lambda_star(lambda_starSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lambda_old(lambda_oldSEXP);
    Rcpp::traits::input_parameter< int >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type interaction_indicator(interaction_indicatorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type off_set(off_setSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_full(eta_fullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type risk_sum(risk_sumSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type temporal_corr_inv_old(temporal_corr_inv_oldSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_lambda_trans(metrop_var_lambda_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_lambda(acctot_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda_update(lambda_star, lambda_old, ind, p, q, m, x, z, interaction_indicator, off_set, w, gamma, beta, eta_full, risk_sum, temporal_corr_inv_old, metrop_var_lambda_trans, acctot_lambda));
    return rcpp_result_gen;
END_RCPP
}
// neg_two_loglike_update
double neg_two_loglike_update(int n, int p, int m, arma::vec y, arma::mat x, arma::mat z, arma::vec off_set, arma::vec tri_als, int likelihood_indicator, int r, double sigma2_epsilon, arma::vec beta, arma::vec eta_full, arma::mat risk_sum);
RcppExport SEXP _CWVSmix_neg_two_loglike_update(SEXP nSEXP, SEXP pSEXP, SEXP mSEXP, SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP off_setSEXP, SEXP tri_alsSEXP, SEXP likelihood_indicatorSEXP, SEXP rSEXP, SEXP sigma2_epsilonSEXP, SEXP betaSEXP, SEXP eta_fullSEXP, SEXP risk_sumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type off_set(off_setSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tri_als(tri_alsSEXP);
    Rcpp::traits::input_parameter< int >::type likelihood_indicator(likelihood_indicatorSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2_epsilon(sigma2_epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_full(eta_fullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type risk_sum(risk_sumSEXP);
    rcpp_result_gen = Rcpp::wrap(neg_two_loglike_update(n, p, m, y, x, z, off_set, tri_als, likelihood_indicator, r, sigma2_epsilon, beta, eta_full, risk_sum));
    return rcpp_result_gen;
END_RCPP
}
// norm_rs
double norm_rs(double a, double b);
RcppExport SEXP _CWVSmix_norm_rs(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(norm_rs(a, b));
    return rcpp_result_gen;
END_RCPP
}
// phi_update
Rcpp::List phi_update(double phi_old, int m, double alpha_phi, double beta_phi, arma::vec w, Rcpp::List temporal_corr_info, double metrop_var_phi_trans, int acctot_phi_trans);
RcppExport SEXP _CWVSmix_phi_update(SEXP phi_oldSEXP, SEXP mSEXP, SEXP alpha_phiSEXP, SEXP beta_phiSEXP, SEXP wSEXP, SEXP temporal_corr_infoSEXP, SEXP metrop_var_phi_transSEXP, SEXP acctot_phi_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type phi_old(phi_oldSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_phi(alpha_phiSEXP);
    Rcpp::traits::input_parameter< double >::type beta_phi(beta_phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type temporal_corr_info(temporal_corr_infoSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_phi_trans(metrop_var_phi_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_phi_trans(acctot_phi_transSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_update(phi_old, m, alpha_phi, beta_phi, w, temporal_corr_info, metrop_var_phi_trans, acctot_phi_trans));
    return rcpp_result_gen;
END_RCPP
}
// r_update
int r_update(int n, int p, int m, arma::vec y, arma::mat x, arma::mat z, arma::vec off_set, int a_r, int b_r, arma::vec beta, arma::vec eta_full, arma::mat risk_sum);
RcppExport SEXP _CWVSmix_r_update(SEXP nSEXP, SEXP pSEXP, SEXP mSEXP, SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP off_setSEXP, SEXP a_rSEXP, SEXP b_rSEXP, SEXP betaSEXP, SEXP eta_fullSEXP, SEXP risk_sumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type off_set(off_setSEXP);
    Rcpp::traits::input_parameter< int >::type a_r(a_rSEXP);
    Rcpp::traits::input_parameter< int >::type b_r(b_rSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_full(eta_fullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type risk_sum(risk_sumSEXP);
    rcpp_result_gen = Rcpp::wrap(r_update(n, p, m, y, x, z, off_set, a_r, b_r, beta, eta_full, risk_sum));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_pgdraw
arma::vec rcpp_pgdraw(arma::vec b, arma::vec c);
RcppExport SEXP _CWVSmix_rcpp_pgdraw(SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_pgdraw(b, c));
    return rcpp_result_gen;
END_RCPP
}
// rho_update
Rcpp::List rho_update(double rho_old, int p, int m, double alpha_rho, double beta_rho, arma::mat lambda_star, Rcpp::List temporal_corr_info, double metrop_var_rho_trans, int acctot_rho_trans);
RcppExport SEXP _CWVSmix_rho_update(SEXP rho_oldSEXP, SEXP pSEXP, SEXP mSEXP, SEXP alpha_rhoSEXP, SEXP beta_rhoSEXP, SEXP lambda_starSEXP, SEXP temporal_corr_infoSEXP, SEXP metrop_var_rho_transSEXP, SEXP acctot_rho_transSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho_old(rho_oldSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_rho(alpha_rhoSEXP);
    Rcpp::traits::input_parameter< double >::type beta_rho(beta_rhoSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lambda_star(lambda_starSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type temporal_corr_info(temporal_corr_infoSEXP);
    Rcpp::traits::input_parameter< double >::type metrop_var_rho_trans(metrop_var_rho_transSEXP);
    Rcpp::traits::input_parameter< int >::type acctot_rho_trans(acctot_rho_transSEXP);
    rcpp_result_gen = Rcpp::wrap(rho_update(rho_old, p, m, alpha_rho, beta_rho, lambda_star, temporal_corr_info, metrop_var_rho_trans, acctot_rho_trans));
    return rcpp_result_gen;
END_RCPP
}
// rnorm_trunc
double rnorm_trunc(double mu, double sigma, double lower, double upper);
RcppExport SEXP _CWVSmix_rnorm_trunc(SEXP muSEXP, SEXP sigmaSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(rnorm_trunc(mu, sigma, lower, upper));
    return rcpp_result_gen;
END_RCPP
}
// sigma2_epsilon_update
double sigma2_epsilon_update(int n, int p, int m, arma::vec y, arma::mat x, arma::mat z, int likelihood_indicator, double a_sigma2_epsilon, double b_sigma2_epsilon, arma::vec beta_old, arma::vec eta_full, arma::mat risk_sum);
RcppExport SEXP _CWVSmix_sigma2_epsilon_update(SEXP nSEXP, SEXP pSEXP, SEXP mSEXP, SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP likelihood_indicatorSEXP, SEXP a_sigma2_epsilonSEXP, SEXP b_sigma2_epsilonSEXP, SEXP beta_oldSEXP, SEXP eta_fullSEXP, SEXP risk_sumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type likelihood_indicator(likelihood_indicatorSEXP);
    Rcpp::traits::input_parameter< double >::type a_sigma2_epsilon(a_sigma2_epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type b_sigma2_epsilon(b_sigma2_epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_full(eta_fullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type risk_sum(risk_sumSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma2_epsilon_update(n, p, m, y, x, z, likelihood_indicator, a_sigma2_epsilon, b_sigma2_epsilon, beta_old, eta_full, risk_sum));
    return rcpp_result_gen;
END_RCPP
}
// temporal_corr_fun
Rcpp::List temporal_corr_fun(int p_z, double phi);
RcppExport SEXP _CWVSmix_temporal_corr_fun(SEXP p_zSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p_z(p_zSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(temporal_corr_fun(p_z, phi));
    return rcpp_result_gen;
END_RCPP
}
// unif_rs
double unif_rs(double a, double b);
RcppExport SEXP _CWVSmix_unif_rs(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(unif_rs(a, b));
    return rcpp_result_gen;
END_RCPP
}
// w1_update
Rcpp::List w1_update(int n, int p, int m, arma::mat x, arma::mat z, arma::vec off_set, arma::vec w, arma::vec gamma, arma::vec beta, arma::vec delta, arma::vec delta_star, arma::vec w2_old, double A11_old, double A22_old, double A21_old, arma::mat risk_sum, arma::mat corr_inv1);
RcppExport SEXP _CWVSmix_w1_update(SEXP nSEXP, SEXP pSEXP, SEXP mSEXP, SEXP xSEXP, SEXP zSEXP, SEXP off_setSEXP, SEXP wSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP deltaSEXP, SEXP delta_starSEXP, SEXP w2_oldSEXP, SEXP A11_oldSEXP, SEXP A22_oldSEXP, SEXP A21_oldSEXP, SEXP risk_sumSEXP, SEXP corr_inv1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type off_set(off_setSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta_star(delta_starSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w2_old(w2_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A11_old(A11_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type risk_sum(risk_sumSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv1(corr_inv1SEXP);
    rcpp_result_gen = Rcpp::wrap(w1_update(n, p, m, x, z, off_set, w, gamma, beta, delta, delta_star, w2_old, A11_old, A22_old, A21_old, risk_sum, corr_inv1));
    return rcpp_result_gen;
END_RCPP
}
// w2_update
arma::vec w2_update(int p, int m, arma::mat z, arma::vec delta_star, arma::vec w1, double A22_old, double A21_old, arma::mat corr_inv2);
RcppExport SEXP _CWVSmix_w2_update(SEXP pSEXP, SEXP mSEXP, SEXP zSEXP, SEXP delta_starSEXP, SEXP w1SEXP, SEXP A22_oldSEXP, SEXP A21_oldSEXP, SEXP corr_inv2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta_star(delta_starSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< double >::type A22_old(A22_oldSEXP);
    Rcpp::traits::input_parameter< double >::type A21_old(A21_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type corr_inv2(corr_inv2SEXP);
    rcpp_result_gen = Rcpp::wrap(w2_update(p, m, z, delta_star, w1, A22_old, A21_old, corr_inv2));
    return rcpp_result_gen;
END_RCPP
}
// w_update
Rcpp::List w_update(int n, int p, int m, arma::vec y, arma::mat x, arma::mat z, arma::vec off_set, arma::vec tri_als, int likelihood_indicator, int r, arma::vec beta_old, arma::vec eta_full, arma::mat risk_sum);
RcppExport SEXP _CWVSmix_w_update(SEXP nSEXP, SEXP pSEXP, SEXP mSEXP, SEXP ySEXP, SEXP xSEXP, SEXP zSEXP, SEXP off_setSEXP, SEXP tri_alsSEXP, SEXP likelihood_indicatorSEXP, SEXP rSEXP, SEXP beta_oldSEXP, SEXP eta_fullSEXP, SEXP risk_sumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type off_set(off_setSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tri_als(tri_alsSEXP);
    Rcpp::traits::input_parameter< int >::type likelihood_indicator(likelihood_indicatorSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta_full(eta_fullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type risk_sum(risk_sumSEXP);
    rcpp_result_gen = Rcpp::wrap(w_update(n, p, m, y, x, z, off_set, tri_als, likelihood_indicator, r, beta_old, eta_full, risk_sum));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CWVSmix_A11_update", (DL_FUNC) &_CWVSmix_A11_update, 16},
    {"_CWVSmix_A21_update", (DL_FUNC) &_CWVSmix_A21_update, 5},
    {"_CWVSmix_A22_update", (DL_FUNC) &_CWVSmix_A22_update, 9},
    {"_CWVSmix_CWVSmix", (DL_FUNC) &_CWVSmix_CWVSmix, 39},
    {"_CWVSmix_beta_update", (DL_FUNC) &_CWVSmix_beta_update, 12},
    {"_CWVSmix_delta_star_update", (DL_FUNC) &_CWVSmix_delta_star_update, 6},
    {"_CWVSmix_delta_update", (DL_FUNC) &_CWVSmix_delta_update, 18},
    {"_CWVSmix_exp_rs", (DL_FUNC) &_CWVSmix_exp_rs, 2},
    {"_CWVSmix_half_norm_rs", (DL_FUNC) &_CWVSmix_half_norm_rs, 2},
    {"_CWVSmix_lambda_update", (DL_FUNC) &_CWVSmix_lambda_update, 18},
    {"_CWVSmix_neg_two_loglike_update", (DL_FUNC) &_CWVSmix_neg_two_loglike_update, 14},
    {"_CWVSmix_norm_rs", (DL_FUNC) &_CWVSmix_norm_rs, 2},
    {"_CWVSmix_phi_update", (DL_FUNC) &_CWVSmix_phi_update, 8},
    {"_CWVSmix_r_update", (DL_FUNC) &_CWVSmix_r_update, 12},
    {"_CWVSmix_rcpp_pgdraw", (DL_FUNC) &_CWVSmix_rcpp_pgdraw, 2},
    {"_CWVSmix_rho_update", (DL_FUNC) &_CWVSmix_rho_update, 9},
    {"_CWVSmix_rnorm_trunc", (DL_FUNC) &_CWVSmix_rnorm_trunc, 4},
    {"_CWVSmix_sigma2_epsilon_update", (DL_FUNC) &_CWVSmix_sigma2_epsilon_update, 12},
    {"_CWVSmix_temporal_corr_fun", (DL_FUNC) &_CWVSmix_temporal_corr_fun, 2},
    {"_CWVSmix_unif_rs", (DL_FUNC) &_CWVSmix_unif_rs, 2},
    {"_CWVSmix_w1_update", (DL_FUNC) &_CWVSmix_w1_update, 17},
    {"_CWVSmix_w2_update", (DL_FUNC) &_CWVSmix_w2_update, 8},
    {"_CWVSmix_w_update", (DL_FUNC) &_CWVSmix_w_update, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_CWVSmix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
