#include "../inst/include/shrinkTVP.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "shrinkTVP.h"
#include "TG_sampling_functions.h"
#include "TG_MH_step.h"
#include "TG_log_ratio_value.h"
#include "sample_TG_TVP.h"
#include "sample_DG_TVP.h"
#include "sample_beta_McCausland.h"
#include "DG_sampling_functions.h"
#include "DG_MH_step.h"
#include "DG_log_ratio_value.h"
#include "cpp_utilities.h"
#include "common_sampling_functions.h"
#include "do_rgig1.h"

using namespace Rcpp;

// pred_dens_mix_approx
arma::vec pred_dens_mix_approx(arma::vec x_test, arma::vec y_test, arma::mat theta_sr, arma::mat beta_mean, arma::vec sig2_samp, bool sv, arma::vec sv_phi, arma::vec sv_mu, arma::vec sv_sigma2, arma::cube chol_C_N_inv_samp, arma::cube m_N_samp, int M, bool log);
RcppExport SEXP _shrinkTVP_pred_dens_mix_approx(SEXP x_testSEXP, SEXP y_testSEXP, SEXP theta_srSEXP, SEXP beta_meanSEXP, SEXP sig2_sampSEXP, SEXP svSEXP, SEXP sv_phiSEXP, SEXP sv_muSEXP, SEXP sv_sigma2SEXP, SEXP chol_C_N_inv_sampSEXP, SEXP m_N_sampSEXP, SEXP MSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x_test(x_testSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y_test(y_testSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta_sr(theta_srSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_mean(beta_meanSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sig2_samp(sig2_sampSEXP);
    Rcpp::traits::input_parameter< bool >::type sv(svSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sv_phi(sv_phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sv_mu(sv_muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sv_sigma2(sv_sigma2SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type chol_C_N_inv_samp(chol_C_N_inv_sampSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type m_N_samp(m_N_sampSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_dens_mix_approx(x_test, y_test, theta_sr, beta_mean, sig2_samp, sv, sv_phi, sv_mu, sv_sigma2, chol_C_N_inv_samp, m_N_samp, M, log));
    return rcpp_result_gen;
END_RCPP
}
// calc_fitted_cpp
arma::mat calc_fitted_cpp(arma::vec y, arma::mat x, Rcpp::List beta);
RcppExport SEXP _shrinkTVP_calc_fitted_cpp(SEXP ySEXP, SEXP xSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_fitted_cpp(y, x, beta));
    return rcpp_result_gen;
END_RCPP
}
// shrinkTVP_cpp
List shrinkTVP_cpp(arma::vec y, arma::mat x, std::string mod_type, int niter, int nburn, int nthin, double c0, double g0, double G0, double d1, double d2, double e1, double e2, bool learn_lambda2_B, bool learn_kappa2_B, double lambda2_B, double kappa2_B, bool learn_a_xi, bool learn_a_tau, double a_xi, double a_tau, bool learn_c_xi, bool learn_c_tau, double c_xi, double c_tau, bool a_eq_c_xi, bool a_eq_c_tau, double a_tuning_par_xi, double a_tuning_par_tau, double c_tuning_par_xi, double c_tuning_par_tau, double beta_a_xi, double beta_a_tau, double alpha_a_xi, double alpha_a_tau, double beta_c_xi, double beta_c_tau, double alpha_c_xi, double alpha_c_tau, bool display_progress, bool sv, double Bsigma_sv, double a0_sv, double b0_sv, double bmu, double Bmu, arma::vec adaptive, arma::vec target_rates, arma::vec max_adapts, arma::ivec batch_sizes, Rcpp::List starting_vals);
RcppExport SEXP _shrinkTVP_shrinkTVP_cpp(SEXP ySEXP, SEXP xSEXP, SEXP mod_typeSEXP, SEXP niterSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP c0SEXP, SEXP g0SEXP, SEXP G0SEXP, SEXP d1SEXP, SEXP d2SEXP, SEXP e1SEXP, SEXP e2SEXP, SEXP learn_lambda2_BSEXP, SEXP learn_kappa2_BSEXP, SEXP lambda2_BSEXP, SEXP kappa2_BSEXP, SEXP learn_a_xiSEXP, SEXP learn_a_tauSEXP, SEXP a_xiSEXP, SEXP a_tauSEXP, SEXP learn_c_xiSEXP, SEXP learn_c_tauSEXP, SEXP c_xiSEXP, SEXP c_tauSEXP, SEXP a_eq_c_xiSEXP, SEXP a_eq_c_tauSEXP, SEXP a_tuning_par_xiSEXP, SEXP a_tuning_par_tauSEXP, SEXP c_tuning_par_xiSEXP, SEXP c_tuning_par_tauSEXP, SEXP beta_a_xiSEXP, SEXP beta_a_tauSEXP, SEXP alpha_a_xiSEXP, SEXP alpha_a_tauSEXP, SEXP beta_c_xiSEXP, SEXP beta_c_tauSEXP, SEXP alpha_c_xiSEXP, SEXP alpha_c_tauSEXP, SEXP display_progressSEXP, SEXP svSEXP, SEXP Bsigma_svSEXP, SEXP a0_svSEXP, SEXP b0_svSEXP, SEXP bmuSEXP, SEXP BmuSEXP, SEXP adaptiveSEXP, SEXP target_ratesSEXP, SEXP max_adaptsSEXP, SEXP batch_sizesSEXP, SEXP starting_valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type mod_type(mod_typeSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< double >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< double >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< double >::type G0(G0SEXP);
    Rcpp::traits::input_parameter< double >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< double >::type d2(d2SEXP);
    Rcpp::traits::input_parameter< double >::type e1(e1SEXP);
    Rcpp::traits::input_parameter< double >::type e2(e2SEXP);
    Rcpp::traits::input_parameter< bool >::type learn_lambda2_B(learn_lambda2_BSEXP);
    Rcpp::traits::input_parameter< bool >::type learn_kappa2_B(learn_kappa2_BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2_B(lambda2_BSEXP);
    Rcpp::traits::input_parameter< double >::type kappa2_B(kappa2_BSEXP);
    Rcpp::traits::input_parameter< bool >::type learn_a_xi(learn_a_xiSEXP);
    Rcpp::traits::input_parameter< bool >::type learn_a_tau(learn_a_tauSEXP);
    Rcpp::traits::input_parameter< double >::type a_xi(a_xiSEXP);
    Rcpp::traits::input_parameter< double >::type a_tau(a_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type learn_c_xi(learn_c_xiSEXP);
    Rcpp::traits::input_parameter< bool >::type learn_c_tau(learn_c_tauSEXP);
    Rcpp::traits::input_parameter< double >::type c_xi(c_xiSEXP);
    Rcpp::traits::input_parameter< double >::type c_tau(c_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type a_eq_c_xi(a_eq_c_xiSEXP);
    Rcpp::traits::input_parameter< bool >::type a_eq_c_tau(a_eq_c_tauSEXP);
    Rcpp::traits::input_parameter< double >::type a_tuning_par_xi(a_tuning_par_xiSEXP);
    Rcpp::traits::input_parameter< double >::type a_tuning_par_tau(a_tuning_par_tauSEXP);
    Rcpp::traits::input_parameter< double >::type c_tuning_par_xi(c_tuning_par_xiSEXP);
    Rcpp::traits::input_parameter< double >::type c_tuning_par_tau(c_tuning_par_tauSEXP);
    Rcpp::traits::input_parameter< double >::type beta_a_xi(beta_a_xiSEXP);
    Rcpp::traits::input_parameter< double >::type beta_a_tau(beta_a_tauSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_a_xi(alpha_a_xiSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_a_tau(alpha_a_tauSEXP);
    Rcpp::traits::input_parameter< double >::type beta_c_xi(beta_c_xiSEXP);
    Rcpp::traits::input_parameter< double >::type beta_c_tau(beta_c_tauSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_c_xi(alpha_c_xiSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_c_tau(alpha_c_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< bool >::type sv(svSEXP);
    Rcpp::traits::input_parameter< double >::type Bsigma_sv(Bsigma_svSEXP);
    Rcpp::traits::input_parameter< double >::type a0_sv(a0_svSEXP);
    Rcpp::traits::input_parameter< double >::type b0_sv(b0_svSEXP);
    Rcpp::traits::input_parameter< double >::type bmu(bmuSEXP);
    Rcpp::traits::input_parameter< double >::type Bmu(BmuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type target_rates(target_ratesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type max_adapts(max_adaptsSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type batch_sizes(batch_sizesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type starting_vals(starting_valsSEXP);
    rcpp_result_gen = Rcpp::wrap(shrinkTVP_cpp(y, x, mod_type, niter, nburn, nthin, c0, g0, G0, d1, d2, e1, e2, learn_lambda2_B, learn_kappa2_B, lambda2_B, kappa2_B, learn_a_xi, learn_a_tau, a_xi, a_tau, learn_c_xi, learn_c_tau, c_xi, c_tau, a_eq_c_xi, a_eq_c_tau, a_tuning_par_xi, a_tuning_par_tau, c_tuning_par_xi, c_tuning_par_tau, beta_a_xi, beta_a_tau, alpha_a_xi, alpha_a_tau, beta_c_xi, beta_c_tau, alpha_c_xi, alpha_c_tau, display_progress, sv, Bsigma_sv, a0_sv, b0_sv, bmu, Bmu, adaptive, target_rates, max_adapts, batch_sizes, starting_vals));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _shrinkTVP_RcppExport_validate(const char* sig) {
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Rcpp::List(*shrinkTVP_cpp)(arma::vec,arma::mat,std::string,int,int,int,double,double,double,double,double,double,double,bool,bool,double,double,bool,bool,double,double,bool,bool,double,double,bool,bool,double,double,double,double,double,double,double,double,double,double,double,double,bool,bool,double,double,double,double,double,arma::vec,arma::vec,arma::vec,arma::ivec,Rcpp::List)");
        signatures.insert("double(*DG_log_ratio_value_marginalBFS)(double,double,double,const arma::vec&,double,double)");
        signatures.insert("void(*DG_sample_local_shrink)(arma::vec&,const arma::vec&,double,double)");
        signatures.insert("double(*DG_sample_global_shrink)(const arma::vec&,double,double,double)");
        signatures.insert("double(*TG_MH_step)(double,double,double,const arma::vec&,const arma::vec&,double,double,bool,double,double,bool,bool,arma::vec&,double&,double,double,int&,int,int&)");
        signatures.insert("double(*TG_log_ratio_value_marginalBFS)(double,double,double,const arma::vec&,const arma::vec&,double,double,double,double,bool)");
        signatures.insert("double(*TG_log_ratio_value_tg)(double,double,double,const arma::vec&,const arma::vec&,double,double,double,double)");
        signatures.insert("void(*TG_sample_prior_var_til)(arma::vec&,const arma::vec&,const arma::vec&,double,double,double)");
        signatures.insert("double(*TG_sample_global_shrink)(const arma::vec&,const arma::vec&,const arma::vec&,double,double,double,bool)");
        signatures.insert("void(*TG_sample_local_shrink)(arma::vec&,const arma::vec&,const arma::vec&,double,double,double)");
        signatures.insert("double(*TG_sample_d2)(double,double,double)");
        signatures.insert("void(*sample_alpha)(arma::vec&,arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&,const arma::vec&)");
        signatures.insert("void(*resample_alpha)(arma::vec&,arma::vec&,const arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&)");
        signatures.insert("void(*sample_sigma2)(arma::vec&,const arma::vec&,const arma::mat&,const arma::mat,const arma::vec,const arma::vec&,double,double)");
        signatures.insert("double(*sample_C0)(const arma::vec&,double,double,double)");
        signatures.insert("void(*res_protector)(double&)");
        signatures.insert("void(*calc_xi2_tau2)(arma::vec&,const arma::vec&,const arma::vec&,double,double,double)");
        signatures.insert("void(*to_CP)(arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&)");
        signatures.insert("void(*to_NCP)(arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&)");
        signatures.insert("arma::mat(*robust_chol)(const arma::mat&)");
        signatures.insert("arma::mat(*robust_chol_nontri)(const arma::mat&)");
        signatures.insert("double(*unur_bessel_k_nuasympt)(double,double,bool,bool)");
        signatures.insert("void(*sample_lin_reg_stab)(arma::vec&,const arma::vec&,const arma::mat&,const arma::vec&,const arma::vec&)");
        signatures.insert("void(*sample_lin_reg_rue)(arma::vec&,const arma::vec&,const arma::mat&,const arma::vec&,const arma::vec&)");
        signatures.insert("void(*sample_lin_reg_rue_homosc)(arma::vec&,const arma::vec&,const arma::mat&,double,const arma::vec&)");
        signatures.insert("void(*sample_lin_reg_bhat)(arma::vec&,const arma::vec&,const arma::mat&,const arma::vec&,const arma::vec&)");
        signatures.insert("void(*sample_DG_TVP)(const arma::vec&,const arma::vec&,arma::vec&,arma::vec&,double&,double&,double&,double,double,double&,double,double,double,double,double,double,bool,bool,bool,bool,double,double,const arma::vec&,arma::mat&,arma::vec&,const arma::vec&,const arma::vec&,arma::ivec&,const arma::ivec&,arma::ivec&,int,bool&,std::string&,int&)");
        signatures.insert("void(*sample_TG_TVP)(const arma::vec&,const arma::vec&,arma::vec&,arma::vec&,arma::vec&,arma::vec&,arma::vec&,arma::vec&,arma::vec&,arma::vec&,double&,double&,double&,double,double,double&,double,double,double&,double&,double&,double&,double,double,double,double,bool,bool,bool,bool,bool,bool,double,double,double,double,bool,bool,const arma::vec&,arma::mat&,arma::vec&,const arma::vec&,const arma::vec&,arma::ivec&,const arma::ivec&,arma::ivec&,int,bool&,std::string&,int&)");
        signatures.insert("void(*sample_beta_McCausland)(arma::mat&,const arma::vec&,const arma::mat&,const arma::vec&,const arma::vec&,const arma::vec&,arma::vec&,arma::mat&)");
        signatures.insert("void(*FFBS)(arma::mat&,const arma::vec&,const arma::mat&,const arma::vec&,const arma::vec&,const arma::vec&,arma::vec&,arma::mat&)");
        signatures.insert("double(*do_rgig1)(double,double,double)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _shrinkTVP_RcppExport_registerCCallable() {
    R_RegisterCCallable("shrinkTVP", "shrinkTVP_cpp", (DL_FUNC)shrinkTVP_cpp);
    R_RegisterCCallable("shrinkTVP", "DG_log_ratio_value_marginalBFS", (DL_FUNC)DG_log_ratio_value_marginalBFS);
    R_RegisterCCallable("shrinkTVP", "DG_sample_local_shrink", (DL_FUNC)DG_sample_local_shrink);
    R_RegisterCCallable("shrinkTVP", "DG_sample_global_shrink", (DL_FUNC)DG_sample_global_shrink);
    R_RegisterCCallable("shrinkTVP", "TG_MH_step", (DL_FUNC)TG_MH_step);
    R_RegisterCCallable("shrinkTVP", "TG_log_ratio_value_marginalBFS", (DL_FUNC)TG_log_ratio_value_marginalBFS);
    R_RegisterCCallable("shrinkTVP", "TG_log_ratio_value_tg", (DL_FUNC)TG_log_ratio_value_tg);
    R_RegisterCCallable("shrinkTVP", "TG_sample_prior_var_til", (DL_FUNC)TG_sample_prior_var_til);
    R_RegisterCCallable("shrinkTVP", "TG_sample_global_shrink", (DL_FUNC)TG_sample_global_shrink);
    R_RegisterCCallable("shrinkTVP", "TG_sample_local_shrink", (DL_FUNC)TG_sample_local_shrink);
    R_RegisterCCallable("shrinkTVP", "TG_sample_d2", (DL_FUNC)TG_sample_d2);
    R_RegisterCCallable("shrinkTVP", "sample_alpha", (DL_FUNC)sample_alpha);
    R_RegisterCCallable("shrinkTVP", "resample_alpha", (DL_FUNC)resample_alpha);
    R_RegisterCCallable("shrinkTVP", "sample_sigma2", (DL_FUNC)sample_sigma2);
    R_RegisterCCallable("shrinkTVP", "sample_C0", (DL_FUNC)sample_C0);
    R_RegisterCCallable("shrinkTVP", "res_protector", (DL_FUNC)res_protector);
    R_RegisterCCallable("shrinkTVP", "calc_xi2_tau2", (DL_FUNC)calc_xi2_tau2);
    R_RegisterCCallable("shrinkTVP", "to_CP", (DL_FUNC)to_CP);
    R_RegisterCCallable("shrinkTVP", "to_NCP", (DL_FUNC)to_NCP);
    R_RegisterCCallable("shrinkTVP", "robust_chol", (DL_FUNC)robust_chol);
    R_RegisterCCallable("shrinkTVP", "robust_chol_nontri", (DL_FUNC)robust_chol_nontri);
    R_RegisterCCallable("shrinkTVP", "unur_bessel_k_nuasympt", (DL_FUNC)unur_bessel_k_nuasympt);
    R_RegisterCCallable("shrinkTVP", "sample_lin_reg_stab", (DL_FUNC)sample_lin_reg_stab);
    R_RegisterCCallable("shrinkTVP", "sample_lin_reg_rue", (DL_FUNC)sample_lin_reg_rue);
    R_RegisterCCallable("shrinkTVP", "sample_lin_reg_rue_homosc", (DL_FUNC)sample_lin_reg_rue_homosc);
    R_RegisterCCallable("shrinkTVP", "sample_lin_reg_bhat", (DL_FUNC)sample_lin_reg_bhat);
    R_RegisterCCallable("shrinkTVP", "sample_DG_TVP", (DL_FUNC)sample_DG_TVP);
    R_RegisterCCallable("shrinkTVP", "sample_TG_TVP", (DL_FUNC)sample_TG_TVP);
    R_RegisterCCallable("shrinkTVP", "sample_beta_McCausland", (DL_FUNC)sample_beta_McCausland);
    R_RegisterCCallable("shrinkTVP", "FFBS", (DL_FUNC)FFBS);
    R_RegisterCCallable("shrinkTVP", "do_rgig1", (DL_FUNC)do_rgig1);
    R_RegisterCCallable("shrinkTVP", "_shrinkTVP_RcppExport_validate", (DL_FUNC)_shrinkTVP_RcppExport_validate);
    return R_NilValue;
}


static const R_CallMethodDef CallEntries[] = {
    {"_shrinkTVP_RcppExport_registerCCallable", (DL_FUNC) &_shrinkTVP_RcppExport_registerCCallable, 0},
    {"_shrinkTVP_pred_dens_mix_approx", (DL_FUNC) &_shrinkTVP_pred_dens_mix_approx, 13},
    {"_shrinkTVP_calc_fitted_cpp", (DL_FUNC) &_shrinkTVP_calc_fitted_cpp, 3},
    {"_shrinkTVP_shrinkTVP_cpp", (DL_FUNC) &_shrinkTVP_shrinkTVP_cpp, 51},
    {NULL, NULL, 0}
};

RcppExport void R_init_shrinkTVP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
