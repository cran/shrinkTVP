#ifndef shrinkTVP_H
#define shrinkTVP_H
#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace shrinkTVP {

using namespace Rcpp;

inline Rcpp::List shrinkTVP_cpp(arma::vec y,
                                arma::mat x,
                                std::string mod_type,
                                int niter,
                                int nburn,
                                int nthin,
                                double c0,
                                double g0,
                                double G0,
                                double d1,
                                double d2,
                                double e1,
                                double e2,
                                bool learn_lambda2_B,
                                bool learn_kappa2_B,
                                double lambda2_B,
                                double kappa2_B,
                                bool learn_a_xi,
                                bool learn_a_tau,
                                double a_xi,
                                double a_tau,
                                bool learn_c_xi,
                                bool learn_c_tau,
                                double c_xi,
                                double c_tau,
                                bool a_eq_c_xi,
                                bool a_eq_c_tau,
                                double a_tuning_par_xi,
                                double a_tuning_par_tau,
                                double c_tuning_par_xi,
                                double c_tuning_par_tau,
                                double beta_a_xi,
                                double beta_a_tau,
                                double alpha_a_xi,
                                double alpha_a_tau,
                                double beta_c_xi,
                                double beta_c_tau,
                                double alpha_c_xi,
                                double alpha_c_tau,
                                bool display_progress,
                                bool sv,
                                double Bsigma_sv,
                                double a0_sv,
                                double b0_sv,
                                double bmu,
                                double Bmu,
                                arma::vec adaptive,
                                arma::vec target_rates,
                                arma::vec max_adapts,
                                arma::ivec batch_sizes,
                                Rcpp::List starting_vals) {
  typedef Rcpp::List(*Func)(arma::vec,arma::mat,std::string,int,int,int,double,double,double,double,double,double,double,bool,bool,double,double,bool,bool,double,double,bool,bool,double,double,bool,bool,double,double,double,double,double,double,double,double,double,double,double,double,bool,bool,double,double,double,double,double,arma::vec,arma::vec,arma::vec,arma::ivec,Rcpp::List);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "shrinkTVP_cpp");
  }
  {
    return func(y,
                x,
                mod_type,
                niter,
                nburn,
                nthin,
                c0,
                g0,
                G0,
                d1,
                d2,
                e1,
                e2,
                learn_lambda2_B,
                learn_kappa2_B,
                lambda2_B,
                kappa2_B,
                learn_a_xi,
                learn_a_tau,
                a_xi,
                a_tau,
                learn_c_xi,
                learn_c_tau,
                c_xi,
                c_tau,
                a_eq_c_xi,
                a_eq_c_tau,
                a_tuning_par_xi,
                a_tuning_par_tau,
                c_tuning_par_xi,
                c_tuning_par_tau,
                beta_a_xi,
                beta_a_tau,
                alpha_a_xi,
                alpha_a_tau,
                beta_c_xi,
                beta_c_tau,
                alpha_c_xi,
                alpha_c_tau,
                display_progress,
                sv,
                Bsigma_sv,
                a0_sv,
                b0_sv,
                bmu,
                Bmu,
                adaptive,
                target_rates,
                max_adapts,
                batch_sizes,
                starting_vals);
  }
}

inline double DG_log_ratio_value_marginalBFS(double proposal,
                                             double old_val,
                                             double scale_par,
                                             const arma::vec& param_vec,
                                             double b1,
                                             double b2) {
  typedef double(*Func)(double,double,double,const arma::vec&,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "DG_log_ratio_value_marginalBFS");
  }
  {
    return func(proposal,
                old_val,
                scale_par,
                param_vec,
                b1,
                b2);
  }
}

inline void DG_sample_local_shrink(arma::vec& local_shrink,
                                   const arma::vec& param_vec,
                                   double global_shrink,
                                   double a) {
  typedef void(*Func)(arma::vec&,const arma::vec&,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "DG_sample_local_shrink");
  }
  {
    func(local_shrink,
         param_vec,
         global_shrink,
         a);
  }
}

inline double DG_sample_global_shrink(const arma::vec& prior_var,
                                      double a,
                                      double hyper1,
                                      double hyper2) {
  typedef double(*Func)(const arma::vec&,double,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "DG_sample_global_shrink");
  }
  {
    return func(prior_var,
                a,
                hyper1,
                hyper2);
  }
}

inline double TG_MH_step(double current_val,
                         double c_tuning_par,
                         double scale_par,
                         const arma::vec& scale_vec,
                         const arma::vec& param_vec,
                         double b,
                         double nu,
                         bool is_c,
                         double scale_scale,
                         double other_hyp,
                         bool common_shrink_par,
                         bool adaptive,
                         arma::vec& batch,
                         double& curr_sd,
                         double target_rate,
                         double max_adapt,
                         int& batch_nr,
                         int batch_size,
                         int& batch_pos) {
  typedef double(*Func)(double,double,double,const arma::vec&,const arma::vec&,double,double,bool,double,double,bool,bool,arma::vec&,double&,double,double,int&,int,int&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "TG_MH_step");
  }
  {
    return func(current_val,
                c_tuning_par,
                scale_par,
                scale_vec,
                param_vec,
                b,
                nu,
                is_c,
                scale_scale,
                other_hyp,
                common_shrink_par,
                adaptive,
                batch,
                curr_sd,
                target_rate,
                max_adapt,
                batch_nr,
                batch_size,
                batch_pos);
  }
}

inline double TG_log_ratio_value_marginalBFS(double proposal,
                                             double old_val,
                                             double scale_par,
                                             const arma::vec& scale_vec,
                                             const arma::vec& param_vec,
                                             double scale_scale,
                                             double c,
                                             double b1,
                                             double b2,
                                             bool common_shrink_par) {
  typedef double(*Func)(double,double,double,const arma::vec&,const arma::vec&,double,double,double,double,bool);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "TG_log_ratio_value_marginalBFS");
  }
  {
    return func(proposal,
                old_val,
                scale_par,
                scale_vec,
                param_vec,
                scale_scale,
                c,
                b1,
                b2,
                common_shrink_par);
  }
}

inline double TG_log_ratio_value_tg(double proposal,
                                    double old_val,
                                    double scale_par,
                                    const arma::vec& scale_vec,
                                    const arma::vec& param_vec,
                                    double scale_scale,
                                    double a,
                                    double b1,
                                    double b2) {
  typedef double(*Func)(double,double,double,const arma::vec&,const arma::vec&,double,double,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "TG_log_ratio_value_tg");
  }
  {
    return func(proposal,
                old_val,
                scale_par,
                scale_vec,
                param_vec,
                scale_scale,
                a,
                b1,
                b2);
  }
}

inline void TG_sample_prior_var_til(arma::vec& prior_var_til,
                                    const arma::vec& param_vec,
                                    const arma::vec& local_shrink,
                                    double global_shrink,
                                    double a,
                                    double c) {
  typedef void(*Func)(arma::vec&,const arma::vec&,const arma::vec&,double,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "TG_sample_prior_var_til");
  }
  {
    func(prior_var_til,
         param_vec,
         local_shrink,
         global_shrink,
         a,
         c);
  }
}

inline double TG_sample_global_shrink(const arma::vec& prior_var_til,
                                      const arma::vec& local_shrink,
                                      const arma::vec& param_vec,
                                      double a,
                                      double c,
                                      double hyper2,
                                      bool common_shrink_par) {
  typedef double(*Func)(const arma::vec&,const arma::vec&,const arma::vec&,double,double,double,bool);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "TG_sample_global_shrink");
  }
  {
    return func(prior_var_til,
                local_shrink,
                param_vec,
                a,
                c,
                hyper2,
                common_shrink_par);
  }
}

inline void TG_sample_local_shrink(arma::vec& local_shrink,
                                   const arma::vec& param_vec,
                                   const arma::vec& prior_var_til,
                                   double global_shrink,
                                   double c,
                                   double a) {
  typedef void(*Func)(arma::vec&,const arma::vec&,const arma::vec&,double,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "TG_sample_local_shrink");
  }
  {
    func(local_shrink,
         param_vec,
         prior_var_til,
         global_shrink,
         c,
         a);
  }
}

inline double TG_sample_d2(double global_shrink,
                           double a,
                           double c) {
  typedef double(*Func)(double,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "TG_sample_d2");
  }
  {
    return func(global_shrink,
                a,
                c);
  }
}

inline void sample_alpha(arma::vec& beta_mean,
                         arma::vec& theta_sr,
                         const arma::vec& y,
                         const arma::mat& x,
                         const arma::mat& beta_nc,
                         const arma::vec& sigma2,
                         const arma::vec& tau2,
                         const arma::vec& xi2) {
  typedef void(*Func)(arma::vec&,arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&,const arma::vec&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_alpha");
  }
  {
    func(beta_mean,
         theta_sr,
         y,
         x,
         beta_nc,
         sigma2,
         tau2,
         xi2);
  }
}

inline void resample_alpha(arma::vec& beta_mean,
                           arma::vec& theta_sr,
                           const arma::mat& beta,
                           const arma::mat& beta_nc,
                           const arma::vec& xi2,
                           const arma::vec& tau2) {
  typedef void(*Func)(arma::vec&,arma::vec&,const arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "resample_alpha");
  }
  {
    func(beta_mean,
         theta_sr,
         beta,
         beta_nc,
         xi2,
         tau2);
  }
}

inline void sample_sigma2(arma::vec& sigma2,
                          const arma::vec& y,
                          const arma::mat& x,
                          const arma::mat beta_nc,
                          const arma::vec beta_mean,
                          const arma::vec& theta_sr,
                          double c0,
                          double C0) {
  typedef void(*Func)(arma::vec&,const arma::vec&,const arma::mat&,const arma::mat,const arma::vec,const arma::vec&,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_sigma2");
  }
  {
    func(sigma2,
         y,
         x,
         beta_nc,
         beta_mean,
         theta_sr,
         c0,
         C0);
  }
}

inline double sample_C0(const arma::vec& sigma2,
                        double g0,
                        double c0,
                        double G0) {
  typedef double(*Func)(const arma::vec&,double,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_C0");
  }
  {
    return func(sigma2,
                g0,
                c0,
                G0);
  }
}

inline void res_protector(double& x) {
  typedef void(*Func)(double&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "res_protector");
  }
  {
    func(x);
  }
}



inline void calc_xi2_tau2(arma::vec& param,
                          const arma::vec& param_til,
                          const arma::vec& loc_shrink_til,
                          double glob_shrink,
                          double c,
                          double a) {
  typedef void(*Func)(arma::vec&,const arma::vec&,const arma::vec&,double,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "calc_xi2_tau2");
  }
  {
    func(param,
         param_til,
         loc_shrink_til,
         glob_shrink,
         c,
         a);
  }
}

inline void calc_xi2_til_tau2_til(arma::vec& param_til,
                                  const arma::vec& param,
                                  const arma::vec& loc_shrink_til,
                                  double glob_shrink,
                                  double c,
                                  double a) {
  typedef void(*Func)(arma::vec&,const arma::vec&,const arma::vec&,double,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "calc_xi2_til_tau2_til");
  }
  {
    func(param_til,
         param,
         loc_shrink_til,
         glob_shrink,
         c,
         a);
  }
}

inline void calc_kappa2_lambda2(arma::vec& param,
                                const arma::vec param_til,
                                double glob_shrink,
                                double c) {
  typedef void(*Func)(arma::vec&,const arma::vec,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "calc_kappa2_lambda2");
  }
  {
    func(param,
         param_til,
         glob_shrink,
         c);
  }
}

inline void calc_kappa2_til_lambda2_til(arma::vec& param_til,
                                        const arma::vec& param,
                                        double glob_shrink,
                                        double c) {
  typedef void(*Func)(arma::vec&,const arma::vec&,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "calc_kappa2_til_lambda2_til");
  }
  {
    func(param_til,
         param,
         glob_shrink,
         c);
  }
}

inline void to_CP(arma::mat& beta,
                  const arma::mat& beta_nc,
                  const arma::vec& theta_sr,
                  const arma::vec& beta_mean) {
  typedef void(*Func)(arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "to_CP");
  }
  {
    func(beta,
         beta_nc,
         theta_sr,
         beta_mean);
  }
}

inline void to_NCP(arma::mat& beta_nc,
                   const arma::mat& beta,
                   const arma::vec& theta_sr,
                   const arma::vec& beta_mean) {
  typedef void(*Func)(arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "to_NCP");
  }
  {
    func(beta_nc,
         beta,
         theta_sr,
         beta_mean);
  }
}

inline arma::mat robust_chol (const arma::mat& V) {
  typedef arma::mat(*Func)(const arma::mat&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "robust_chol");
  }
  {
    return func(V);
  }
}

inline arma::mat robust_chol_nontri (const arma::mat& V) {
  typedef arma::mat(*Func)(const arma::mat&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "robust_chol_nontri");
  }
  {
    return func(V);
  }
}

inline double unur_bessel_k_nuasympt(double x,
                                     double nu,
                                     bool islog,
                                     bool expon_scaled) {
  typedef double(*Func)(double,double,bool,bool);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "unur_bessel_k_nuasympt");
  }
  {
    return func(x,
                nu,
                islog,
                expon_scaled);
  }
}

inline void sample_lin_reg_stab(arma::vec& param_vec,
                                const arma::vec& y,
                                const arma::mat& X,
                                const arma::vec& sigma2,
                                const arma::vec& prior_var) {
  typedef void(*Func)(arma::vec&,const arma::vec&,const arma::mat&,const arma::vec&,const arma::vec&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_lin_reg_stab");
  }
  {
    func(param_vec,
         y,
         X,
         sigma2,
         prior_var);
  }
}

inline void sample_lin_reg_rue(arma::vec& param_vec,
                               const arma::vec& y,
                               const arma::mat& X,
                               const arma::vec& sigma2,
                               const arma::vec& prior_var) {
  typedef void(*Func)(arma::vec&,const arma::vec&,const arma::mat&,const arma::vec&,const arma::vec&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_lin_reg_rue");
  }
  {
    func(param_vec,
         y,
         X,
         sigma2,
         prior_var);
  }
}

inline void sample_lin_reg_rue_homosc(arma::vec& param_vec,
                                      const arma::vec& Xty,
                                      const arma::mat& XtX,
                                      double sigma2,
                                      const arma::vec& prior_var) {
  typedef void(*Func)(arma::vec&,const arma::vec&,const arma::mat&,double,const arma::vec&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_lin_reg_rue_homosc");
  }
  {
    func(param_vec,
         Xty,
         XtX,
         sigma2,
         prior_var);
  }
}

inline void sample_lin_reg_bhat(arma::vec& param_vec,
                                const arma::vec& y,
                                const arma::mat& x,
                                double sigma2,
                                const arma::vec& prior_var) {
  typedef void(*Func)(arma::vec&,const arma::vec&,const arma::mat&,double,const arma::vec&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_lin_reg_bhat");
  }
  {
    func(param_vec,
         y,
         x,
         sigma2,
         prior_var);
  }
}

inline void sample_DG_TVP(const arma::vec& beta_mean,
                          const arma::vec& theta_sr,
                          arma::vec& tau2,
                          arma::vec& xi2,
                          double& lambda2_B,
                          double& kappa2_B,
                          double& a_xi,
                          double beta_a_xi,
                          double alpha_a_xi,
                          double& a_tau,
                          double beta_a_tau,
                          double alpha_a_tau,
                          double d1,
                          double d2,
                          double e1,
                          double e2,
                          bool learn_kappa2_B,
                          bool learn_lambda2_B,
                          bool learn_a_xi,
                          bool learn_a_tau,
                          double a_tuning_par_xi,
                          double a_tuning_par_tau,
                          const arma::vec& adaptive,
                          arma::mat& batches,
                          arma::vec& curr_sds,
                          const arma::vec& target_rates,
                          const arma::vec& max_adapts,
                          arma::ivec& batch_nrs,
                          const arma::ivec& batch_sizes,
                          arma::ivec& batch_pos,
                          int j,
                          bool& succesful,
                          std::string& fail,
                          int& fail_iter) {
  typedef void(*Func)(const arma::vec&,const arma::vec&,arma::vec&,arma::vec&,double&,double&,double&,double,double,double&,double,double,double,double,double,double,bool,bool,bool,bool,double,double,const arma::vec&,arma::mat&,arma::vec&,const arma::vec&,const arma::vec&,arma::ivec&,const arma::ivec&,arma::ivec&,int,bool&,std::string&,int&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_DG_TVP");
  }
  {
    func(beta_mean,
         theta_sr,
         tau2,
         xi2,
         lambda2_B,
         kappa2_B,
         a_xi,
         beta_a_xi,
         alpha_a_xi,
         a_tau,
         beta_a_tau,
         alpha_a_tau,
         d1,
         d2,
         e1,
         e2,
         learn_kappa2_B,
         learn_lambda2_B,
         learn_a_xi,
         learn_a_tau,
         a_tuning_par_xi,
         a_tuning_par_tau,
         adaptive,
         batches,
         curr_sds,
         target_rates,
         max_adapts,
         batch_nrs,
         batch_sizes,
         batch_pos,
         j,
         succesful,
         fail,
         fail_iter);
  }
}

inline void sample_TG_TVP(const arma::vec& beta_mean,
                          const arma::vec& theta_sr,
                          arma::vec& tau2,
                          arma::vec& xi2,
                          arma::vec& tau2_til,
                          arma::vec& xi2_til,
                          arma::vec& lambda2_til,
                          arma::vec& kappa2_til,
                          double& lambda2_B,
                          double& kappa2_B,
                          double& a_xi,
                          double beta_a_xi,
                          double alpha_a_xi,
                          double& a_tau,
                          double beta_a_tau,
                          double alpha_a_tau,
                          double& d2,
                          double& e2,
                          double& c_xi,
                          double& c_tau,
                          double beta_c_xi,
                          double alpha_c_xi,
                          double beta_c_tau,
                          double alpha_c_tau,
                          bool learn_kappa2_B,
                          bool learn_lambda2_B,
                          bool learn_a_xi,
                          bool learn_a_tau,
                          bool learn_c_xi,
                          bool learn_c_tau,
                          double a_tuning_par_xi,
                          double a_tuning_par_tau,
                          double c_tuning_par_xi,
                          double c_tuning_par_tau,
                          bool a_eq_c_xi,
                          bool a_eq_c_tau,
                          const arma::vec& adaptive,
                          arma::mat& batches,
                          arma::vec& curr_sds,
                          const arma::vec& target_rates,
                          const arma::vec& max_adapts,
                          arma::ivec& batch_nrs,
                          const arma::ivec& batch_sizes,
                          arma::ivec& batch_pos,
                          int j,
                          bool& succesful,
                          std::string& fail,
                          int& fail_iter) {
  typedef void(*Func)(const arma::vec&,const arma::vec&,arma::vec&,arma::vec&,arma::vec&,arma::vec&,arma::vec&,arma::vec&,double&,double&,double&,double,double,double&,double,double,double&,double&,double&,double&,double,double,double,double,bool,bool,bool,bool,bool,bool,double,double,double,double,bool,bool,const arma::vec&,arma::mat&,arma::vec&,const arma::vec&,const arma::vec&,arma::ivec&,const arma::ivec&,arma::ivec&,int,bool&,std::string&,int&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_TG_TVP");
  }
  {
    func(beta_mean,
         theta_sr,
         tau2,
         xi2,
         tau2_til,
         xi2_til,
         lambda2_til,
         kappa2_til,
         lambda2_B,
         kappa2_B,
         a_xi,
         beta_a_xi,
         alpha_a_xi,
         a_tau,
         beta_a_tau,
         alpha_a_tau,
         d2,
         e2,
         c_xi,
         c_tau,
         beta_c_xi,
         alpha_c_xi,
         beta_c_tau,
         alpha_c_tau,
         learn_kappa2_B,
         learn_lambda2_B,
         learn_a_xi,
         learn_a_tau,
         learn_c_xi,
         learn_c_tau,
         a_tuning_par_xi,
         a_tuning_par_tau,
         c_tuning_par_xi,
         c_tuning_par_tau,
         a_eq_c_xi,
         a_eq_c_tau,
         adaptive,
         batches,
         curr_sds,
         target_rates,
         max_adapts,
         batch_nrs,
         batch_sizes,
         batch_pos,
         j,
         succesful,
         fail,
         fail_iter);
  }
}

inline void sample_beta_McCausland(arma::mat& beta_nc_samp,
                                   const arma::vec& y,
                                   const arma::mat& x,
                                   const arma::vec& theta_sr,
                                   const arma::vec& sigma2,
                                   const arma::vec& beta_mean,
                                   arma::vec& m_N,
                                   arma::mat& chol_C_N_inv) {
  typedef void(*Func)(arma::mat&,const arma::vec&,const arma::mat&,const arma::vec&,const arma::vec&,const arma::vec&,arma::vec&,arma::mat&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "sample_beta_McCausland");
  }
  {
    func(beta_nc_samp,
         y,
         x,
         theta_sr,
         sigma2,
         beta_mean,
         m_N,
         chol_C_N_inv);
  }
}

inline void FFBS(arma::mat& beta_nc,
                 const arma::vec& y,
                 const arma::mat& x,
                 const arma::vec& theta_sr,
                 const arma::vec& beta_mean,
                 const arma::vec& sigma2,
                 arma::vec& m_N,
                 arma::mat& chol_C_N_inv) {
  typedef void(*Func)(arma::mat&,const arma::vec&,const arma::mat&,const arma::vec&,const arma::vec&,const arma::vec&,arma::vec&,arma::mat&);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "FFBS");
  }
  {
    func(beta_nc,
         y,
         x,
         theta_sr,
         beta_mean,
         sigma2,
         m_N,
         chol_C_N_inv);
  }
}

inline double do_rgig1(double lambda,
                       double chi,
                       double psi) {
  typedef double(*Func)(double,double,double);
  static Func func = NULL;
  if (func == NULL) {
    func = (Func)R_GetCCallable("shrinkTVP", "do_rgig1");
  }
  {
    return func(lambda,
                chi,
                psi);
  }
}



}


#endif // RCPP_shrinkTVP_H_GEN_
