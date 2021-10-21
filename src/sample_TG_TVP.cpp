#include <RcppArmadillo.h>
#include "cpp_utilities.h"
#include "TG_sampling_functions.h"
#include "TG_MH_step.h"
#include <math.h>
using namespace Rcpp;


void sample_TG_TVP(const arma::vec& beta_mean,
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
                   int& fail_iter){

  // Intermediary storage objects for current adaptive MH batch
  arma::vec curr_batch;
  bool is_adaptive = (arma::sum(adaptive) > 0);

  // Sample a_xi/a_tau using MH
  if (learn_a_xi){
    try {
      if (is_adaptive) {
        curr_batch = batches.col(0);
      }
      a_xi = TG_MH_step(a_xi,
                        a_tuning_par_xi,
                        kappa2_B,
                        kappa2_til,
                        theta_sr,
                        beta_a_xi,
                        alpha_a_xi,
                        false,
                        d2,
                        c_xi,
                        a_eq_c_xi,
                        adaptive(0),
                        curr_batch,
                        curr_sds(0),
                        target_rates(0),
                        max_adapts(0),
                        batch_nrs(0),
                        batch_sizes(0),
                        batch_pos(0));
      if (is_adaptive) {
        batches.col(0) = curr_batch;
      }
    } catch(...) {
      a_xi = nanl("");
      if (succesful == true){
        fail = "sample a_xi";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    if (a_eq_c_xi == true) {
      c_xi = a_xi;
    }
  }


  if (learn_a_tau){
    try {
      if (is_adaptive) {
        curr_batch = batches.col(1);
      }
      a_tau = TG_MH_step(a_tau,
                         a_tuning_par_tau,
                         lambda2_B,
                         lambda2_til,
                         beta_mean,
                         beta_a_tau,
                         alpha_a_tau,
                         false,
                         e2,
                         c_tau,
                         a_eq_c_tau,
                         adaptive(1),
                         curr_batch,
                         curr_sds(1),
                         target_rates(1),
                         max_adapts(1),
                         batch_nrs(1),
                         batch_sizes(1),
                         batch_pos(1));
      if (is_adaptive) {
        batches.col(1) = curr_batch;
      }
    } catch(...){
      a_tau = nanl("");
      if (succesful == true){
        fail = "sample a_tau";
        fail_iter = j + 1;
        succesful = false;
      }
    }
  }

  // Sample tau2_til
  try {
    TG_sample_prior_var_til(tau2_til,
                            beta_mean,
                            lambda2_til,
                            lambda2_B,
                            a_tau,
                            c_tau);
  } catch(...) {
    tau2_til.fill(nanl(""));
    if (succesful == true){
      fail = "sample tau2";
      fail_iter = j + 1;
      succesful = false;
    }
  }

  // Sample xi2_til
  try {
    TG_sample_prior_var_til(xi2_til,
                            theta_sr,
                            kappa2_til,
                            kappa2_B,
                            a_xi,
                            c_xi);
  } catch(...) {
    xi2_til.fill(nanl(""));
    if (succesful == true){
      fail = "sample xi2";
      fail_iter = j + 1;
      succesful = false;
    }
  }

  // sample c_xi and c_tau with MH
  if (learn_c_xi && a_eq_c_xi == false){
    try {
      if (is_adaptive) {
        curr_batch = batches.col(2);
      }
      c_xi = TG_MH_step(c_xi,
                        c_tuning_par_xi,
                        kappa2_B,
                        xi2_til,
                        theta_sr,
                        beta_c_xi,
                        alpha_c_xi,
                        true,
                        d2,
                        a_xi,
                        a_eq_c_xi,
                        adaptive(2),
                        curr_batch,
                        curr_sds(2),
                        target_rates(2),
                        max_adapts(2),
                        batch_nrs(2),
                        batch_sizes(2),
                        batch_pos(2));
      if (is_adaptive) {
        batches.col(2) = curr_batch;
      }
    } catch(...) {
      c_xi = nanl("");
      if (succesful == true){
        fail = "sample c_xi";
        fail_iter = j + 1;
        succesful = false;
      }
    }
  } else if (a_eq_c_xi == true) {
    c_xi = a_xi;
  }


  if (learn_c_tau && a_eq_c_tau == false){
    try {
      if (is_adaptive) {
        curr_batch = batches.col(3);
      }
      c_tau = TG_MH_step(c_tau,
                         c_tuning_par_tau,
                         lambda2_B,
                         tau2_til,
                         beta_mean,
                         beta_c_tau,
                         alpha_c_tau,
                         true,
                         e2,
                         a_tau,
                         a_eq_c_tau,
                         adaptive(3),
                         curr_batch,
                         curr_sds(3),
                         target_rates(3),
                         max_adapts(3),
                         batch_nrs(3),
                         batch_sizes(3),
                         batch_pos(3));
      if (is_adaptive) {
        batches.col(3) = curr_batch;
      }
    } catch(...){
      c_tau = nanl("");
      if (succesful == true){
        fail = "sample c_tau";
        fail_iter = j + 1;
        succesful = false;
      }
    }
  } else if (a_eq_c_tau == true) {
    c_tau = a_tau;
  }

  // Have to immediately sample kappa2_til/lambda2_til after the MH steps (as it was marginalized out)
  try{
    TG_sample_local_shrink(kappa2_til,
                           theta_sr,
                           xi2_til,
                           kappa2_B,
                           c_xi,
                           a_xi);
  } catch (...) {
    kappa2_til = nanl("");
    if (succesful == true){
      fail = "sample kappa2_til";
      fail_iter = j + 1;
      succesful = false;
    }
  }

  try{
    TG_sample_local_shrink(lambda2_til,
                           beta_mean,
                           tau2_til,
                           lambda2_B,
                           c_tau,
                           a_tau);
  } catch (...) {
    lambda2_til = nanl("");
    if (succesful == true){
      fail = "sample lambda2_til";
      fail_iter = j + 1;
      succesful = false;
    }
  }

  // sample kappa2_B and lambda2_B, if the user specified it
  if (learn_kappa2_B) {
    try {
      d2 = TG_sample_d2(kappa2_B,
                        a_xi,
                        c_xi);
    } catch(...){
      d2 = nanl("");
      if (succesful == true){
        fail = "sample d2";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    try {
      kappa2_B = TG_sample_global_shrink(xi2_til,
                                         kappa2_til,
                                         theta_sr,
                                         a_xi,
                                         c_xi,
                                         d2,
                                         a_eq_c_xi);
    } catch (...) {
      kappa2_B = nanl("");
      if (succesful == true){
        fail = "sample kappa2";
        fail_iter = j + 1;
        succesful = false;
      }
    }
  }

  if (learn_lambda2_B) {
    try {
      e2 = TG_sample_d2(lambda2_B,
                        a_tau,
                        c_tau);
    } catch(...){
      e2 = nanl("");
      if (succesful == true){
        fail = "sample e2";
        fail_iter = j + 1;
        succesful = false;
      }
    }
    try {
      lambda2_B = TG_sample_global_shrink(tau2_til,
                                          lambda2_til,
                                          beta_mean,
                                          a_tau,
                                          c_tau,
                                          e2,
                                          a_eq_c_tau);
    } catch (...) {
      lambda2_B = nanl("");
      if (succesful == true){
        fail = "sample lambda2";
        fail_iter = j + 1;
        succesful = false;
      }
    }
  }

  // Update xi2/tau2 (needed for sammpling theta_sr and beta_mean)
  calc_xi2_tau2(xi2,
                xi2_til,
                kappa2_til,
                kappa2_B,
                c_xi,
                a_xi);

  calc_xi2_tau2(tau2,
                tau2_til,
                lambda2_til,
                lambda2_B,
                c_tau,
                a_tau);
}
