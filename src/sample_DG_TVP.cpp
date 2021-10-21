#include <RcppArmadillo.h>
#include "cpp_utilities.h"
#include "DG_sampling_functions.h"
#include "DG_MH_step.h"
#include <math.h>
using namespace Rcpp;


void sample_DG_TVP(const arma::vec& beta_mean,
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

  // Intermediary storage objects for current adaptive MH batch
  arma::vec curr_batch;
  bool is_adaptive = (arma::sum(adaptive) > 0);

  // Sample a_xi and a_tau with MH
  if (learn_a_xi){
    try {
      if (is_adaptive) {
        curr_batch = batches.col(0);
      }
      a_xi = DG_MH_step(a_xi,
                        a_tuning_par_xi,
                        kappa2_B,
                        theta_sr,
                        beta_a_xi,
                        alpha_a_xi,
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
  }


  if (learn_a_tau){
    try {
      if (is_adaptive) {
        curr_batch = batches.col(1);
      }
      a_tau = DG_MH_step(a_tau,
                         a_tuning_par_tau,
                         lambda2_B,
                         beta_mean,
                         beta_a_tau,
                         alpha_a_tau,
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

  // sample tau2 and xi2
  try {
    DG_sample_local_shrink(tau2,
                           beta_mean,
                           lambda2_B,
                           a_tau);
  } catch(...) {
    tau2.fill(nanl(""));
    if (succesful == true){
      fail = "sample tau2";
      fail_iter = j + 1;
      succesful = false;
    }
  }
  try {
    DG_sample_local_shrink(xi2,
                           theta_sr,
                           kappa2_B,
                           a_xi);
  } catch(...) {
    xi2.fill(nanl(""));
    if (succesful == true){
      fail = "sample xi2";
      fail_iter = j + 1;
      succesful = false;
    }
  }

  // sample kappa2 and lambda2, if the user specified it
  if (learn_kappa2_B){
    try {
      kappa2_B = DG_sample_global_shrink(xi2,
                                         a_xi,
                                         d1,
                                         d2);
    } catch (...) {
      kappa2_B = nanl("");
      if (succesful == true){
        fail = "sample kappa2_B";
        fail_iter = j + 1;
        succesful = false;
      }
    }
  }

  if (learn_lambda2_B){
    try {
      lambda2_B = DG_sample_global_shrink(tau2,
                                          a_tau,
                                          e1,
                                          e2);
    } catch (...) {
      lambda2_B = nanl("");
      if (succesful == true){
        fail = "sample lambda2_B";
        fail_iter = j + 1;
        succesful = false;
      }
    }
  }
}
