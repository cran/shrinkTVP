#include <RcppArmadillo.h>
#include <math.h>
#include "TG_log_ratio_value.h"
#include "cpp_utilities.h"
using namespace Rcpp;

double TG_MH_step(double current_val,
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
                  int& batch_pos){

  //parameters of beta distribution
  double b1 = nu;
  //double b2 = nu * b;
  double b2 = b;

  double sd;
  if (adaptive) {

    if (batch_pos == (batch_size - 1)){

      double delta = std::min(max_adapt, 1.0/std::sqrt(batch_nr));
      double acc_rate = arma::accu(batch.rows(0, batch_size - 1))/batch_size;

      if (acc_rate > target_rate) {
        sd = std::exp(std::log(curr_sd) + delta);
      } else {
        sd = std::exp(std::log(curr_sd) - delta);
      }

      curr_sd = sd;
      batch_nr += 1;

    } else {
      sd = curr_sd;
    }

  } else {
    sd = c_tuning_par;
  }


  double old_value = current_val;

  double logit_prop = R::rnorm(std::log(old_value/(0.5 - old_value)), sd);
  double proposal = 0.5* std::exp(logit_prop)/(1.0 + std::exp(logit_prop));

  double unif = R::runif(0, 1);

  double log_R;

  if (is_c){
    log_R = TG_log_ratio_value_tg(proposal,
                                  old_value,
                                  scale_par,
                                  scale_vec,
                                  param_vec,
                                  scale_scale,
                                  other_hyp,
                                  b1,
                                  b2);
  } else {
    log_R = TG_log_ratio_value_marginalBFS(proposal,
                                           old_value,
                                           scale_par,
                                           scale_vec,
                                           param_vec,
                                           scale_scale,
                                           other_hyp,
                                           b1,
                                           b2,
                                           common_shrink_par);
  }

  double res;
  if (std::log(unif) < log_R){
    res = proposal;
  } else {
    res = old_value;
  }

  res_protector(res);

  if (adaptive){
    batch(batch_pos) = (res == old_value) ? 0 : 1;
    batch_pos = (batch_pos + 1) % batch_size;
  }

  return res;

}
