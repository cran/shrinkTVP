// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include "DG_log_ratio_value.h"
#include "cpp_utilities.h"
using namespace Rcpp;


double DG_MH_step(double current_val,
                  double a_tuning_par,
                  double scale_par,
                  const arma::vec& param_vec,
                  double b,
                  double nu,
                  bool adaptive,
                  arma::vec& batch,
                  double& curr_sd,
                  double target_rate,
                  double max_adapt,
                  int& batch_nr,
                  int batch_size,
                  int& batch_pos){

  double b1 = nu;
  double b2 = nu * b;

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
    sd = a_tuning_par;
  }

  double old_value = current_val;
  double log_prop = R::rnorm(std::log(old_value), sd);
  double proposal = std::exp(log_prop);

  double unif = R::runif(0, 1);

  double log_R = DG_log_ratio_value_marginalBFS(proposal,
                                                old_value,
                                                scale_par,
                                                param_vec,
                                                b1,
                                                b2);

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
