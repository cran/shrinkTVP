// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include "log_ratio_value.h"
#include "sample_parameters.h"
using namespace Rcpp;


double MH_step(double current_val, double c_tuning_par, int d, double scale_par, arma::vec param_vec, double b, double nu,
               double hyp1, double hyp2){

  double b1 = nu;
  double b2 = nu * b;

  double old_value = current_val;
  double log_prop = R::rnorm(std::log(old_value), c_tuning_par);
  double proposal = std::exp(log_prop);

  double unif = R::runif(0, 1);

  double log_R = log_ratio_value_marginalBFS(d, proposal, old_value, scale_par, param_vec, b1, b2);

  double res;
  if (std::log(unif) < log_R){
    res = proposal;
  } else {
    res = old_value;
  }

  res_protector(res);

  return res;

}
