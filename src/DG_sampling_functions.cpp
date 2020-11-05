#include <RcppArmadillo.h>
#include "do_rgig1.h"
#include "cpp_utilities.h"
#include <math.h>
using namespace Rcpp;

void DG_sample_local_shrink(arma::vec& local_shrink,
                            const arma::vec& param_vec,
                            double global_shrink,
                            double a){
  int d = local_shrink.n_elem;

  arma::vec param_vec2 = arma::pow(param_vec, 2);

  double p1 = a - 0.5;
  double p2 = a * global_shrink;

  for (int j = 0; j < d; j++){
    double p3 = param_vec2(j);

    local_shrink(j) = do_rgig1(p1, p3, p2);
  }

  std::for_each(local_shrink.begin(), local_shrink.end(), res_protector);
}


double DG_sample_global_shrink(const arma::vec& prior_var,
                               double a,
                               double hyper1,
                               double hyper2) {

  int d = prior_var.n_elem;

  double hyper1_full = hyper1 + a * d;
  double hyper2_full = hyper2 + arma::mean(prior_var) * a * d * 0.5;

  double res = R::rgamma(hyper1_full, 1.0/hyper2_full);

  res_protector(res);
  return res;
}
