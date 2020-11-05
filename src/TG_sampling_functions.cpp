#include <RcppArmadillo.h>
#include "do_rgig1.h"
#include "cpp_utilities.h"
#include <math.h>
using namespace Rcpp;


void TG_sample_prior_var_til(arma::vec& prior_var_til,
                             const arma::vec& param_vec,
                             const arma::vec& local_shrink,
                             double global_shrink,
                             double a,
                             double c) {

  int d = param_vec.n_elem;
  arma::vec param_vec2 = arma::pow(param_vec, 2);

  double p1 = a - 0.5;
  // a/c in prior_var
  double p2 = 2;
  //double p2 = 2 * a;

  for (int j = 0; j < d; j++){
    double p3 = local_shrink(j) * global_shrink * param_vec2(j) * 0.5 * a/c;

    //double p3 = local_shrink(j) * global_shrink * param_vec2(j)/2;
    double res = do_rgig1(p1, p3, p2);

    res_protector(res);

    prior_var_til(j) = res;
  }
}

double TG_sample_global_shrink(const arma::vec& prior_var_til,
                               const arma::vec& local_shrink,
                               const arma::vec& param_vec,
                               double a,
                               double c,
                               double hyper2,
                               bool common_shrink_par){
  int d = prior_var_til.n_elem;
  arma::vec param_vec2 = arma::pow(param_vec, 2);

  double hyper1_full = a + d * 0.5;
  // a/c in prior variance
  // double hyper2_full = hyper2 + 0.25 * arma::sum(param_vec2 % local_shrink / prior_var_til);
  double hyper2_full = hyper2 + 0.25 * a/c * arma::sum(param_vec2 % local_shrink / (prior_var_til));

  if (std::isnan(hyper2_full) | std::isinf(hyper2_full)){
    hyper2_full = hyper2 + 0.25 * a/c * arma::sum(arma::exp(arma::log(param_vec2) + arma::log(local_shrink) - arma::log(prior_var_til)));
  }

  double global_shrink = R::rgamma(hyper1_full, 1.0/hyper2_full);

  res_protector(global_shrink);
  return(global_shrink);
}


void TG_sample_local_shrink(arma::vec& local_shrink,
                            const arma::vec& param_vec,
                            const arma::vec& prior_var_til,
                            double global_shrink,
                            double c,
                            double a){

  int d = local_shrink.n_elem;
  arma::vec param_vec2 = arma::pow(param_vec, 2);

  double a_full = c + 0.5;

  for (int j = 0; j < d; j++){
    // a/c in prior variance
    double b_full = 1 + global_shrink * param_vec2(j) * a / (4 * prior_var_til(j) * c);
    // double b_full = c + global_shrink * param_vec2(j) / (4 * prior_var_til(j));
    local_shrink(j) = R::rgamma(a_full, 1.0/b_full);
  }

  std::for_each(local_shrink.begin(), local_shrink.end(), res_protector);
}

double TG_sample_d2(double global_shrink,
                    double a,
                    double c){

  double p1 = a + c;
  double res = R::rgamma(p1, 1.0/(global_shrink + 2*c/a));

  res_protector(res);
  return(res);
}


