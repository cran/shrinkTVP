#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>
#include "cpp_utilities.h"
using namespace Rcpp;

//this is the log ratio value integrating out the xi2_til
double TG_log_ratio_value_marginalBFS(double proposal,
                                      double old_val,
                                      double scale_par,
                                      const arma::vec& scale_vec,
                                      const arma::vec& param_vec,
                                      double scale_scale,
                                      double c,
                                      double b1,
                                      double b2,
                                      bool common_shrink_par) {

  int d = param_vec.n_elem;

  arma::vec param_vec2 = arma::pow(param_vec, 2);

  // Change of variable
  // with new prior on (0, 0.5)
  double part1 = std::log(proposal) - std::log(old_val) +
    std::log(0.5 - proposal) - std::log(0.5 - old_val);

  // Prior on a
  // beta prior
  double part2 = (b1 - 1) * (std::log(2*proposal) - std::log(2*old_val)) +
    (b2 - 1) * (std::log(1 -2*proposal) - std::log(1 - 2*old_val));

  // kappa2 part
  double part3 = -(R::lbeta(proposal, c) - R::lbeta(old_val, c)) +
    (proposal*(std::log(proposal) + std::log(scale_par/(2.0*c))) - old_val*(std::log(old_val) + std::log(scale_par/(2.0*c)))) -
    (std::log(proposal) - std::log(old_val)) - ((proposal + c)*std::log(1 + proposal*scale_par/(2.0*c)) -
    (old_val + c)*std::log(1 + old_val*scale_par/(2.0*c)));

  // kappa2_til_j part only if a and c are equal
  double parteq;
  if(common_shrink_par == TRUE){
    parteq = -d*(std::lgamma(proposal + 1)- std::lgamma(old_val + 1)) + d*(std::log(proposal) - std::log(old_val)) +
      arma::sum(arma::log(scale_vec))*(proposal - old_val);
  }else{
    parteq = 0;
  }

  // First part of prior on sqrt(theta)_j/beta_j
  double part4 = (-d * std::log(2) + d * 0.5 * std::log(scale_par) - d * 0.5 * std::log(c) +
                  0.5 * arma::sum(arma::log(scale_vec)) + 0.5 * arma::sum(arma::log(param_vec2))) * (proposal - old_val);

  // Second part of prior on sqrt(theta)_j/beta_j
  double part5 = 5 * 0.25 * d * (std::log(proposal) - std::log(old_val)) + d * 0.5 * (proposal * std::log(proposal) -
                                 old_val * std::log(old_val)) - d * (std::lgamma(proposal + 1) - std::lgamma(old_val + 1));

  // Create final part of prior on sqrt(theta)_j/beta_j
  arma::vec besselK_vals_proposal(d, arma::fill::zeros);
  arma::vec besselK_vals_old_val(d, arma::fill::zeros);

  double bessel_nu_proposal = std::abs(proposal - 0.5);
  double bessel_nu_old_val = std::abs(old_val - 0.5);

  double bessel_arg_proposal;
  double bessel_arg_old_val;

  for (int j = 0; j < d; j++){
    bessel_arg_proposal = std::exp(0.5*std::log(proposal) - 0.5*std::log(c) + 0.5*std::log(scale_par)
                                     + 0.5 * std::log(scale_vec(j)) + std::log(std::abs(param_vec(j))));
    bessel_arg_old_val = std::exp(0.5*std::log(old_val) - 0.5*std::log(c) + 0.5*std::log(scale_par)
                                    + 0.5 * std::log(scale_vec(j)) + std::log(std::abs(param_vec(j))));

    if(bessel_nu_proposal < 50 and bessel_arg_proposal < 50){
      besselK_vals_proposal(j) = std::log(R::bessel_k(bessel_arg_proposal, bessel_nu_proposal, true)) - bessel_arg_proposal;
    }else{
      besselK_vals_proposal(j) = unur_bessel_k_nuasympt(bessel_arg_proposal, bessel_nu_proposal, true, false);
    }

    if(bessel_nu_old_val < 50 and bessel_arg_old_val < 50){
      besselK_vals_old_val(j) = std::log(R::bessel_k(bessel_arg_old_val, bessel_nu_old_val, true)) - bessel_arg_old_val;
    }else{
      besselK_vals_old_val(j) = unur_bessel_k_nuasympt(bessel_arg_old_val, bessel_nu_old_val, true, false);
    }
  }

  double part6 = arma::sum(besselK_vals_proposal) - arma::sum(besselK_vals_old_val);

  double res = part1 + part2 + part3 + part4 + part5 + part6 + parteq;

  return(res);
}


double TG_log_ratio_value_tg(double proposal,
                             double old_val,
                             double scale_par,
                             const arma::vec& scale_vec,
                             const arma::vec& param_vec,
                             double scale_scale,
                             double a,
                             double b1,
                             double b2){

  int d = param_vec.n_elem;

  arma::vec param_vec2 = arma::pow(param_vec, 2);

  // Change of variable
  double part1 = std::log(proposal) - std::log(old_val) +  std::log(0.5 - proposal) -
    std::log(0.5 - old_val);

  // Prior on a
  // beta prior
  double part2 = (b1 - 1) * (std::log(2*proposal) - std::log(2*old_val)) +
    (b2 - 1) * (std::log(1 -2*proposal) - std::log(1 - 2*old_val));

  // Likelihood
  double part3 = d * (std::lgamma(proposal + 0.5) - std::lgamma(old_val + 0.5)) -
    d * (std::lgamma(proposal + 1) - std::lgamma(old_val + 1)) +
    d * 0.5 * (std::log(proposal) - std::log(old_val)) -
    (proposal + 0.5) * arma::sum(arma::log(param_vec2 * scale_par*a  + 4 * proposal * scale_vec) -  arma::log(4 * proposal * scale_vec)) +
    (old_val + 0.5) * arma::sum(arma::log(param_vec2 * scale_par*a  + 4 * old_val * scale_vec) - arma::log(4 * old_val * scale_vec));


  double part4 = - (R::lbeta(a, proposal) - R::lbeta(a, old_val))  -
    (a-1)*(std::log(proposal) - std::log(old_val)) -
    ((a + proposal)*(std::log(1 + a*scale_par/(2.0*proposal))) - (a+ old_val)*(std::log(1 + a*scale_par/(2.0*old_val))));


  double res = part1 + part2 + part3 + part4;
  return res;
}
