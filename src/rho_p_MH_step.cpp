#include <RcppArmadillo.h>
#include <math.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>
#include "cpp_utilities.h"
using namespace Rcpp;

double rho_p_log_ratio_value_marg_oeverything(double proposal,
                                              double old_val,
                                              const arma::vec& psi,
                                              double psi_0,
                                              double a,
                                              double c,
                                              double a_rho,
                                              double b_rho,
                                              double alpha,
                                              double beta){

  int N = psi.n_elem;



  // Precalculation
  double lprop = std::log(proposal);
  double lold = std::log(old_val);
  double lone_minus_prop = std::log(1 - proposal);
  double lone_minus_old = std::log(1 - old_val);

  // Change of variable
  // with prior on (0, b)
  double part1 = lprop - lold +
    std::log(b_rho - proposal) - std::log(b_rho - old_val);

  // Prior on rho
  if ((proposal > b_rho) || (proposal < 0)) return -std::numeric_limits<double>::infinity();
  double part2 = (a_rho * alpha - 1) * (lprop - lold) +
    (beta - 1) * (std::log(std::pow(b_rho, a_rho) - std::pow(proposal, a_rho)) -
    std::log(std::pow(b_rho, a_rho) - std::pow(old_val, a_rho)));



  // Likelihood
  double arg = a*a*proposal*psi(0)*psi_0/((a*psi(0) + c*(1-proposal))*(a*psi_0 + c*(1 - proposal)));
  double part3 = std::log(gsl_sf_hyperg_2F1(a + c, a + c, a, arg)) -
    (a + c) * (std::log(c + a*psi(0)/(1 - proposal)) + std::log(a*psi_0 + c*(1 - proposal))) +
    N*c * std::log(1 - proposal); // This part appears the same in all, hence *N

  for (int t = 1; t < N; t++) {
    arg =  a*a*proposal*psi(t)*psi(t - 1)/((a*psi(t) + c*(1-proposal))*(a*psi(t - 1) + c*(1 - proposal)));
    part3 += std::log(gsl_sf_hyperg_2F1(a + c, a + c, a, arg)) -
      (a + c) * (std::log(c + a*psi(t)/(1 - proposal)) + std::log(a*psi(t-1) + c*(1 - proposal)));
  }

  arg = a*a*old_val*psi(0)*psi_0/((a*psi(0) + c*(1-old_val))*(a*psi_0 + c*(1 - old_val)));
  double part4 = -(std::log(gsl_sf_hyperg_2F1(a + c, a + c, a, arg)) -
                   (a + c) * (std::log(c + a*psi(0)/(1 - old_val)) + std::log(a*psi_0 + c*(1 - old_val))) +
                   N*c * std::log(1 - old_val)); // This part appears the same in all, hence *N

  for (int t = 1; t < N; t++) {
    arg =  a*a*old_val*psi(t)*psi(t - 1)/((a*psi(t) + c*(1-old_val))*(a*psi(t - 1) + c*(1 - old_val)));
    part3 -= std::log(gsl_sf_hyperg_2F1(a + c, a + c, a, arg)) -
      (a + c) * (std::log(c + a*psi(t)/(1 - old_val)) + std::log(a*psi(t-1) + c*(1 - old_val)));
  }

  return part1 + part2 + part3 + part4;

}

double rho_p_MH_step_marg_oeverything(double current_val,
                                      const arma::vec& psi,
                                      double psi_0,
                                      double a,
                                      double c,
                                      double a_rho,
                                      double b_rho,
                                      double alpha,
                                      double beta,
                                      double rho_tuning_par,
                                      bool adaptive,
                                      arma::vec& batch,
                                      double& curr_sd,
                                      double target_rate,
                                      double max_adapt,
                                      int& batch_nr,
                                      int batch_size,
                                      int& batch_pos){

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
    sd = rho_tuning_par;
  }


  double old_value = current_val;

  double logit_prop = R::rnorm(std::log(old_value/(b_rho - old_value)), sd);
  double proposal = b_rho * std::exp(logit_prop)/(1.0 + std::exp(logit_prop));



  double log_R = rho_p_log_ratio_value_marg_oeverything(proposal,
                                                        old_value,
                                                        psi,
                                                        psi_0,
                                                        a,
                                                        c,
                                                        a_rho,
                                                        b_rho,
                                                        alpha,
                                                        beta);




  double unif = R::runif(0, 1);
  double res;
  if (std::log(unif) < log_R){
    res = proposal;
  } else {
    res = old_value;
  }


  if (adaptive){



    batch(batch_pos) = (res == old_value) ? 0 : 1;

    batch_pos = (batch_pos + 1) % batch_size;

  }

  return res;
}

