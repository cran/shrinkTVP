// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>
#include "unur_bessel_k_nuasympt.h"
using namespace Rcpp;

double log_ratio_value_marginalBFS(int d, double proposal, double old_val, double scale_par, arma::vec param_vec, double b1, double b2) {
  arma::vec besselKvalue_partA(d, arma::fill::none);
  arma::vec besselKvalue_partB(d, arma::fill::none);
  long double par1A = std::abs(proposal - 0.5);
  long double par1B = std::abs(old_val - 0.5);
  for (int j = 0; j < d; j++){
    double par2A = std::exp(0.5*std::log(proposal) + 0.5*std::log(scale_par) + std::log(std::abs(param_vec(j))));
    double par2B = std::exp(0.5*std::log(old_val) + 0.5*std::log(scale_par) + std::log(std::abs(param_vec(j))));


    if(par1A < 50 and par2A < 50){
      besselKvalue_partA(j) = std::log(R::bessel_k(par2A, par1A, true)) - par2A;
    }else{
      besselKvalue_partA(j) = unur_bessel_k_nuasympt(par2A, par1A, true, false);
    }

    if(par1B < 50 and par2B < 50){
      besselKvalue_partB(j) = std::log(R::bessel_k(par2B, par1B, true)) - par2B;
    }else{
      besselKvalue_partB(j) = unur_bessel_k_nuasympt(par2B, par1B, true, false);
    }

  }

  //gamma prior
  double partA = (b1 - 1 + 1 + d/4.0)*(std::log(proposal) - std::log(old_val));
  double partB = (d/2.0*std::log(scale_par) - d*std::log(2) -
                  b2 + arma::as_scalar(arma::sum(arma::log(arma::abs(param_vec)))))*(proposal - old_val);
  double partC = d/2.0*(std::log(proposal)*proposal - std::log(old_val)*old_val);
  double partD =  - d*(std::lgamma(proposal + 1) - log(proposal) - std::lgamma(old_val +1) + log(old_val));
  double partE = arma::as_scalar(arma::sum(besselKvalue_partA -besselKvalue_partB));
  double res = partA + partB + partC +partD + partE;

  return res;
}

