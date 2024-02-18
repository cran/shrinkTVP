#include <RcppArmadillo.h>
#include <math.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>
#include "cpp_utilities.h"
using namespace Rcpp;

double sample_lambda_iid(double a,
                         double c,
                         double psi) {
  double p1 = a + c;
  double p2 = a/c + 1/psi;

  double lambda_samp = R::rgamma(p1, 1.0/p2);

  res_protector(lambda_samp);

  return lambda_samp;
}

double sample_psi_iid(double c,
                      double lambda,
                      double w) {

  double p1 = c + 0.5;
  double p2 = lambda + std::pow(w, 2) * 0.5;

  double psi_samp = 1/R::rgamma(p1, 1.0/p2);

  res_protector(psi_samp);

  return psi_samp;
}


double sample_lambda_0(double kappa_1,
                       double a,
                       double c,
                       double rho) {
  double p1 = a + kappa_1;
  double p2 = a/c * 1.0/(1.0 - rho);
  double lambda_0_samp = R::rgamma(p1, 1.0/p2);

  res_protector(lambda_0_samp);

  return lambda_0_samp;
}

arma::rowvec sample_lambda(const arma::vec & kappa,
                           const arma::vec & psi,
                           double a,
                           double c,
                           double rho){


  int N = kappa.n_elem;
  arma::vec lambda(N, arma::fill::zeros);

  for (int t = 0; t < (N - 1); t++) {
    double p1 = a + c + kappa(t) + kappa(t + 1);
    double p2 = a/c * (1.0 + rho)/(1.0 - rho) + 1.0/psi(t);
    lambda(t) = R::rgamma(p1, 1.0/p2);
  }

  double p1 = a + c + kappa(N - 1);
  double p2 = a/c * 1.0/(1.0 - rho) + 1.0/psi(N-1);
  lambda(N-1)  = R::rgamma(p1, 1.0/p2);

  std::for_each(lambda.begin(), lambda.end(), res_protector);

  return lambda.t();
}

arma::rowvec sample_psi(const arma::vec& lambda,
                        const arma::vec& beta_nc,
                        double c){

  int N = lambda.n_elem;
  arma::vec res(N, arma::fill::zeros);

  double p1 = c + 0.5;
  double p2 = lambda(0) + std::pow(beta_nc(0), 2) * 0.5;
  res(0) = 1/R::rgamma(p1, 1/p2);

  for (int t = 1; t < N; t++) {
    p2 = lambda(t) + std::pow(beta_nc(t) - beta_nc(t-1), 2) * 0.5;
    res(t) = 1/R::rgamma(p1, 1/p2);
  }

  std::for_each(res.begin(), res.end(), res_protector);


  return res.t();
}

void gen_P(const arma::vec& psi,
           const arma::vec& kappa,
           double a,
           double c,
           double rho,
           arma::uvec index,
           unsigned int batch_no,
           arma::vec& a_t,
           arma::vec& b_t,
           arma::vec& z_t,
           arma::vec& P,
           const arma::vec& pi_t) {

  int T = psi.n_elem;

  for (int t = 0; t < T; t++) {

    if (index(t) == batch_no) {
      if (t == 0) {
        a_t(t) = a;
        b_t(t) = a + c + kappa(1);
        z_t(t) = rho * pi_t(t) * a * psi(t) /
          (a * psi(t) + c * (1 - rho));

        P(t) = 1/gsl_sf_hyperg_2F1(a_t(t), b_t(t), a, z_t(t));

      } else if (t == (T - 1)) {
        a_t(t) = a + c + kappa(t-1);
        b_t(t) = a + c;
        z_t(t) = (1 - pi_t(t-1)) * a * psi(t) /
          (a * psi(t) + c * (1 - rho));

        P(t) = 1/gsl_sf_hyperg_2F1(a_t(t), b_t(t), a, z_t(t));

      } else {
        a_t(t) = a + c + kappa(t-1);
        b_t(t) = a + c + kappa(t+1);
        z_t(t) = (1 - pi_t(t-1)) * pi_t(t) * a * psi(t) /
          (a * psi(t) + c * (1 - rho));

        P(t) = 1/gsl_sf_hyperg_2F1(a_t(t), b_t(t), a, z_t(t));
      }
    }
  }
}

void sample_kappa_fast_marg_alternating(arma::vec& kappa,
                                        const arma::vec& psi,
                                        double a,
                                        double c,
                                        double rho) {

  // Modify function so that it first samples all even entries and then all odd entries
  int T = psi.n_elem;

  // Generate vectors to hold intermediate values
  arma::vec pi_t = (a * psi + c * (1 - rho)) / ((1 + rho) * a * psi + c * (1 - rho));
  arma::vec a_t(T, arma::fill::zeros);
  arma::vec b_t(T, arma::fill::zeros);
  arma::vec z_t(T, arma::fill::zeros);


  // Generate index for even and odd entries, as to sample them separately
  arma::vec index = arma::linspace(0, T-1, T);
  arma::uvec batch_ind = (index - arma::floor(index/2.0)*2.0) == 0.0;

  // Generate one vector of uniform random variables for inverse CDF sampling
  arma::vec U = Rcpp::runif(T, 0, 1);

  for(int i=0; i<2; i++) {

    // Generate vector of initial probabilities
    // This function is indexed, so it only generates them for the current batch
    arma::vec P(T, arma::fill::ones);
    gen_P(psi,
          kappa,
          a,
          c,
          rho,
          batch_ind,
          i,
          a_t,
          b_t,
          z_t,
          P,
          pi_t);

    // Set running index to 0 and reset kappa to 0 for the current batch
    int k = 0;
    kappa(arma::find(batch_ind == i)).fill(0);

    // Generate dummy index variable, just so the while loop runs at least once
    arma::uvec index(T, arma::fill::zeros);
    while(index.n_elem > 0) {
      k += 1;

      // index now contains the indices of the elements that are larger than the probabilities
      // This means that kappa for these entries is at least 1
      index = U > P;
      kappa = kappa + index;

      // index is then overwritten with the indices of the elements that are larger than the probabilities
      index = arma::find(index);

      // Update residual probabilities, only for the current batch
      U(index) = U(index) - P(index);

      // Calculate probabilites for next iteration
      P(index) = P(index) % (a_t(index) + double(k) - 1) % (b_t(index) + k - 1)/(a + double(k) - 1) % z_t(index)/double(k);

      // Realistically, a kappa of 500 is so large that it should not happen all too often
      // Nice thing about this implementation is it only calculates up to the current realisation
      // So it should be quite fast
      if (k > 500) break;
    }


  }


}

