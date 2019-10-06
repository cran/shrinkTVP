// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <math.h>
using namespace Rcpp;

//[[Rcpp::export]]
arma::vec pred_dens_mix_approx(arma::vec x_test, arma::vec y_test,
                               arma::mat theta_sr, arma::mat beta_mean, arma::vec sig2_samp,
                               bool sv, arma::vec sv_phi, arma::vec sv_mu, arma::vec sv_sigma2,
                               arma::cube chol_C_N_inv_samp, arma::cube m_N_samp,
                               int M, bool log){


  int nobs = y_test.n_elem;
  arma::mat F(arma::size(x_test), arma::fill::zeros);
  arma::mat LF(arma::size(x_test), arma::fill::zeros);
  double S;
  double mu;
  arma::vec sig2_pred(M, arma::fill::zeros);

  arma::vec log_dens(M);
  double max;

  arma::vec res(nobs);


  if (sv == false){
    sig2_pred = sig2_samp;
  } else {
    for (int m = 0; m < M; m++){
      sig2_pred(m) = std::exp(R::rnorm(sv_mu(m) + sv_phi(m) * (std::log(sig2_samp(m)) - sv_mu(m)), std::sqrt(sv_sigma2(m))));
    }
  }



  // Loop over observations (if there are more than 1)
  for (int j = 0; j < nobs; j++){
    for (int m = 0; m < M; m++){




      // Create mean and sd
      F = x_test % (theta_sr.row(m)).t();
      LF = arma::solve(arma::trimatl(chol_C_N_inv_samp.slice(m)), F);
      S = arma::as_scalar(LF.t() * LF + F.t() * F + sig2_pred(m));
      mu = arma::as_scalar(x_test.t() * (beta_mean.row(m)).t() + F.t() * m_N_samp.slice(m));

      // Calculate analytical log densities for stability trick
      log_dens(m) = -0.5 * (std::log(2*M_PI) + std::log(S) + std::pow((y_test(j) - mu), 2)/S);

    }

    // Use Sylvias numerical trick to stabilize
    max = log_dens.max();
    res(j) =  -std::log(M) + max + std::log(arma::as_scalar(arma::sum(arma::exp(log_dens - max))));

  }

  if (log == true){
    return(res);
  } else {
    return(arma::exp(res));
  }


}
//sum_pred_y += R::dnorm(y_test, mu, std::sqrt(S), 0);
