#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <math.h>
#include <float.h>
using namespace Rcpp;

void res_protector(double& x){
  if (std::abs(x) < DBL_MIN * std::pow(10, 280)){
    double sign = std::copysign(1, x);
    x = DBL_MIN * std::pow(10, 280) * sign;
  }

  if (std::abs(x) > DBL_MAX * std::pow(10, -280)){
    double sign = std::copysign(1, x);
    x = DBL_MAX * std::pow(10, -280) * sign;
  }

  if (std::isnan(x)){
    throw 1;
  }
}


void calc_xi2_tau2(arma::vec& param,
                   const arma::vec& param_til,
                   const arma::vec& loc_shrink_til,
                   double glob_shrink,
                   double c,
                   double a){

  param = 2.0 * param_til * c / (loc_shrink_til * glob_shrink * a);

  if (param.has_inf() || param.has_nan() || !all(param != 0)){
    param = arma::exp(std::log(2) + arma::log(param_til) + std::log(c) - arma::log(loc_shrink_til) - std::log(glob_shrink) - std::log(a));
  }

  std::for_each(param.begin(), param.end(), res_protector);

}


void to_CP(arma::mat& beta,
           const arma::mat& beta_nc,
           const arma::vec& theta_sr,
           const arma::vec& beta_mean) {
  beta = (beta_nc.each_col() % theta_sr).each_col() + beta_mean;
}


void to_NCP(arma::mat& beta_nc,
            const arma::mat& beta,
            const arma::vec& theta_sr,
            const arma::vec& beta_mean) {
  beta_nc = (beta.each_col() - beta_mean).each_col() / theta_sr;
}


arma::mat robust_chol (const arma::mat& V){
  // Note that this function returns the lower triangular form,
  // i.e. cholV * cholV.t() = V;

  arma::mat cholV;
  bool chol_success = arma::chol(cholV, V);

  if (chol_success == false) {
    int max_tries = 1000;
    int num_tries = 1;
    double jitter = 1e-12 * arma::mean(V.diag());
    // Jitter main diagonal to (potentially) solve numerical issues
    while ((chol_success == false) & (num_tries <= max_tries) & !std::isinf(jitter)){
      chol_success = arma::chol(cholV, V + jitter * arma::eye(arma::size(V)));

      jitter *= 1.1;
      num_tries += 1;
    }
  }

  if (chol_success == true) {
    return cholV.t();
  } else {
    return arma::mat(arma::size(V), arma::fill::none);
  }
}

arma::mat robust_chol_nontri (const arma::mat& V){
  // Note that this function is only suitable when the solution can be
  // non lower triangular, as it falls back on an eigen decomposition
  //

  arma::mat cholV;
  bool chol_success = arma::chol(cholV, V);

  // Fall back on eigen decomposition in case of failure of chol
  if (chol_success == false){
    arma::vec vals;
    arma::mat vecs;
    arma::mat V_sym = symmatl(V);
    arma::eig_sym(vals, vecs, V_sym, "std");
    // Set negative eigenvalues to 0 since we know that this is the case for positive, semidefinite matrices
    vals.transform( [](double x) { return (x < 0.0) ? 0.0 : x; } );
    arma::mat sqrtLambda = arma::diagmat(arma::sqrt(vals));
    cholV = vecs * sqrtLambda;
  } else {
    cholV = cholV.t();
  }

  return cholV;
}

arma::mat robust_solve(arma::mat A,
                       arma::mat B) {

  arma::mat res;

  bool solve_success = arma::solve(res, A, B, arma::solve_opts::equilibrate);
  if (!solve_success) {
    // Import Rs solve function
    Environment base = Environment("package:base");
    Function Rsolve = base["solve"];
    NumericMatrix res_R = Rsolve(A, B);
    res = arma::mat(res_R.begin(), res_R.nrow(), res_R.ncol(), false);
  }

  return res;
}



/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Original implementation in R code (R package "Bessel" v. 0.5-3) by        */
/*   Martin Maechler, Date: 23 Nov 2009, 13:39                               */
/*                                                                           */
/* Translated into C code by Kemal Dingic, Oct. 2011.                        */
/*                                                                           */
/* Modified by Josef Leydold on Tue Nov  1 13:22:09 CET 2011                 */
/*                                                                           */
/* Translated into C++ code by Peter Knaus, Mar. 2019.                       */
/*                                                                           */
/*---------------------------------------------------------------------------*/

double unur_bessel_k_nuasympt(double x,
                              double nu,
                              bool islog,
                              bool expon_scaled){

  double M_LNPI = 1.14472988584940017414342735135;

  double z;
  double sz, t, t2, eta;
  double d, u1t,u2t,u3t,u4t;
  double res;


  z = x / nu;

  sz = hypot(1,z);
  t = 1. / sz;
  t2 = t*t;

  if (expon_scaled){
    eta = (1./(z + sz));
  } else {
    eta = sz;
  }

  eta += std::log(z) - log1p(sz);
  u1t = (t * (3. - 5.*t2))/24.;
  u2t = t2 * (81. + t2*(-462. + t2 * 385.))/1152.;
  u3t = t*t2 * (30375. + t2 * (-369603. + t2 * (765765. - t2 * 425425.)))/414720.;
  u4t = t2*t2 * (4465125.
                   + t2 * (-94121676.
                   + t2 * (349922430.
                   + t2 * (-446185740.
                   + t2 * 185910725.)))) / 39813120.;
                   d = (-u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;

                   res = std::log(1.+d) - nu*eta - 0.5*(std::log(2.*nu*sz) - M_LNPI);

                   if (islog){
                     return res;
                   } else {
                     return std::exp(res);
                   }
}


void sample_lin_reg_stab(arma::vec& param_vec,
                         const arma::vec& y,
                         const arma::mat& X,
                         const arma::vec& sigma2,
                         const arma::vec& prior_var) {
  int dim = X.n_cols;
  // Re-writing standard results from linear regression into more numerically stable form (see BFS)
  arma::mat X_til = X.t() * arma::diagmat(1.0/sigma2);
  arma::mat A0_sr = arma::diagmat(arma::sqrt(prior_var));
  arma::mat a = X_til * y; // posterior mean 2nd part
  arma::mat Omega_star = A0_sr * X_til * X * A0_sr + arma::eye(dim, dim);

  // Exploit two different forms, as both can be stable in different cases
  arma::mat A_t; // Posterior variance
  arma::mat A_t_til;
  bool solved = arma::solve(A_t_til, Omega_star, A0_sr);
  if (solved == 1){
    A_t = A0_sr * A_t_til;
  } else {
    arma::mat A0 = arma::diagmat(prior_var);
    A_t = arma::inv(X_til * X + arma::diagmat(1.0/arma::diagvec(A0)));
  }


  arma::mat cholA = robust_chol_nontri(A_t);
  /*
   POSTERIOR
   */
  arma::vec v = rnorm(dim);
  param_vec = A_t * a + cholA * v;
}


// Version with time varying sigma2
void sample_lin_reg_rue(arma::vec& param_vec,
                        const arma::vec& y,
                        const arma::mat& X,
                        const arma::vec& sigma2,
                        const arma::vec& prior_var) {

  int dim = X.n_cols;
  int N = X.n_rows;

  arma::mat X_til = X.t();
  for (int t = 0; t < N; t++) {
    X_til.col(t) *= 1.0/sigma2(t);
  }

  arma::mat XtX = X_til * X;
  arma::vec Xty = X_til * y;

  arma::mat L = arma::chol(XtX + arma::diagmat(1.0/prior_var), "lower");
  arma::mat v = robust_solve(arma::trimatl(L), Xty);
  arma::vec mu = robust_solve(arma::trimatu(L.t()), v);

  arma::vec eps = Rcpp::rnorm(dim, 0, 1);
  param_vec = mu + robust_solve(arma::trimatu(L.t()), eps);
}

// Version with homoscedastic sigma2
// For this version, X'X and X'y can be calculated once and recycled (savings can be quite large)!
// Note that this only works for priors where variance of beta does not depend on sigma2
void sample_lin_reg_rue_homosc(arma::vec& param_vec,
                               const arma::vec& Xty,
                               const arma::mat& XtX,
                               double sigma2,
                               const arma::vec& prior_var) {

  int dim = XtX.n_cols;

  arma::mat L = arma::chol((1/sigma2 * XtX + arma::diagmat(1.0/prior_var)), "lower");
  arma::mat v = robust_solve(arma::trimatl(L), Xty/sigma2);
  arma::vec mu = robust_solve(arma::trimatu(L.t()), v);

  arma::vec eps = Rcpp::rnorm(dim, 0, 1);
  param_vec = mu + robust_solve(arma::trimatu(L.t()), eps);
}

// Note that this only works for priors where variance of beta does not depend on sigma2
void sample_lin_reg_bhat(arma::vec& param_vec,
                         const arma::vec& y,
                         const arma::mat& x,
                         double sigma2,
                         const arma::vec& prior_var){

  //This is the algorithm of Bhattacharya et al. (2016)
  int N = x.n_rows;
  int p = x.n_cols;

  double sigma_inv = 1.0/std::sqrt(sigma2);
  arma::vec alpha = y * sigma_inv;
  arma::mat Phi = x * sigma_inv;

  arma::vec eps = Rcpp::rnorm(p, 0, 1);
  arma::vec u = arma::sqrt(prior_var) % eps;
  arma::vec delta = Rcpp::rnorm(N, 0, 1);
  arma::vec v = Phi * u + delta;

  arma::mat DPhi_t = Phi.t();


  for (int i = 0; i < p; i++) {
    DPhi_t.row(i) *= prior_var(i);
  }

  arma::mat W = Phi * DPhi_t + arma::eye(N, N);


  arma::mat L = arma::chol(W, "lower");
  arma::mat vv = robust_solve(arma::trimatl(L), alpha - v);

  arma::mat w = robust_solve(arma::trimatu(L.t()), vv);
  param_vec = u + DPhi_t * w;
}




double samp_disc_given(arma::vec to_sample,
                       arma::vec probs) {
  return arma::as_scalar(RcppArmadillo::sample(to_sample, 1, 1, probs));
}

