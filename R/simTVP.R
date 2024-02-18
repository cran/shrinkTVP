#' Generate synthetic data from a time-varying parameter model
#'
#' \code{simTVP} generates synthetic data from a time-varying parameter model. The covariates are always
#' generated i.i.d. from a Normal(0,1) distribution.
#'
#' @param N integer > 2. Indicates the length of the time series to be
#' generated. The default value is 200.
#' @param d positive integer. Indicates the number of covariates to simulate.
#' The default value is 3.
#' @param sv logical value. If set to \code{TRUE}, the data will be generated with
#' stochastic volatility for the errors of the observation equation using \code{\link{svsim}}. The default value is FALSE.
#' @param DTG logical value. If set to \code{TRUE}, the betas will be generated as dynamic triple gamma processes.
#' The default value is FALSE.
#' @param sigma2 positive real number. Determines the variance of the errors of the observation
#' equation. Ignored if sv is \code{TRUE}. The default value is 1.
#' @param theta \emph{(optional)} vector containing positive real numbers. If supplied, these determine the variances of the
#' innovations of the state equation. Otherwise, the elements of \code{theta} are generated from a X^2(1) distribution.
#' Has to be of length \code{d} or an error will be thrown.
#' @param beta_mean \emph{(optional)} vector containing real numbers. If supplied, these determine the mean of the
#' initial value of the state equation. Otherwise, the elements of \code{beta_mean} are generated from a Normal(0,1) distribution.
#' Has to be of length \code{d} or an error will be thrown.
#' @param sv_mu real number. Determines the mean of the logarithm of the volatility. Ignored if \code{sv} is \code{FALSE}.
#' The default value is 0.
#' @param sv_phi real number between -1 and 1. Determines the persistence of the SV process. Ignored if \code{sv} is \code{FALSE}.
#' The default value is 0.98.
#' @param sv_sigma2 positive, real number. Determines the variance of the innovations of the logarithm of the volatility.
#' Ignored if \code{sv} is \code{FALSE}. The default value is 0.2.
#' @param a_psi positive, real number. Determines the pole parameter of the dynamic triple gamma process.
#' Ignored if \code{DTG} is \code{FALSE}. The default value is 0.5.
#' @param c_psi positive, real number. Determines the tail parameter of the dynamic triple gamma process.
#' Ignored if \code{DTG} is \code{FALSE}. The default value is 2.
#' @param rho real number between 0 and 1. Determines the persistence of the dynamic triple gamma process.
#' Ignored if \code{DTG} is \code{FALSE}. The default value is 0.9
#'
#' @return The value returned is a list object containing:
#' \item{\code{data}}{data frame that holds the simulated data.}
#' \item{\code{true_vals}}{list object containing:
#'   \itemize{
#'   \item \code{theta}: the values of theta used in the data generating process.
#'   \item \code{beta_mean}: the values of beta_mean used in the data generating process.
#'   \item \code{beta}: the true paths of beta used for the data generating process.
#'   \item \code{sigma2}: the value(s) of sigma2 used in the data generating process.
#'   \item \code{lambda_p}: the true paths of lambda_p used for the data generating process. Not returned if DTG is \code{FALSE}.
#'   \item \code{lambda_p_0}: the values of lambda_p_0 used for the data generating process. Not returned if DTG is \code{FALSE}.
#'   \item \code{kappa_p}: the true paths of kappa_p used for the data generating process. Not returned if DTG is \code{FALSE}.
#'   \item \code{psi}: the true paths of psi used for the data generating process. Not returned if DTG is \code{FALSE}.
#'   }
#' }
#'
#' @examples
#' # Generate a time series of length 300
#' res <- simTVP(N = 300)
#'
#' # Extract the generated data
#' data <- res$data
#'
#' # Now with stochastic volatility
#' res_sv <- simTVP(N = 300, sv = TRUE)
#'
#' # Now with dynamic triple gamma process
#' res_DTG <- simTVP(N = 300, DTG = TRUE, c_psi = 1)
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @export
simTVP <- function(N = 200, d = 3, sv = FALSE, DTG = FALSE, sigma2 = 1, theta, beta_mean, sv_mu = 0, sv_phi = 0.98, sv_sigma2 = 0.2, a_psi = 0.5, c_psi = 2, rho = 0.9){

  # Check all inputs
  if (int_input_bad(N)){
    stop("N has to be a single, positive integer")
  }

  if (N < 2){
    stop("N can not be smaller than 2")
  }

  if (int_input_bad(d)){
    stop("d has to be a single, positive integer")
  }

  if (d == 0){
    stop("d can not be 0")
  }

  if (bool_input_bad(sv)){
    stop("sv has to be a single logical value")
  }

  # Check inputs associated with dynamic triple gamma prior
  if (bool_input_bad(DTG)){
    stop("DTG has to be a single logical value")
  }

  if (DTG == TRUE){
    if (numeric_input_bad(a_psi)){
      stop("a_psi has to be a positive real number")
    }

    if (numeric_input_bad(c_psi)){
      stop("c_psi has to be a positive real number")
    }

    if (numeric_input_bad_(rho) || rho < 0 || rho > 1){
      stop("rho has to be a real number between 0 and 1")
    }
  }

  if (sv == FALSE){
    if (numeric_input_bad(sigma2)){
      stop("sigma2 has to be a positive real number")
    }
  } else {
    if (numeric_input_bad_(sv_mu)) {
      stop("sv_mu has to be a real number")
    }

    if (numeric_input_bad_(sv_phi) || abs(sv_phi) >= 1 ) {
      stop("sv_phi has to be a real number between -1 and 1")
    }

    if (numeric_input_bad(sv_sigma2)) {
      stop("sv_sigma2 has to be a positive real number")
    }
  }


  # Generate matrix x
  x <- matrix(NA, N, d)
  x[,1]<- rep(1,N)

  if (d > 1){
    x[, 2:d] <- matrix(rnorm(N * (d - 1)), N, d - 1)
  }

  # Take user supplied beta_mean/theta or generate them
  if (missing(beta_mean)){
    beta_mean <- rnorm(d, 0, 1)
  } else {
    if (is.vector(beta_mean) == FALSE){
      stop("beta_mean has to be a vector")
    }

    if (length(beta_mean) != d){
      stop("beta_mean has to be of length d")
    }
    if (any(sapply(beta_mean, numeric_input_bad_))){
      stop("all elements of beta_mean have to be single numbers")
    }
  }

  if (missing(theta)){
    theta <- rnorm(d, 0, 1)^2
  } else {

    if (is.vector(theta) == FALSE){
      stop("theta has to be a vector")
    }

    if (length(theta) != d){
      stop("theta has to be of length d")
    }
    if (any(sapply(theta, numeric_input_bad_zer))){
      stop("all elements of theta have to be scalars that are >= 0")
    }
  }

  # If dynamic, simulate from DTG prior
  if (DTG){
    kappa_p <- matrix(NA_real_, d, N)
    lambda_p <- matrix(NA_real_, d, N)

    lambda_p_0 <- rgamma(d, a_psi, a_psi/c_psi)
    kappa_p[, 1] <- rpois(d, a_psi/c_psi * rho/(1 - rho) * lambda_p_0)
    lambda_p[, 1] <- rgamma(d, a_psi + kappa_p[, 1], a_psi/c_psi * 1/(1 - rho))

    for(t in 2:N){
      kappa_p[, t] <- rpois(d, a_psi/c_psi * rho/(1 - rho) * lambda_p[, t-1])
      lambda_p[, t] <- rgamma(d, a_psi + kappa_p[, t], a_psi/c_psi * 1/(1 - rho))
    }
    psi0 <- 1/rgamma(d, c_psi, lambda_p_0)
    psi <- matrix(1/rgamma(d*N, c_psi, lambda_p), d, N)
  } else {
    psi0 <- rep(1, d)
    psi <- matrix(1, d, N)
  }

  # Simulate intial state from multivariate normal and simulate states forward

  beta <- matrix(NA, d, N)
  beta_0 <- rnorm(d, beta_mean, sqrt(theta*psi0))
  beta[, 1] <- beta_0 + rnorm(d, 0, sqrt(theta*psi[, 1])) #beta_init + cholQ %*% rnorm(d)
  for(t in 2:N){
    beta[, t]<-  beta[, t-1] + rnorm(d, 0, sqrt(theta*psi[, t])) #cholQ %*% rnorm(d)
  }

  # Create noise depending on sv input, use svsim from stochvol
  if (sv){
    sigma2 <- stochvol::svsim(N, mu = sv_mu, phi = sv_phi, sigma = sqrt(sv_sigma2))$vol^2
  }
  err <- rnorm(N, 0, sqrt(sigma2))

  # Simulate observation equation
  y <- rep(NA, N)
  for(t in 1:N){
    y[t] <- x[t, ]%*%beta[, t] + err[t]
  }

  # Create output
  if (d == 1){
    colnames(x) <- "Intercept"
  } else {
    colnames(x) <- c("Intercept", paste0("x", 1:(ncol(x)-1)))
  }

  res <- list()
  res$data <- data.frame(y, x)

  res$true_vals <- list()
  res$true_vals$theta <- theta
  res$true_vals$beta_mean <- beta_mean
  res$true_vals$beta <- beta
  res$true_vals$sigma2 <- sigma2

  if (DTG){
    res$true_vals$lambda_p <- lambda_p
    res$true_vals$lambda_p_0 <- lambda_p_0
    res$true_vals$kappa_p <- kappa_p
    res$true_vals$psi <- psi
  }

  return(res)

}
