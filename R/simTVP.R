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
#' @param sigma2 positive real number. Determines the variance on the errors of the observation
#' equation. Ignored if sv is \code{TRUE}. The default value is 1.
#' @param theta \emph{(optional)} vector containing positive real numbers. If supplied, these determine the variances of the
#' innovations of the state equation. Otherwise, the elements of \code{theta} are generated from a ChiSq(1) distribution.
#' Has to be of length \code{d} or an error will be thrown.
#' @param beta_mean \emph{(optional)} vector containing real numbers. If supplied, these determine the mean of the
#' initial value of the state equation. Otherwise, the elements of \code{beta_mean} are generated from a Normal(0,1) distribution.
#' Has to be of length \code{d} or an error will be thrown.
#'
#' @return The value returned is a list object containing:
#' \item{\code{data}}{a data frame that holds the simulated data.}
#' \item{\code{true_vals}}{a list object containing:
#'   \itemize{
#'   \item \code{theta}: the values of theta used in the data generating process.
#'   \item \code{beta_mean}: the values of beta_mean used in the data generating process.
#'   \item \code{beta}: the true paths of beta used for the data generating process.
#'   \item \code{sigma2}: the value(s) of sigma2 used in the data generating process.
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
#' @export
simTVP <- function(N = 200, d = 3, sv = FALSE, sigma2 = 1, theta, beta_mean){

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

  if (sv == FALSE){
    if (numeric_input_bad(sigma2)){
      stop("sigma2 has to be a positive scalar")
    }
  }


  # Generate matrix x
  x <- matrix(NA, N, d)
  x[,1]<- rep(1,N)

  if (d > 1){
    for(i in 2:d){
      x[,i] <- rnorm(N)
    }
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

  # Simulate intial state from multivariate normal and simulate states forward

  if (d == 1){
    Q <- matrix(theta)
  } else {
    Q <- diag(theta)
  }

  cholQ <- suppressWarnings(chol(Q, pivot = TRUE))
  pivot <- attr(cholQ, "pivot")
  oo <- order(pivot)
  cholQ <- t(Q[, oo])

  beta <- matrix(NA, d, N)
  beta_init <- beta_mean + cholQ %*% rnorm(d)
  beta[, 1] <- beta_init + cholQ %*% rnorm(d)
  for(t in 2:N){
    beta[, t]<- beta[, t-1] + cholQ %*% rnorm(d)
  }

  # Create noise depending on sv input, use svsim from stochvol
  if (sv){
    sigma2 <- stochvol::svsim(N, mu = 0)$vol
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

  return(res)

}
