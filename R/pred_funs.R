#' Evaluate the one-step ahead predictive density of a fitted TVP model
#'
#' \code{eval_pred_dens} evaluates the one-step ahead predictive density of a fitted TVP model resulting from a call to
#' shrinkTVP at the points supplied in \code{x}. For details on the approximation of the one-step ahead predictive density used,
#' see the vignette.
#'
#' @param x a real number or a vector of real numbers, taken to be the points at which the predictive density will be evaluated.
#' @param mod an object of class \code{shrinkTVP}, containing the fitted model for which the predictive density should be evaluated.
#' @param data_test a data frame with one row, containing the one-step ahead covariates. The names of the covariates have to match the
#' names of the covariates used during model estimation in the call to \code{shrinkTVP}.
#' @param log a single logical value determining whether the density should be evaluated on the log scale or not.
#'
#' @return The value returned is a vector of length \code{length(x)}, containing the values of the predictive density evaluated
#' at the points supplied in \code{x}.
#'
#' @examples
#' \donttest{
#' # Simulate data
#'set.seed(123)
#'sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#'data <- sim$data
#'
#'# Estimate model
#'res <- shrinkTVP(y ~ x1 + x2, data = data[1:199, ])
#'
#'# Create sequence of x values where the density is to be evaluated
#'x_vals <- seq(0, 12, by = 0.1)
#'
#'# Evaluate density and plot
#'dens <- eval_pred_dens(x_vals, res, data[200, ])
#'plot(x_vals, dens, type = "l")
#'
#'# Add vertical line where true value of the one-step ahead y lies
#'abline(v = data$y[200])
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family prediction functions
#'@export
eval_pred_dens <- function(x, mod, data_test, log = FALSE){

  if (!"shrinkTVP" %in% class(mod)){
    stop("mod has to be of class 'shrinkTVP'")
  }

  if (bool_input_bad(log)){
    stop("log has to be a single logical value")
  }

  if (is.data.frame(data_test) == FALSE || nrow(data_test) > 1){
    stop("data_test has to be a data frame with one row")
  }

  if (is.vector(x) == FALSE || any(sapply(x, numeric_input_bad_))){
    stop("x has to be a vector of numbers")
  }

  nsamp <- nrow(mod$beta_mean)
  N <- length(mod$model$y)
  d <- ncol(mod$beta_mean)

  # For dynamic shrinkage, psis have to be predicted
  if ("shrinkDTVP" %in% class(mod)){

    psi_future <- matrix(NA_real_, nsamp, d)

    a_psi <- attr(mod, "a_psi")
    c_psi <- attr(mod, "c_psi")



    # In the iid case, the psis are drawn from an F distribution
    if (attr(mod, "iid")) {

      for (j in 1:d) {
        if (attr(mod, "shrink_inter") == FALSE & attr(mod, "inter_column") == j) {
          psi_future[,j] <- 1
        } else {
          psi_future[,j] <- rf(nsamp, 2 * a_psi[j], 2 * c_psi[j])
        }
      }

    } else {
      # In the non-iid case, they are simulated from their respective law of motion
      for (j in 1:d) {
        if (attr(mod, "shrink_inter") == FALSE & attr(mod, "inter_column") == j) {
          psi_future[,j] <- 1
        } else {
          kappa_future <- rpois(nsamp, mod$lambda_p[[j]][, N] * a_psi[j]/c_psi[j] * mod$rho_p[,j]/(1 - mod$rho_p[,j]))
          lambda_future <- rgamma(nsamp, a_psi[j] + kappa_future, a_psi[j]/c_psi[j] * 1/(1 - mod$rho_p[,j]))
          psi_future[,j] <- 1/rgamma(nsamp, c_psi[j], lambda_future)
        }
      }
    }
  } else {
    # For static shrinkage, psis are set to 1, as they are essentially a special case of dynamic shrinkage
    psi_future <- matrix(1, nsamp, d)
  }

  terms <- delete.response(mod$model$terms)
  m <- model.frame(terms, data = data_test, xlev = mod$model$xlevels)
  x_test <- model.matrix(terms, m)


  if (any(is.na(x_test))){
    stop("No NA values are allowed in covariates")
  }

  if (nrow(x_test) > 1) {
    stop("Mismatch in names of data_test")
  }


  # Add autoregressive terms
  p <- attr(mod, "p")
  if (p > 0) {
    N <- length(mod$model$y)
    ar_terms <- matrix(mod$model$y[(N - p + 1):N], ncol = p)
    colnames(ar_terms) <- paste0("ar", 1:p)
    x_test <- cbind(x_test, ar_terms)
  }

  colnames(x_test)[colnames(x_test) == "(Intercept)"] <- "Intercept"

  if (attr(mod, "sv") == TRUE){
    sig2 <- mod$sigma2[, ncol(mod$sigma2)]
    sv_phi <- mod$sv_phi
    sv_mu <- mod$sv_mu
    sv_sigma2 <- mod$sv_sigma2
  } else {
    sig2 <- mod$sigma2
    sv_phi = c(1)
    sv_mu = c(1)
    sv_sigma2 = c(1)
  }

  res <- pred_dens_mix_approx_dyn(theta_sr = as.matrix(mod$theta_sr),
                                  beta_mean = as.matrix(mod$beta_mean),
                                  sig2_samp = as.matrix(mod$sigma2),
                                  psi_future = psi_future,
                                  sv = attr(mod, "sv"), sv_phi = sv_phi, sv_mu = sv_mu, sv_sigma2 = sv_sigma2,
                                  x_test = x_test,
                                  y_test = x,
                                  chol_C_N_inv_samp = mod$internals$chol_C_N_inv,
                                  m_N_samp = mod$internals$m_N,
                                  M = nrow(mod$theta_sr), log = log)

  return(as.vector(res))

}


#' Calculate the Log Predictive Density Score for a fitted TVP model
#'
#'
#' \code{LPDS} calculates the one-step ahead Log Predictive Density Score (LPDS) of a fitted TVP model resulting from a call to
#' \code{shrinkTVP} For details on the approximation of the one-step ahead predictive density used, see the vignette.
#'
#' @param mod an object of class \code{shrinkTVP}, containing the fitted model for which the LPDS should be calculated.
#' @param data_test a data frame with one row, containing the one-step ahead covariates and response. The names of the covariates
#' and the response have to match the names used during model estimation in the call to \code{shrinkTVP}.
#'
#' @return A real number equaling the calculated LPDS.
#'
#' @examples
#' \donttest{
#' # Simulate data
#'set.seed(123)
#'sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#'data <- sim$data
#'
#'# Estimate model
#'res <- shrinkTVP(y ~ x1 + x2, data = data[1:199, ])
#'
#'# Calculate LPDS
#' LPDS(res, data[200,])
#' }
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family prediction functions
#'
#'@export
LPDS <- function(mod, data_test){

  if (!"shrinkTVP" %in% class(mod)){
    stop("mod has to be of class 'shrinkTVP'")
  }

  if (is.data.frame(data_test) == FALSE || nrow(data_test) > 1){
    stop("data_test has to be a data frame with one row")
  }


  # Create Vector y
  terms <- mod$model$terms
  m <- model.frame(terms, data = data_test, xlev = mod$model$xlevels)
  y <- model.response(m, "numeric")

  return(as.numeric(eval_pred_dens(y, mod, data_test, log = TRUE)))

}

#' Draw from posterior predictive density of a fitted TVP model
#'
#' \code{forecast_shrinkTVP} draws from the posterior predictive distribution of a fitted TVP model resulting from a call to
#' \code{shrinkTVP}.
#'
#' @param mod an object of class \code{shrinkTVP}, containing the fitted model.
#' @param newdata a data frame containing the future covariates. The names of the covariates
#' have to match the names used during model estimation in the call to \code{shrinkTVP}.
#' @param n.ahead a single, positive integer indicating the forecasting horizon, i.e. how many time-points into the future
#' the posterior predictive distribution should be sampled from. Can not be larger than the number of rows in \code{newdata}.
#'
#' @return The value returned is a list object of class \code{shrinkTVP_forc} containing the samples from the
#' posterior predictive density.
#'
#' @examples
#' \donttest{
#' # Simulate data
#'set.seed(123)
#'sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#'data <- sim$data
#'
#'# Estimate model
#'res <- shrinkTVP(y ~ x1 + x2, data = data[1:190, ])
#'
#'# Forecast
#'forc <- forecast_shrinkTVP(res, data[191:200, ])
#'
#'# Plot
#'plot(forc)
#' }
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family prediction functions
#'
#' @export
forecast_shrinkTVP <- function(mod, newdata, n.ahead){

  # Check if mod is of class shrinkTVP
  if (!inherits(mod, "shrinkTVP")) {
    stop("mod has to be an object of class shrinkTVP")
  }

  if (!missing(n.ahead)) {
    if (int_input_bad(n.ahead) || n.ahead == 0) {
      stop("n.ahead has to be a single, positive integer")
    }
  }

  # Check x is data frame and that either newdata or n.ahead is supplied
  if (!missing(newdata)){
    if (is.data.frame(newdata) == FALSE){
      stop("newdata has to be a data frame")
    }

    if (missing(n.ahead)){
      n.ahead <- nrow(newdata)
    } else {
      if (n.ahead > nrow(newdata)) {
        stop("n.ahead can not be larger than the number of rows in newdata")
      }
    }
  } else {
    # Create a dummy newdata data frame, so the formula interface knows the dimensions
    newdata <- data.frame(dummy = rep(NA, n.ahead))
  }

  terms <- delete.response(mod$model$terms)
  m <- model.frame(terms, data = newdata, xlev = mod$model$xlevels)
  x_test <- model.matrix(terms, m)

  # Number of saved posterior draws, number of covariates, number of observations, number of ar coefficients
  nsamp <- nrow(mod$beta_mean)
  N <- length(mod$model$y)
  d <- ncol(mod$beta_mean)
  p <- attr(mod, "p")


  # Get final beta
  if(is.list(mod$beta)){
    beta_prev <- sapply(mod$beta, function(x){
      return(x[, ncol(x)])
    })
  } else {
    beta_prev <- as.matrix(mod$beta[, ncol(mod$beta)])
  }

  # If dynamic (and not iid), get final lambda and other hyperparameters
  if ("shrinkTVP_dyn" %in% class(mod)){

    a_psi <- attr(mod, "a_psi")
    c_psi <- attr(mod, "c_psi")

    if (!attr(mod, "iid")) {
      if(is.list(mod$lambda_p)){
        lambda_prev <- sapply(mod$lambda_p, function(x){
          return(x[, ncol(x)])
        })
      } else {
        lambda_prev <- as.matrix(mod$lambda_p[, ncol(mod$lambda_p)])
      }
    }
  }

  # Storage mods
  beta_pred <- matrix(NA, nrow = nsamp, ncol = d)
  y_pred <- as.data.frame(matrix(0, ncol = n.ahead, nrow = nsamp))

  if (attr(mod, "sv") == TRUE){
    ht_prev <- log(mod$sigma2[, ncol(mod$sigma2)])
    ht_store <- data.frame(matrix(NA, ncol = n.ahead, nrow = nsamp))
  } else {
    sigma2 <- mod$sigma2
  }

  for (ti in 1:n.ahead){

    for (j in 1:d){

      if ("shrinkTVP_dyn" %in% class(mod)){

        # Forecast psi into the future (except for inter_column if it is not shrunken)
        if (attr(mod, "shrink_inter") == FALSE & attr(mod, "inter_column") == j) {
          psi_pred <- rep(1, nsamp)
        } else {
          if (attr(mod, "iid")) {
            psi_pred <- rf(nsamp, 2*a_psi[j], 2*c_psi[j])
          } else {
            kappa_pred <- rpois(nsamp, lambda_prev[,j] * a_psi[j]/c_psi[j] * mod$rho_p[,j]/(1 - mod$rho_p[,j]))
            lambda_pred <- rgamma(nsamp, a_psi[j] + kappa_pred, a_psi[j]/c_psi[j] * 1/(1 - mod$rho_p[,j]))
            psi_pred <- 1/rgamma(nsamp, c_psi[j], lambda_pred)
            lambda_prev <- lambda_pred
          }
        }
        beta_pred[, j] <- beta_prev[, j] + rnorm(nsamp, 0, abs(mod$theta_sr[, j]) * sqrt(psi_pred))

      } else {
        beta_pred[, j] <- beta_prev[, j] + rnorm(nsamp, 0, abs(mod$theta_sr[, j]))
      }
    }


    if (attr(mod, "sv") == TRUE){
      mean_ht <- mod$sv_mu + mod$sv_phi * (ht_prev - mod$sv_mu)
      ht_pred <- rnorm(nsamp, mean_ht, sqrt(mod$sv_sigma2))
      sigma2 <- exp(ht_pred)

      ht_store[, ti] <- ht_pred
      ht_prev <- ht_pred
    }

    # Add regular regressors (excluding AR components)
    if (! 0 %in% dim(x_test)) {
      y_pred[, ti] <- t(x_test[ti, ] %*% t(beta_pred[, 1:(d - p)]))
    }

    # Add AR components
    if (p != 0) {
      for (curr_p in 1:p){
        # Differentiate between case where the previous y is unknown and where it is known
        curr_prev_t <- ti - curr_p

        if (curr_prev_t <= 0) {
          curr_prev_y <- mod$model$y[length(mod$model$y) + curr_prev_t]
        } else {
          curr_prev_y <- y_pred[, ti - curr_p]
        }

        y_pred[, ti] <- y_pred[, ti] + curr_prev_y * beta_pred[, (d - p) + curr_p]
      }
    }

    # Add error term
    y_pred[, ti] <- y_pred[, ti] + rnorm(nsamp, 0, sqrt(sigma2))

    beta_prev <- beta_pred
  }

  res <- list(y_pred = y_pred, y_orig = mod$model$y)

  if (attr(mod, "sv") == TRUE){
    res$n.ahead <- ht_store
  }

  class(res) <- c("shrinkTVP_forc")
  attr(res, "index") <- attr(mod, "index")

  return(res)
}
