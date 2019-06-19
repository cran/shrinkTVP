#' @export
summary.shrinkTVP_res <- function(object, digits = 3, showprior = TRUE, ...) {

  # Check if digits is scalar, integer and positive
  if (is.scalar(digits) == FALSE |
      is.numeric(digits) == FALSE |
      digits %% 1 != 0 |
      digits < 0){
    stop("digits has to be a single, positive integer")
  }

  # Check if showprior is a logical
  if (bool_input_bad(showprior)){
    stop("showprior has to be a single logical value")
  }

  ret <- attributes(object)
  class(ret) <- c("summary.shrinkTVP_res")
  ret$priorvals <- object$priorvals
  ret$summaries <- object$summaries
  ret$types <- lapply(object, function(x) return(attr(x, "type")))
  ret$digits <- digits
  ret$showprior <- showprior
  ret
}

#' @method print summary.shrinkTVP_res
#' @export
print.summary.shrinkTVP_res <- function(x, ...) {
  cat("\nSummary of ", (x$niter - x$nburn), " MCMC draws after burn-in of ", x$nburn, ".\n", sep = "")

  if(x$showprior == TRUE){
    cat("\nPrior distributions:\n\n")

    # Print prior distributions for beta that react to user choices
    cat(" \u03b2\u2c7c|\u03c4\u00b2\u2c7c \t ~ Normal(0, \u03c4\u00b2\u2c7c)\n")
    cat(" \u03c4\u00b2\u2c7c|a^\u03c4,\u03bb\u00b2 \t ~ Gamma(",
        ifelse(x$learn_a_tau, "a^\u03c4, ", paste0(x$priorvals["a_tau"], ", ")),
        ifelse(x$learn_lambda2,
               ifelse(x$learn_a_tau,
                      "a^\u03c4 \u03bb\u00b2/2",
                      paste0("\u03bb\u00b2", x$priorvals["a_tau"] / 2)),
               ifelse(x$learn_a_tau,
                      paste0( x$priorvals["lambda2"] / 2, "a^\u03c4"),
                      x$priorvals["lambda2"] * x$priorvals["a_tau"] / 2)),
        ")\n",
        sep = "")
    if (x$learn_a_tau == TRUE){
      cat(" a^\u03c4 \t \t ~ Gamma(", x$priorvals["nu_tau"], ", ", x$priorvals["nu_tau"] * x$priorvals["b_tau"], ")\n", sep = "")
    }
    if (x$learn_lambda2 == TRUE){
      cat(" \u03bb\u00b2 \t \t ~ Gamma(", x$priorvals["e1"], ", ", x$priorvals["e2"], ")\n", sep = "")
    }

    cat("\n")

    # Print prior distributions for theta that react to user choices
    cat(" \u221a\u03b8\u2c7c|\u03be\u00b2\u2c7c \t ~ Normal(0, \u03be\u00b2\u2c7c)\n")
    cat(" \u03be\u00b2\u2c7c|a^\u03be,\u03ba\u00b2 \t ~ Gamma(",
        ifelse(x$learn_a_xi, "a^\u03be, ", paste0(x$priorvals["a_xi"], ", ")),
        ifelse(x$learn_kappa2,
               ifelse(x$learn_a_xi,
                      "a^\u03be \u03ba\u00b2/2",
                      paste0(x$priorvals["a_xi"] / 2, "\u03ba\u00b2")),
               ifelse(x$learn_a_xi,
                      paste0(x$priorvals["kappa2"] / 2, "a^\u03be"),
                      x$priorvals["kappa2"] * x$priorvals["a_xi"] / 2)),
        ")\n",
        sep = "")
    if (x$learn_a_xi == TRUE){
      cat(" a^\u03be \t \t ~ Gamma(", x$priorvals["nu_xi"], ", ", x$priorvals["nu_xi"] * x$priorvals["b_xi"], ")\n", sep = "")
    }
    if (x$learn_kappa2 == TRUE){
      cat(" \u03ba\u00b2 \t \t ~ Gamma(", x$priorvals["d1"], ", ", x$priorvals["d2"], ")\n", sep = "")
    }

    # Prior distributions for sigma2
    if (x$sv == TRUE){
      cat("\nPrior distributions on stochastic volatility parameters:\n")
      cat(" \u03bc \t \t ~ Normal(", x$priorvals["bmu"], ", ", sqrt(x$priorvals["Bmu"]), ")\n", sep = "")
      cat(" (\u03c6 + 1)/2 \t ~ Beta(", x$priorvals["a0_sv"], ", ", x$priorvals["b0_sv"], ")\n", sep = "")
      cat(" \u03c3\u209b\u00b2 \t \t ~ Gamma(0.5, ", 0.5 * x$priorvals["Bsigma_sv"], ")\n", sep = "")
    } else {
      cat("\n")
      cat(" \u03c3\u00b2|C\u2080 \t \t ~ InvGamma(", x$priorvals["c0"], ", ", "C\u2080)\n", sep = "")
      cat(" C\u2080 \t \t ~ Gamma(", x$priorvals["g0"], ", ", x$priorvals["G0"], ")\n", sep = "")
    }

    cat("\n")

  }

  # Posterior summaries
  cat("\nStatistics of posterior draws of parameters (thinning = ", x$nthin, "):\n\n", sep = "")

  # The posterior summaries are printed by printing a dataframe where all NA values are
  # set to "", thereby not showing up in the console. This is done to make sure everything
  # is nicely aligned
  ind <- 1
  for (i in x$summaries){
    if (ind == 1){
      post_sum <- rbind(round(i, x$digits), matrix(NA, ncol = ncol(i), nrow = 1))
    } else {
      post_sum <- rbind(post_sum, round(i, x$digits), matrix(NA, ncol = ncol(i), nrow = 1))
    }

    ind <- ind + 1
  }
  post_sum <- as.data.frame(cbind(rownames(post_sum), post_sum), stringsAsFactors = FALSE, row.names = FALSE)
  post_sum[is.na(post_sum)] <- ""

  # Adjust column names
  colnames(post_sum) <- c("param", "mean", "sd", "median", "HPD 2.5%", "HPD 97.5%", "ESS")

  # Remove all automatically generated row names (X.1, X.2, etc.)
  # post_sum$param[grep("X\\.\\d", rownames(post_sum))] <- ""

  # The summaries of theta_sr are based on the absolute value, this has to be reflected in the param name
  post_sum$param <- ifelse(grepl("theta_sr", post_sum$param), paste0("abs(", post_sum$param, ")"), post_sum$param)

  # Finally print it
  print(post_sum, row.names = FALSE, right = FALSE)

  invisible(x)
}


#' Graphical summary of posterior distribution for a time-varying parameter
#'
#' \code{plot.mcmc.tvp} plots empirical posterior quantiles for a time-varying parameter.
#'
#' @param x a mcmc.tvp object
#' @param probs numeric vector of quantiles to plot for each point in time, with values in [0,1].
#' The value 0.5 is added to \code{probs} if the user does not include it.
#' The default value is \code{c(0.025, 0.25, 0.5, 0.75, 0.975)}.
#' @param ... further arguments to be passed to \code{plot}.
#' @return Called for its side effects and returns invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#' data <- sim$data
#'
#' output <- shrinkTVP(y ~ x1 + x2, data)
#' plot(output$beta$beta_x1)
#' }
#' @export
plot.mcmc.tvp <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...){

  if (is.vector(probs) == FALSE | is.list(probs) == TRUE){
    stop("probs has to be a vector")
  }

  # Check if all probs are specified correctly
  if (any(probs > 1) | any(probs < 0) | any(sapply(probs, numeric_input_bad_zer))){
    stop("all elements of probs have to be numbers between >= 0 and <= 1")
  }

  # Add median prob if the user did not
  if (!0.5 %in% probs){
    warning("Adding median to probs", immediate. = TRUE)
    probs <- c(probs, 0.5)
  }
  # Sort by ascending size
  probs <- sort(probs)

  # Calculate all quantiles
  if (length(probs) == 1){
    quants <- matrix()
  }
  quants <- apply(x, 2, quantile, probs)

  # Create x_vec for plotting
  x_vec <- 1:ncol(x)

  # Extract all user specified args and add standard ones that are missing
  args <- list(...)
  standard_args <- list(ylim = c(min(quants), max(quants)), type = "l", xlab = "T", ylab = "")
  missing_args <- names(standard_args)[!names(standard_args) %in% names(args)]
  args[missing_args] <- standard_args[missing_args]

  # Add x and y to arguments
  args$x <- x_vec
  if (length(probs) == 1){
    args$y <- quants
  } else {
    args$y <- quants["50%", ]
  }


  # Plot
  do.call(plot, args)

  # Plot all quantiles as red dashed lines
  for (i in rownames(quants)){
    if (i != "50%"){
      lines(quants[i, ], lty = 2, col = "red")
    }
  }
}

#' Graphical summary of posterior distribution
#'
#' \code{plot.shrinkTVP_res} generates plots visualizing the posterior distribution.
#'
#' @param x a \code{shrinkTVP_res} object.
#' @param pars a character vector containing the names of the parameters to be visualized.
#' The names have to coincide with the names of the list elements of the \code{shrinkTVP_res}
#' object. Throws an error if any element of \code{pars} does not fulfill this criterium.
#' The default is \code{c("beta")}.
#' @param nplot positive integer that indicates the number of tvp plots to display on a single
#' page before a new page is generated. The default value is 3.
#' @param ... further arguments to be passed to the respective plotting functions.
#' @return Called for its side effects and returns invisibly.
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#' data <- sim$data
#'
#' output <- shrinkTVP(y ~ x1 + x2, data)
#' plot(output)
#' }
#'
#' ## Will produce an error because 'hello' is not a parameter in the model
#' \dontrun{
#' plot(output, pars = c("beta", "hello"))
#' }
#'
#' @export
plot.shrinkTVP_res <- function(x, pars = c("beta"), nplot = 3, ...){


  # Check that any pars are specified
  if (length(pars) == 0){
    stop("Please specify pars to plot.")
  }

  # Check if pars is of type "char"
  if (is.character(pars) == FALSE){
    stop("chars has to be of type characater")
  }

  # Check if all specified par values are allowed
  sample_names <- names(x)[sapply(x, function(y) attr(y, "type")) == "sample"]
  if (any(!pars %in% sample_names)){
    stop(paste0("pars can only contain the following values: ", paste(sample_names, collapse = ", ")))
  }

  # Check if nplot is okay
  if (int_input_bad(nplot)){
    stop("nplot has to be a single, positive integer")
  }

  if (nplot == 0){
    stop("nplot cant be 0")
  }


  # Save mfrow value so that it can be restored after plotting
  prev_mfrow <- par()$mfrow


  # Create sublist containing only parameters specified in par and copy attributes
  obj <- x[pars]
  mostattributes(obj) <- attributes(x)
  names(obj) <- names(x[pars])

  # Idea: loop over all list elements and apply plot method specific to that object (mcmc or mcmc.tvp)
  for (i in 1:length(obj)){

    # TVP parameters will be a list, so we differentiate between TVP and non-TVP this way
    if (is.list(obj[[i]]) == TRUE){

      # Plot a max of three objects per plot
      if (length(obj[[i]]) > nplot){

        # Change mfrow
        par(mfrow = c(nplot, 1))

        # Apply plotting method in a loop
        for (j in 1:length(obj[[i]])){

          # Extract all user specified args
          args <- list(...)
          # Augment user specified arguments with standard ones
          standard_args <- list(main = paste(names(obj)[i], "of", attr(obj, "colnames")[j]))
          missing_args <- names(standard_args)[!names(standard_args) %in% names(args)]
          args[missing_args] <- standard_args[missing_args]

          # Add object to be plotted to args list
          args$x <- obj[[i]][[j]]

          # Plot
          do.call(plot, args)

          # Wait for user to move on
          if (j %% nplot == 0 & j != length(obj[[i]])){
            readline("Hit <Return> to see next plot: ")
          }
        }

        # Move to next series of plots if current parameter has been thoroughly plotted
        if (i < length(obj)){
          readline("Hit <Return> to see next plot: ")
        }

      } else {

        # If there are less than nplot parameters to plot, adjust mfrow to nr of parameters
        par(mfrow = c(length(obj[[i]]), 1))
        for (j in 1:length(obj[[i]])){

          # Extract all user specified args
          args <- list(...)
          # Augment user specified arguments with standard ones
          standard_args <- list(main = paste(names(obj)[i], "of", attr(obj, "colnames")[j]))
          missing_args <- names(standard_args)[!names(standard_args) %in% names(args)]
          args[missing_args] <- standard_args[missing_args]

          # Add object to be plotted to args list
          args$x <- obj[[i]][[j]]

          # Plot
          do.call(plot, args)

        }

        if (i < length(obj)){
          readline("Hit <Return> to see next plot: ")
        }

      }


    } else {

      # In non-TVP case, simply use coda's plot method
      plot(obj[[i]], ...)

      if (i < length(obj)){
        readline("Hit <Return> to see next plot: ")
      }
    }
  }

  # Reset mfrow
  par(mfrow = prev_mfrow)
}


# Small convenience function to check if something is a scalar
is.scalar <- function(x) is.atomic(x) && length(x) == 1

# Small input checkers
numeric_input_bad <- function(x) {
  if (is.scalar(x) == TRUE){
    return(is.na(x) | x <= 0 | is.numeric(x) == FALSE )
  } else {
    return(TRUE)
  }
}

numeric_input_bad_zer <- function(x) {
  if (is.scalar(x) == TRUE){
    return(is.na(x) | x < 0 | is.numeric(x) == FALSE )
  } else {
    return(TRUE)
  }
}

numeric_input_bad_ <- function(x) {
  if (is.scalar(x) == TRUE){
    return(is.na(x) | is.numeric(x) == FALSE )
  } else {
    return(TRUE)
  }
}

int_input_bad <- function(x) {
  if (is.scalar(x) == TRUE){
    if (is.numeric(x) == TRUE){
      return(is.na(x) | x < 0 | x %% 1 != 0)
    } else {
      return(TRUE)
    }
  } else {
    return(TRUE)
  }
}

bool_input_bad <- function(x){
  if (is.scalar(x) == TRUE){
    return(is.na(x) | is.logical(x) == FALSE)
  } else {
    return(TRUE)
  }
}
