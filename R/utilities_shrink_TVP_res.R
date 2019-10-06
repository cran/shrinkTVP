
#' Nicer printing of shrinkTVP objects
#'
#' @param x A \code{shrinkTVP} object
#' @param ... Currently ignored.
#'
#' @return Called for its side effects and returns invisibly.
#' @export
print.shrinkTVP <- function(x, ...){
  ind <- attr(x, "index")
  cat(paste0("Object containing a fitted TVP model ", ifelse(attr(x, "sv"), "with stochastic volatility ", ""), "with:\n",
             " - ", formatC(length(attr(x, "colnames")), width = 7), " covariates\n",
             " - ", formatC(length(x$model$y), width = 7), " timepoints, running from ", min(ind), " to ", max(ind), "\n",
             " - ", formatC(attr(x, "niter"), width = 7), " MCMC draws\n",
             " - ", formatC(attr(x, "nburn"), width = 7), " burn-n\n",
             " - ", formatC(attr(x, "nthin"), width = 7), " thinning"))
  invisible(x)
}

#' @export
summary.shrinkTVP <- function(object, digits = 3, showprior = TRUE, ...) {

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
  class(ret) <- c("summary.shrinkTVP")
  ret$priorvals <- object$priorvals
  ret$summaries <- object$summaries
  ret$types <- lapply(object, function(x) return(attr(x, "type")))
  ret$digits <- digits
  ret$showprior <- showprior
  ret
}

#' @method print summary.shrinkTVP
#' @export
print.summary.shrinkTVP <- function(x, ...) {
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
#' @param x \code{mcmc.tvp} object
#' @param probs vector of boundaries for credible intervals to plot for each point in time, with values in [0,1].
#' The largest and smallest value form the outermost credible interval, the second smallest and second largest the second outermost and so forth.
#' The default value is \code{c(0.025, 0.25, 0.75, 0.975)}. Note that there have to be the same number of probs
#' < 0.5 as there are > 0.5.
#' @param shaded single logical value or a vector of logical values, indicating whether or not to shade the area between the pointwise
#' credible intervals. If a vector is given, the first value given is used to determine if the area between the outermost credible interval
#' is shaded, the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the
#' number of quantile pairs. The default value is \code{TRUE}.
#' @param quantlines single logical value or a vector of logical values, indicating whether or not to draw borders along the pointwise
#' credible intervals. If a vector is given, the first value given is used to determine whether the outermost credible interval
#' is marked by lines, the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the
#' number of credible intervals. The defualt value is \code{FALSE}.
#' @param shadecol single character string or a vector of character strings. Determines the color of the shaded areas that represent
#' the credible intervals. If a vector is given, the first color given is used for the outermost area,
#' the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the number
#' of shaded areas. Has no effect if \code{shaded = FALSE}. The default value is \code{"skyblue"}.
#' @param shadealpha real number between 0 and 1 or a vector of real numbers between 0 and 1.
#' Determines the level of transparency of the shaded areas that represent
#' the credible intervals. If a vector is given, the first value
#' given is used for the outermost area, the second for the second outermost and so forth. Recycled in the usual fashion
#' if the vector is shorter than the number of shaded areas. Has no effect if \code{shaded = FALSE}.
#' The default value is \code{0.5}.
#' @param quantlty either a single integer in [0,6] or one of the character strings \code{"blank",
#' "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"} or a vector containing these. Determines the line type of the borders
#' drawn around the shaded areas that represent the credible intervals. Note that if a vector is supplied the elements have to either
#' be all integers or all character strings. If a vector is given, the first value given is used for the outermost area, the second for
#' the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the number of shaded areas.
#' Has no effect if \code{quantlines = FALSE}. The default value is \code{2}.
#' @param quantcol single character string or a vector of character strings. Determines the color of the borders drawn around the shaded
#' areas that represent the credible intervals. If a vector is given, the first color given is used for borders of the outermost area,
#' the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the number
#' of shaded areas. Has no effect if \code{quantlines = FALSE}. The default value is \code{"red"}.
#' @param quantlwd single real, positive number or a vector of real, positive numbers. Determines the line width of the borders
#' drawn around the shaded areas that represent the credible intervals. If a vector is given, the first number given is used for
#' the borders of the outermost area, the second for the second outermost and so forth. Recycled in the usual fashion if the vector
#' is shorter than the number of shaded areas. Has no effect if \code{quantlines = FALSE}. The default value is \code{1}.
#' @param drawzero single logical value determining whether to draw a horizontal line at zero or not. The default value is \code{TRUE}.
#' @param zerolty single integer in [0,6] or one of the character strings \code{"blank",
#' "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"}. Determines the line type of the horizontal line at zero. Has no effect
#' if \code{drawzero = FALSE}. The default value is \code{2}.
#' @param zerolwd single real, positive number. Determines the line width of the horizontal line at zero. Has no effect
#' if \code{drawzero = FALSE}. The default value is \code{1}.
#' @param zerocol single character string. Determines the color of the horizontal line at zero. Has no effect if \code{drawzero = FALSE}.
#' The default value is \code{"grey"}.
#' @param ... further arguments to be passed to \code{plot}.
#' @return Called for its side effects and returns invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#' data <- sim$data
#'
#' res <- shrinkTVP(y ~ x1 + x2, data)
#' plot(res$beta$beta_x1)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @export
plot.mcmc.tvp <- function(x, probs = c(0.025, 0.25, 0.75, 0.975),
                          shaded = TRUE, quantlines = FALSE,
                          shadecol = "skyblue", shadealpha = 0.5,
                          quantlty = 2, quantcol = "black", quantlwd = 0.5,
                          drawzero = TRUE, zerolty = 2, zerolwd = 1, zerocol = "grey", ...){

  # Input checking
  if (length(probs) == 0 || is.vector(probs) == FALSE | is.list(probs) == TRUE ){
    stop("probs has to be a vector")
  }

  # Check if all probs are specified correctly
  if (any(probs > 1) | any(probs < 0) | any(sapply(probs, numeric_input_bad_zer))){
    stop("all elements of probs have to be numbers between >= 0 and <= 1")
  }

  if (length(probs[probs < 0.5]) != length(probs[probs > 0.5])){
    stop ("There has to be an equal amount of probs above and below 0.5")
  }

  if (length(quantlines) == 0 || any(sapply(quantlines, bool_input_bad))){
    stop("all elements of quantlines have to be logical values")
  }

  if (length(shaded) == 0 || any(sapply(shaded, bool_input_bad))){
    stop("all elements of shaded have to be logical values")
  }

  if (is.character(shadecol)){
    if (any(sapply(shadecol, char_input_bad))){
      stop("shadecol has to be a string or a positive number or a vector of strings or positive numbers")
    }
  } else if (is.numeric(shadecol)){
    if (any(sapply(shadecol, numeric_input_bad))){
      stop("shadecol has to be a string or a positive number or a vector of strings or positive numbers")
    }
  } else {
    stop("shadecol has to be a string or a positive number or a vector of strings or positive numbers")
  }


  if (length(shadealpha) == 0 || any(shadealpha > 1) | any(shadealpha < 0) | any(sapply(shadealpha, numeric_input_bad_zer))){
    stop("all elements of shadealpha have to be numbers between >= 0 and <= 1")
  }

  if (length(quantlty) == 0 || any(sapply(quantlty, lty_input_bad))) {
    stop("every element of quantlty has to be an integer from 0 to 6 or one of 'blank', 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', or 'twodash'")
  }

  if (is.character(quantcol)){
    if (any(sapply(quantcol, char_input_bad))){
      stop("quantcol has to be a string or a positive number or a vector of strings or positive numbers")
    }
  } else if (is.numeric(quantcol)){
    if (any(sapply(quantcol, numeric_input_bad))){
      stop("quantcol has to be a string or a positive number or a vector of strings or positive numbers")
    }
  } else {
    stop("quantcol has to be a string or a positive number or a vector of strings or positive numbers")
  }

  if (length(quantlwd) == 0 || any(sapply(quantlwd, numeric_input_bad))){
    stop("quantlwd has to be either a positive number or a string of positive numbers")
  }

  if (bool_input_bad(drawzero)){
    stop("drawzero has to be a logical")
  }

  if (lty_input_bad(zerolty)){
    stop("zerolty has to be an integer from 0 to 6 or one of 'blank', 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', or 'twodash'")
  }

  if (numeric_input_bad(zerolwd)){
    stop("zerolwd has to be a positive number")
  }

  if (is.character(zerocol)){
    if (char_input_bad(zerocol)){
      stop("zerocol has to be a string or a positive number")
    }
  } else if (is.numeric(zerocol)){
    if (numeric_input_bad(zerocol)){
      stop("zerocol has to be a string or a positive number")
    }
  } else {
    stop("zerocol has to be a string or a positive number")
  }


  # Sort by ascending size
  probs <- c(probs, 0.5)
  probs <- sort(probs)

  # Calculate all quantiles
  if (length(probs) == 1){
    quants <- matrix()
  }

  # Create x_vec for plotting
  x_vec <- attr(x, "index")

  startpoint <- ifelse(length(x_vec) == ncol(x), 1, 2)
  quants <- apply(x[, startpoint:ncol(x)], 2, quantile, probs)

  # Extract all user specified args and add standard ones that are missing
  args <- list(...)
  standard_args <- list(ylim = c(min(quants), max(quants)), type = "l", xlab = "", ylab = "")
  missing_args <- names(standard_args)[!names(standard_args) %in% names(args)]
  args[missing_args] <- standard_args[missing_args]

  # Add x and y to arguments
  args$x <- x_vec
  if (length(probs) == 1){
    args$y <- quants
  } else {
    args$y <- quants[median(1:nrow(quants)), ]
  }

  # Create plot
  do.call(plot, args)


  # Add horizontal line at 0
  if (drawzero){
    abline(h = 0, lty = zerolty, col = zerocol, lwd = zerolwd)
  }


  nint <- ((nrow(quants) - 1)/2)
  for (i in 1:nint){
    currcol <- ifelse(length(shadecol) > 1, rep_len(shadecol, nint)[i], shadecol)
    curralpha <- ifelse(length(shadealpha) > 1,  rep_len(shadealpha, nint)[i], shadealpha)

    currlty <-  ifelse(length(quantlty) > 1, rep_len(quantlty, nint)[i], quantlty)
    currbordcol <- ifelse(length(quantcol) > 1, rep_len(quantcol, nint)[i], quantcol)
    currlwd <-  ifelse(length(quantlwd) > 1, rep_len(quantlwd, nint)[i], quantlwd)
    currbord <- ifelse(length(quantlines) > 1, rep_len(quantlines, nint)[i], quantlines)

    polygon(c(x_vec, rev(x_vec)), c(quants[i, ], rev(quants[nrow(quants) + 1 - i, ])),
            col = adjustcolor(currcol, alpha.f = curralpha), border = FALSE)

    if (currbord == TRUE){
      lines(y = quants[i, ], x = x_vec, col = currbordcol, lwd = currlwd, lty = currlty)
      lines(y = quants[nrow(quants) + 1 - i, ], x = x_vec, col = currbordcol, lwd = currlwd, lty = currlty)
    }
  }

  # re paint line to be on top
  do.call(lines, args[names(args) %in% c("y", "x", "lwd", "lty", "col")])

}

#' Graphical summary of posterior distribution
#'
#' \code{plot.shrinkTVP} generates plots visualizing the posterior distribution.
#'
#' @param x a \code{shrinkTVP} object.
#' @param pars a character vector containing the names of the parameters to be visualized.
#' The names have to coincide with the names of the list elements of the \code{shrinkTVP}
#' object. Throws an error if any element of \code{pars} does not fulfill this criterium.
#' The default is \code{c("beta")}.
#' @param nplot positive integer that indicates the number of tvp plots to display on a single
#' page before a new page is generated. The default value is 3.
#' @param mar A numerical vector of the form \code{c(bottom, left, top, right)} which gives the number of lines of margin to be
#' specified on the four sides of the plot, as in \code{\link{par}}. The default is c(2, 4, 1, 2) + 0.1.
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
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @export
plot.shrinkTVP <- function(x, pars = c("beta"), nplot = 3, mar = c(2, 4, 1, 2) + 0.1, ...){


  # Check that any pars are specified
  if (length(pars) == 0){
    stop("Please specify pars to plot.")
  }

  # Check if pars is of type "char"
  if (is.character(pars) == FALSE){
    stop("pars has to be of type characater")
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


  # Save mfrow and mar value so that it can be restored after plotting
  prev_mfrow <- par()$mfrow
  prev_mar <- par()$mar


  # Create sublist containing only parameters specified in par and copy attributes
  obj <- x[pars]
  mostattributes(obj) <- attributes(x)
  names(obj) <- names(x[pars])

  # Idea: loop over all list elements and apply plot method specific to that object (mcmc or mcmc.tvp)
  for (i in 1:length(obj)){

    # TVP parameters will be a list, so we differentiate between TVP and non-TVP this way
    if (is.list(obj[[i]]) == TRUE){

      # Plot a max of nplot objects per plot
      if (length(obj[[i]]) > nplot){

        # Change mfrow and mar
        par(mfrow = c(nplot, 1))

        currpage <- 1
        p_left <- length(obj[[i]])

        # Apply plotting method in a loop
        for (j in 1:length(obj[[i]])){

          if ((j-1) %% nplot == 0){
            # The first plot on a new page
            currmar <- c(0, 1, 1, 1) * mar
            xaxt <- "n"
          } else if ((j-1) %% nplot == (nplot - 1)){
            # The last plot on a page
            currmar <- c(1, 1, 0, 1) * mar
            xaxt <- "s"
          } else if (p_left < nplot & j == length(obj[[i]])){
            # The last plot on the last page
            currmar <- c(1, 1, 0, 1) * mar
            xaxt <- "s"
          } else {
            # Any plot in between two others
            currmar <- c(0, 1, 0, 1) * mar
            xaxt <- "n"
          }
          par(mar = currmar)

          # Extract all user specified args
          args <- list(...)
          # Augment user specified arguments with standard ones
          standard_args <- list(ylab = paste(names(obj)[i], "of", attr(obj, "colnames")[j]), xaxt = xaxt)
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

          p_left <- p_left - currpage * nplot
          currpage <- currpage + 1
        }

        # Move to next series of plots if current parameter has been thoroughly plotted
        if (i < length(obj)){
          readline("Hit <Return> to see next plot: ")
        }

      } else {

        # If there are less than nplot parameters to plot, adjust mfrow to nr of parameters
        par(mfrow = c(nplot, 1))

        for (j in 1:length(obj[[i]])){

          # Adjust mar
          if (j == 1){
            currmar <- c(0, 1, 1, 1) * mar
            xaxt <- "n"
          } else if (j == length(obj[[i]])) {
            currmar <- c(1, 1, 0, 1) * mar
            xaxt <- "s"
          } else {
            currmar <- c(0, 1, 0, 1) * mar
            xaxt <- "n"
          }
          par(mar = currmar)


          # Extract all user specified args
          args <- list(...)
          # Augment user specified arguments with standard ones
          standard_args <- list(ylab = paste(names(obj)[i], "of", attr(obj, "colnames")[j]), xaxt = xaxt)
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
  par(mfrow = prev_mfrow, mar = prev_mar)
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

char_input_bad <- function(x){
  if (is.scalar(x) == TRUE){
    return(is.na(x) | is.character(x) == FALSE)
  } else {
    return(TRUE)
  }
}

lty_input_bad <- function(x){
  if (is.scalar(x) == TRUE){
    return((x %in% 0:6 | x %in% c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) == FALSE)
  } else {
    return(TRUE)
  }
}
