# Wrapper around C++ function that deals with special case where beta is a matrix in the univariate case
calc_fitted <- function(y, x, beta){
  if (is.list(beta) == FALSE){
    beta <- list(beta = beta)
  }
  calc_fitted_cpp(y, x, beta)
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

# A function to produce lags
mlag <- function(X, p){
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

# A convenience function for printing the posterior distribution in the summary method
dist_parser <- function(char_rep, summ_obj){

  dict <- list(beta_mean = "\u03b2\u2c7c",
               tau2 = "\u03c4\u00b2\u2c7c",
               lambda2 = "\u03bb\u00b2\u2c7c",
               a_tau = "a^\u03c4",
               alpha_a_tau = "\u03B1_a\u03c4",
               beta_a_tau = "\u03b2_a\u03c4",
               c_tau = "c^\u03c4",
               alpha_c_tau = "\u03B1_c\u03c4",
               beta_c_tau = "\u03b2_c\u03c4",
               lambda2_B = "\u03bb\u00b2_B",
               theta_sr = "\u221a\u03b8\u2c7c",
               xi2 = "\u03be\u00b2\u2c7c",
               kappa2 = "\u03ba\u00b2\u2c7c",
               a_xi = "a^\u03be",
               alpha_a_xi = "\u03B1_a\u03be",
               beta_a_xi = "\u03b2_a\u03be",
               c_xi = "c^\u03be",
               alpha_c_xi = "\u03B1_c\u03be",
               beta_c_xi = "\u03b2_c\u03be",
               kappa2_B = "\u03ba\u00b2_B",
               sigma2 = "\u03c3\u00b2",
               C0 = "C\u2080",
               sv_mu = "\u03bc\u209B\u1D65",
               sv_phi = "\u03c6\u209B\u1D65",
               sv_sigma2 = "\u03c3\u00b2\u209B\u1D65",

               # New dynamic parts
               w2 = "w\u00b2\u2c7c\u209c",
               a_psi = "a\u2c7c",
               c_psi = "c\u2c7c",
               rho_p = "\u03c1\u2c7c",
               theta = "\u03b8\u2c7c",
               iid = "iid ",
               a_rho_sym = "a\u1D68",
               b_rho_sym = "b\u1D68",
               alpha_rho_sym = "\u03B1\u1D68",
               beta_rho_sym = "\u03B2\u1D68",

               # Operators
               `,` = ", ",
               `+` = " + ",
               `~` = " ~ ")

  split <- unlist(strsplit(char_rep, split = " "))

  # Have to hard code in w2, as it is technically not "sampled" (only indirectly)
  if ((split[min(which(split %in% names(dict)))] %in% c(summ_obj$names, "w2"))) {

    # Select all parameters that potentially appear in priors and that were not sampled
    select_static_param <- split %in% names(summ_obj$priorvals) & !split %in% summ_obj$names

    # Preserve all parameters in lhs of the equation (these should stay symbolic)
    if ("|" %in% split & "~" %in% split){
      lhs_params <- rep(TRUE, length(split))
      lhs_params[1:which(split == "~")] <- FALSE
      select_static_param <- select_static_param & lhs_params
    }

    # Everything else should be replaced by its dict representation
    to_replace <- split[split %in% names(dict)]

    # Replace non-sampled parameters with prior values
    for (i in seq_along(split)){
      if (select_static_param[i] == TRUE){
        split[i] <- round(unlist(summ_obj$priorvals[names(summ_obj$priorvals) == split[i]]), summ_obj$digits)
      }
    }

    # Replace everything else with formula representation
    for (i in to_replace){
      split[split == i] <- dict[[i]]
    }

    # Add a line break at the end
    split <- c(split, "\n")

    # And print!
    cat(paste(split, collapse = ""))
  }
}

# Merges user and default values of named list inputs
list_merger <- function(default, user) {

  # Check that user and sv_param are a list
  if (is.list(user) == FALSE | is.data.frame(user)){
    stop(paste0(deparse(substitute(user)), " has to be a list"))
  }

  stand_nam <- names(default)
  user_nam <- names(user)

  # Give out warning if an element of the parameter list is misnamed
  if (any(!user_nam %in% stand_nam)){
    wrong_nam <- user_nam[!user_nam %in% stand_nam]
    warning(paste0(paste(wrong_nam, collapse = ", "),
                   ifelse(length(wrong_nam) == 1, " has", " have"),
                   " been incorrectly named in ", deparse(substitute(user)), " and will be ignored"),
            immediate. = TRUE)
  }

  # Merge users' and default values and ignore all misnamed values
  missing_param <- stand_nam[!stand_nam %in% user_nam]
  user[missing_param] <- default[missing_param]
  user <- user[stand_nam]

  return(user)
}
