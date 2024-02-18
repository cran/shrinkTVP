context("test shrinkDTVP")

test_bed <- function(args, dummies = FALSE, transf = FALSE) {

  set.seed(123)
  full_dat <- data.frame(y = runif(16), x1 = rnorm(16), x2 = runif(16), x3 = factor(rep(c(1, 2), 8)))
  test <- full_dat[11:16, ]
  args$data <- full_dat[1:10, ]

  args$formula <- y ~ x1 + x2 + 1

  args$a_psi <- 3
  args$c_psi <- 5

  if (dummies == TRUE) {
    args$formula <- update.formula(args$formula, "~ . + x3")
  }

  if (transf == TRUE) {
    args$formula <- update.formula(args$formula, "~ . + I(x1^2) + I(log(x2))")
  }

  res <- do.call(shrinkDTVP, args)

  expect_s3_class(res, "shrinkTVP")

  # Test methods
  expect_invisible(plot(res, nplot = 7))
  expect_visible(summary(res))
  expect_invisible(plot(res$beta[[1]]))
  expect_visible(res)

  # Test prediction related functions/methods
  expect_type(LPDS(res, test[1,]), "double")
  expect_type(eval_pred_dens(-5:5, res, test[1,]), "double")
  expect_length(eval_pred_dens(-5:5, res, test[1,]), length(-5:5))

  expect_s3_class(forecast_shrinkTVP(res, test), "shrinkTVP_forc")
  expect_s3_class(predict(res), "shrinkTVP_pred")
  expect_s3_class(fitted(res), "shrinkTVP_fitted")
  expect_s3_class(residuals(res), "shrinkTVP_resid")
}

mod_type = c("triple", "double", "ridge")
transformations <- c(TRUE, FALSE)
dummies <- c(TRUE, FALSE)
#p <- c(0, 1, 2)
scenarios <- expand.grid(mod_type, transformations, dummies)#, p)
names(scenarios) <- c("mod_type", "transf", "dummies")#, "p")

params <- c(
  "learn_a_xi",
  "learn_a_tau",
  "learn_c_xi",
  "learn_c_tau",
  "a_eq_c_xi",
  "a_eq_c_tau",
  "learn_kappa2_B",
  "learn_lambda2_B",
  "sv",
  "a_xi_adaptive",
  "c_xi_adaptive",
  "a_tau_adaptive",
  "c_tau_adaptive",
  "iid",
  "shrink_inter"
)

for(i in length(scenarios)) {

  for (j in params) {

    args <- formals(shrinkDTVP)
    args <- args[sapply(args, function(x) x != "")]

    if (grepl("adaptive", j)) {
      args$MH_tuning$temp <- FALSE
      names(args$MH_tuning) <- j
    } else {
      args[[j]] <- !args[[j]]
    }

    #args$p <- scenarios$p[i]
    args$mod_type <- as.character(scenarios$mod_type[i])
    args$display_progress <- FALSE
    args$niter <- 10
    args$nburn <- 0

    test_that(paste0("scenario: ", i,  ", mod_type: ", scenarios$mod_type[i],
                     ", transformations: ", scenarios$transf[i],
                     ", dummies: ", scenarios$dummies[i],
                     #", p: ", scenarios$p[i],
                     ", toggled: ", j), {
                       test_bed(args, scenarios$dummies[i], scenarios$transf[i])
                     })

  }
}
