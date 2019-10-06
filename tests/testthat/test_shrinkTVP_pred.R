context("test shrinkTVP prediction functions")
library(shrinkTVP)



test_that("shrinkTVP LPDS and eval_pred_dens function", {

  set.seed(123)
  data <- data.frame(y = runif(5), x1 = rnorm(5), x2 = rnorm(5))
  niter <- 5

  res <- shrinkTVP(y ~ ., data[1:4,], niter = niter, display_progress = FALSE)

  expect_type(LPDS(res, data[5,]), "double")
  expect_type(eval_pred_dens(-5:5, res, data[5,]), "double")
  expect_length(eval_pred_dens(-5:5, res, data[5,]), length(-5:5))
})

