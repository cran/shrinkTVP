context("test shrinkTVP")
library(shrinkTVP)



test_that("shrinkTVP executes", {

  data <- data.frame(y = runif(2), x1 = rnorm(2), x2 = rnorm(2))
  niter <- 5

  expect_s3_class(shrinkTVP(y ~ ., data, niter = niter, display_progress = FALSE), "shrinkTVP_res")

})

