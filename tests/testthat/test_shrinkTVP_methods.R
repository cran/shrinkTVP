context("test shrinkTVP methods")
library(shrinkTVP)



test_that("shrinkTVP methods function", {

  set.seed(123)
  data <- data.frame(y = runif(2), x1 = rnorm(2), x2 = rnorm(2))
  niter <- 5

  res <- shrinkTVP(y ~ ., data, niter = niter, display_progress = FALSE)


  expect_invisible(plot(res))
  expect_visible(summary(res))
  expect_invisible(plot(res$beta[[1]]))
  expect_visible(res)
})

