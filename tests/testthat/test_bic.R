library(quantdr)

test_that("the BIC criterion correctly estimates the dimension", {
  set.seed(1234)
  n <- 100; p <- 10; tau <- 0.5
  x <- matrix(rnorm(n * p), n, p); error <- rnorm(n)
  y <- 3 * x[, 1] + x[, 2] + error
  out <- cqs(x, y, tau)
  expect_equal(bic_d(out$qvalues, n), 1)
})
