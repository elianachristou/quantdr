library(quantdr)

test_that("the BIC criterion correctly estimates the dimension", {
  set.seed(1234)
  n <- 100; p <- 10; tau <- 0.5
  x <- matrix(rnorm(n * p), n, p); error <- rnorm(n)
  y <- 3 * x[, 1] + x[, 2] + error
  out <- cqs(x, y, tau)
  expect_equal(bic_d(out$qvalues, n), 1)
})


test_that("the BIC criterion correctly estimates the dimension", {
  lambdas <- c(1, 2, 0.05, 0.06)
  n <- 100
  expect_equal(bic_d(lambdas, n), 2)
})
