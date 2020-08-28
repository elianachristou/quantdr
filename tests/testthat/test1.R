library(quantdr)

test_that("the function gives an error message when the dimension of
          x0 is less than that of x", {
  set.seed(1234)
  n <- 100; p <- 2
  x <- matrix(rnorm(n * p), n, p); error <- rnorm(n); y <- x[, 1]^2 + error
  tau <- 0.5
  x0 <- 1
  expect_error(llqr(x, y, tau = tau, x0 = x0)$ll_est, )
})
