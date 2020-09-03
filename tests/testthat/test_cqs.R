library(quantdr)

test_that("the function returns an error when y is multivariate response", {
  set.seed(1234)
  n <- 100; p <- 2
  x <- matrix(rnorm(n * p), n, p); error <- rnorm(n); y <- x^2 + error
  tau <- 0.5
  expect_error(cqs(x, y, tau = tau), )
})

test_that("the function returns an error when the number of observations
          for y and x differ", {
            set.seed(1234)
            n <- 100; p <- 2
            x <- matrix(rnorm(n * p), n, p); error <- rnorm(n); y <- x[, 1]^2 + error
            tau <- 0.5
            expect_error(cqs(x, y[1:(n / 2)], tau = tau), )
          })

test_that("the quantile level needs to be strictly between 0 and 1", {
  set.seed(1234)
  n <- 100
  x <- rnorm(n); error <- rnorm(n); y <- (x + 1)^3 + 0.1 * (x - 2)^3 + error
  tau <- -0.1
  expect_error(cqs(x, y, tau = tau), )
})

test_that("the function returns an error for NA", {
  set.seed(1234)
  n <- 100
  x <- rnorm(n); error <- rnorm(n); y <- (x + 1)^3 + 0.1 * (x - 2)^3 + error
  y[1] <- 0 / 0
  tau <- 0.5
  expect_error(cqs(x, y, tau = tau), )
})

test_that("quantile level needs to be one number", {
  set.seed(1234)
  n <- 100
  x <- rnorm(n); error <- rnorm(n); y <- (x + 1)^3 + 0.1 * (x - 2)^3 + error
  tau <- c(0.2, 0.5)
  expect_error(cqs(x, y, tau = tau), )
})

test_that("n should be greater than p", {
  set.seed(1234)
  n <- 20; p <- 40
  x <- matrix(rnorm(n * p), n, p); error <- rnorm(n); y <- x[, 1]^2 + error
  tau <- 0.5
  expect_error(cqs(x, y, tau = tau), )
})

test_that("dtau needs to be an integer between 1 and p",{
  set.seed(1234)
  n <- 100; p <- 10
  x <- matrix(rnorm(n * p), n, p); error <- rnorm(n)
  y <- 3 * x[, 1] + x[, 2] + error
  tau <- 0.5
  expect_error(cqs(x, y, tau = tau, dtau = -1), )
  expect_error(cqs(x, y, tau = tau, dtau = -1.5), )
  expect_equal(cqs(x, y, tau = tau, dtau = 1)$dtau, 1)
  expect_error(cqs(x, y, tau = tau, dtau = 2.5), )
  expect_error(cqs(x, y, tau = tau, dtau = 20), )
})

test_that("x should be more than one-dimensional", {
  set.seed(1234)
  n <- 100; p <- 1
  x <- matrix(rnorm(n * p), n, p); error <- rnorm(n); y <- x[, 1]^2 + error
  tau <- 0.5
  expect_error(cqs(x, y, tau = tau), )
})

