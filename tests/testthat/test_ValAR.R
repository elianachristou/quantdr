library(quantdr)

test_that("the function returns an error when y is multivariate response", {
  data(edhec, package = "PerformanceAnalytics")
  y <- cbind(as.vector(edhec[, 1]), as.vector(edhec[, 2]))
  tau <- 0.05
  p <- 5
  expect_error(ValAR(y, p, tau = tau), )
})

test_that("the function returns an error when n is less than p", {
  data(edhec, package = "PerformanceAnalytics")
  y <- as.vector(edhec[1:4, 2])
  tau <- 0.05
  p <- 5
  expect_error(ValAR(y, p, tau = tau), )
})

test_that("the function returns an error when tau is more than one-dimension", {
  data(edhec, package = "PerformanceAnalytics")
  y <- as.vector(edhec[, 2])
  tau <- c(0.01, 0.05)
  p <- 5
  expect_error(ValAR(y, p, tau = tau), )
})

test_that("the function returns an error when tau is outside the (0, 1) interval", {
  data(edhec, package = "PerformanceAnalytics")
  y <- as.vector(edhec[, 2])
  tau <- 1.1
  p <- 5
  expect_error(ValAR(y, p, tau = tau), )
})

test_that("the function returns an error when the moving window is not an integer", {
  data(edhec, package = "PerformanceAnalytics")
  y <- as.vector(edhec[, 2])
  tau <- 0.05
  p <- 5
  movwind <- 150.5
  expect_error(ValAR(y, p, tau = tau, movwind = movwind), )
})

test_that("the function returns an error when the moving window exceeds the number of observations", {
  data(edhec, package = "PerformanceAnalytics")
  y <- as.vector(edhec[, 2])
  tau <- 0.05
  p <- 5
  movwind <- 400
  expect_error(ValAR(y, p, tau = tau, movwind = movwind), )
})

