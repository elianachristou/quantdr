library(quantdr)

test_that("Local linear estimation", {
  set.seed(1234)
  n <- 100
  x <- rnorm(100); error <- rnorm(100); y <- x^2 + error
  tau <- c(0.25, 0.5, 0.75)
  par(mfrow=c(2,2))
  for (i in 1:length(tau)) {
    plot(x, y, main=print(tau[i]))
    points(x, llqr(x, y, tau=tau[i])$ll_est, pch=16, col='red')
  }
})

#test_that("Local linear estimation", {
#  set.seed(1234)
#  n <- 100
#  x <- rnorm(100); error <- rnorm(100); y <- x^2 + error
#  tau <- 0.5
#  expect_equal(llqr(x, y, tau = tau, method = 'CV')$h, 0.3981072)
#})

test_that("the function gives an error message when the dimension of
          x0 is less than that of x", {
  set.seed(1234)
  n <- 100
  x <- cbind(rnorm(n), rnorm(n)); error <- rnorm(100); y <- x[, 1]^2 + error
  tau <- 0.5
  x0 <- 1
  expect_error(llqr(x, y, tau=tau, x0=x0)$ll_est, )
})
