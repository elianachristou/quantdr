#' Cross-Validation for bandwidth selection of local linear quantile regression
#'
#' \code{llqrcv} estimates the bandwidth necessary for the local linear fit of
#' the \eqn{\tau}th conditional quantile of \code{y} given \code{x}.  The
#' estimation is performed using the Cross-Validation criterion.
#'
#' A grid of bandwidth values is created and the local linear fit is estimated
#' using all the data points except for one point, which is used to make the
#' prediction.  This procedure is repeated \code{n} times, where \code{n} is the
#' number of observations.  Then, the bandwidth is selected as the one with the
#' smallest average error.
#'
#' @param x A design matrix.  The rows represent observations and the columns
#'   represent predictor variables.
#' @param y A vector of the response variable.
#' @param tau A quantile level, a number strictly between 0 and 1.
#' @return \code{llqrcv} returns the optimal bandwidth selected using
#'   Cross-Validation criterion for the local linear fit of the \eqn{\tau}th
#'   conditional quantile of \code{y} given \code{x}.
#' @include llqr.R
#' @examples
#' set.seed(1234)
#' n <- 100
#' x <- rnorm(100); error <- rnorm(100); y <- x^2 + error
#' tau <- 0.5
#' llqrcv(x, y, tau = tau)
#' @export
llqrcv <- function(x, y, tau=0.5) {
  n <- length(y)
  # create grid values for bandwidth
  h_lower <- min(n^ (- 1 / 5), min(2, sd(y)) * n^ (- 1 / 5))
  h_upper <- max(n^ (- 1 / 5), min(2, sd(y)) * n^ (- 1 / 5))
  h <- seq(h_lower, h_upper, by = 0.1)
  x <- as.matrix(x)
  qhat <- matrix(0, n, length(h))
  cv <- as.null(length(h))

  for (i in 1:length(h)) {
    for (j in 1:n) {
      qhat[j, i] <- llqr(x = x[-j, ], y = y[-j], tau = tau, h = h[i],
                         x0 = x[j, ])$ll_est
    }
    cv[i] <- mean((y - qhat[, i]) * (tau - (y < qhat[, i])))
  }
  h_opt <- h[which(cv == min(cv))]
  return(h_opt)
}
