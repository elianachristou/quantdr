#' Cross-Validation for bandwidth selection for local linear conditional
#' quantile estimation
#'
#' @param x a design matrix.The rows represent observations and the columns
#'   represent predictor variables
#' @param y a vector of the respons variable
#' @param tau the quantile level, a number strictly between 0 and 1.
#' @return \code{llqrcv} returns the optimal bandwidth selected using
#'   Cross-Validation criterion for the estimation of the \code{tau}-th
#'   conditional quantile of \code{y} given \code{x}.
#' @include llqr.R
#' @examples
#' n <- 100
#' x <- rnorm(100); error <- rnorm(100); y <- x^2 + error
#' tau <- 0.5
#' llqrcv(x, y, tau = tau)
llqrcv <- function(x, y, tau=0.5) {
  n <- length(y)
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
