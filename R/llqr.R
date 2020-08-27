#' Local linear conditional quantile estimation
#'
#' \code{llqr} estimates the \code{tau}-th conditional quantile of \code{y}
#' given \code{x} based on a local linear estimation. The estimation is
#' performed either at a point x0, if specified, or at the entire design matrix
#' x.
#'
#' The local linear estimation of the \code{tau}-th conditional quantile
#' estimator requires a kernel and a bandwidth.  This function uses the normal
#' kernel and it allows the bandwidth to be given by the user, or estimated
#' using either the rule-of-thumb bandwidth of Yu and Jones (1994) or the
#' Cross-Validation criterion. \eqn{a+b}
#'
#' @param x A design matrix.The rows represent observations and the columns
#'   represent predictor variables.
#' @param y A vector of the respons variable.
#' @param h A univariate bandwidth.
#' @param tau A quantile level, a number strictly between 0 and 1.
#' @param method A method to select the bandwidth, if it is missing. There are
#'   two methods considered: the rule-of-thumb bandwidth of Yu and Jones (1994)
#'   and the cross-validation criterion.
#' @param x0 A single observation for which to perform the estimation. It needs
#'   to be a singular value, for a univariate predictor, or a vector, for a
#'   multivariate predictor. If \code{x0} is missing, the estimation will be
#'   performed on the entire design matrix \code{x}.
#' @return \code{llqr} returns the local linear \code{tau}-th conditional
#'   quantile estimator of \code{y} given \code{x}, along with the bandwidth
#'   used.
#' @include llqrcv.R
#' @import stats
#' @examples
#' n <- 100
#' x <- rnorm(100); error <- rnorm(100); y <- x^2 + error
#' tau <- c(0.25, 0.5, 0.75)
#' par(mfrow=c(2,2))
#' for (i in 1:length(tau)) {
#'   plot(x, y, main=print(tau[i]))
#'   points(x, llqr(x, y, tau=tau[i])$ll_est, pch=16, col='red')
#' }
#'
#' \dontrun{
#' n <- 100
#' x <- cbind(rnorm(n), rnorm(n)); error <- rnorm(100); y <- x[, 1]^2 + error
#' x0=1
#' tau=0.5
#' llqr(x, y, tau=tau, x0=x0)$ll_est
#' }
#' @export
llqr <- function(x, y, tau=0.5, h, method="rule", x0) {

  x <- as.matrix(x)
  n <- length(y)
  p <- dim(x)[2]

  if (missing(h)) {
    if (method == "CV") {
      h <- llqrcv(x, y, tau)
    }
    if (method == "rule") {
      h <- KernSmooth::dpill(x, y);
      h <- h * (tau * (1 - tau) / (dnorm(qnorm(tau)))^2)^.2
    }
  } else {
    h <- h
  }

  if (missing(x0)) {
    ll_est <- as.null(n)

    if (p == 1) {
      for (i in 1:dim(x)[1]) {
        z <- x - x[i, 1]
        w <- dnorm(z / h)
        q <- quantreg::rq(y ~ z, weights = w, tau = tau, ci = FALSE)
        ll_est[i] <- q$coef[1]
      }
    } else {
      for (i in 1:dim(x)[1]) {
        z <- matrix(0, n, p)
        z <- x - t(matrix(rep(x[i, ], p * n), p, n))
        w <- mvtnorm::dmvnorm(z / h)
        q <- quantreg::rq(y ~ z, weights = w, tau = tau, ci = FALSE)
        ll_est[i] <- q$coef[1]
      }
    }
  } else {
    if (p == 1) {
      z <- x - x0
      w <- dnorm(z / h)
      q <- quantreg::rq(y ~ z, weights = w, tau = tau, ci = FALSE)
      ll_est <- q$coef[1]
    } else {
      if (length(x0) != p) {
        print("The dimension of x0 needs to be the same as the number of columns of x")
        stop()
      }
      z <- matrix(0, n, p)
      z <- x - t(matrix(rep(x0, p * n), p, n))
      w <- mvtnorm::dmvnorm(z / h)
      q <- quantreg::rq(y ~ z, weights = w, tau = tau, ci = FALSE)
      ll_est <- q$coef[1]
    }
  }
  list(ll_est = ll_est, h = h)
}

#' @describeIn llqr Cross-Validation for bandwidth selection for local linear
#'   conditional quantile estimation
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
