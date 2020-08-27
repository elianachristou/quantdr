#' Local linear quantile estimation
#'
#' \code{llqr} estimates the \eqn{\tau}-th conditional quantile of \code{y}
#' given \code{x} based on a local linear fitting.  The estimation is performed
#' at each of the design points or, if specified, at a single observation point.
#'
#' The function computes the local linear quantile regression fit for a specified
#' quantile \eqn{\tau} at the design points of the matrix \code{x} or at a
#' pre-specified point \code{x0}.  The estimation is based on a standard normal
#' kernel and a univariate bandwidth, which, if not specified by the user, it is
#' defined using either the rule-of-thumb given by Yu and Jones (1994) or the
#' cross-validation criterion.
#'
#' The estimation applies to univariate and multivariate predictor variables.
#' For the later, the function uses the multivariate standard normal kernel.
#' Note that, if the estimation is performed at a pre-specified point \code{x0},
#' then \code{x0} should be a scaler (for univariate predictor) or a vector (for
#' multivariate predictor).
#'
#' @param x A design matrix.  The rows represent observations and the columns
#'   represent predictor variables.
#' @param y A vector of the response variable.
#' @param h A univariate bandwidth.
#' @param tau A quantile level, a number strictly between 0 and 1.
#' @param method A method to select the bandwidth, if it is missing.  Use
#' \code{rule} for the rule-of-thumb bandwidth of Yu and Jones (1994) or
#' \code{CV} for the method of cross-validation.
#' @param x0 A single observation for which to perform the estimation.  It needs
#'   to be a singular value, for a univariate predictor, or a vector, for a
#'   multivariate predictor.  If \code{x0} is missing, the estimation will be
#'   performed on the design matrix \code{x}.
#' @return \code{llqr} computes the local linear \eqn{\tau}-th conditional
#'   quantile estimator of \code{y} given \code{x}, and returns:
#'
#'   \code{ll_est}   The estimated function value at the design points \code{x}
#'   or, if specified, at the point \code{x0}.
#'
#'   \code{h}  The bandwidth.
#'
#' @references Yu, K. and Jones, M.C. (1998), Local linear quantile regression.
#' Journal of the American Statistical Association, 93, 228-237.
#' @import stats
#' @examples
#'
#' set.seed(1234)
#' n <- 100
#' x <- rnorm(100); error <- rnorm(100); y <- (x + 1)^3 + 0.1 * (x - 2)^3 + error
#' tau <- 0.5
#' plot(x, y)
#' points(x, llqr(x, y, tau=tau)$ll_est, col='red', pch=16)
#' x0=1
#' llqr(x, y, tau=tau, x0=x0)
#'
#' n <- 100; p <- 2
#' x <- matrix(rnorm(p * n), n, p); error <- rnorm(100); y <- exp(x[,1]+x[,2]) + error
#' tau <- 0.5
#' x0=c(1,2)
#' llqr(x, y, tau=tau, x0=x0)
#'
#' require(MASS)
#' data(mcycle)
#' attach(mcycle)
#' plot(times,accel,xlab = "milliseconds", ylab = "acceleration (in g)")
#' taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' for(i in 1:length(taus)) {
#'  fit <- llqr(times, accel, tau = taus[i])$ll_est
#'  lines(times, fit ,lty=i)
#' }
#' legend(38,-50,c("tau=0.1","tau=0.25","tau=0.5","tau=0.75", "tau=0.9"), lty=1:length(taus))
#'
#' plot(times,accel,xlab = "milliseconds", ylab = "acceleration (in g)")
#' tau <- 0.5
#'  lines(times, llqr(times, accel, tau = taus[i], h=0.5)$ll_est ,lty=1)
#'  lines(times, llqr(times, accel, tau = taus[i], method="rule")$ll_est ,lty=2)
#'  lines(times, llqr(times, accel, tau = taus[i], method="CV")$ll_est ,lty=3)
#' legend(40,-70,c("h=0.5","h=rule","h=CV"), lty=1:3)
#'
#' set.seed(1234)
#' n <- 100
#' x <- rnorm(100); error <- rnorm(100); y <- (x + 1)^3 + 0.1 * (x - 2)^3 + error
#' tau <- c(0.25, 0.5, 0.75)
#' plot(x, y)
#' for (i in 1:length(tau)) {
#'   points(x, llqr(x, y, tau=tau[i])$ll_est, col = i)
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
