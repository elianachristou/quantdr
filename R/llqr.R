#' Local linear quantile regression
#'
#' \code{llqr} estimates the \eqn{\tau}th conditional quantile of \code{y} given
#' \code{x} based on a local linear fit.  The estimation is performed at each of
#' the design points or, if specified, at a single observation point \code{x0}.
#'
#' The function computes the local linear quantile regression fit for a specified
#' quantile level \eqn{\tau} at the design points of the matrix \code{x} or at a
#' pre-specified point \code{x0}.  The estimation is based on a standard normal
#' kernel and a univariate bandwidth.  The bandwidth, if not specified by the
#' user, is defined using either the rule-of-thumb given by Yu and Jones (1994)
#' or the cross-validation criterion.
#'
#' The estimation applies to univariate and multivariate predictor variables. For
#' the latter, the local linear fit uses the multivariate standard normal kernel.
#' Note that if the estimation is performed at a pre-specified point \code{x0},
#' then \code{x0} should be a scalar (for univariate predictor) or a vector (for
#' multivariate predictor).
#'
#' @section Warning:
#' The user needs to be careful about the bandwidth selection. When the
#' dimension of the predictor variable is large compared to the sample size,
#' local linear fitting meets the 'curse of dimensionality' problem. In
#' situations like that, the bandwidth selected by the rule-of-thumb or the
#' cross- validation criterion might be small and lead to a sparse
#' neighborhood. This will cause the function to fail. For these cases, we
#' advice the user to specify a bandwidth in the function. See the last example
# below.
#'
#' @param x A design matrix (n x p).  The rows represent observations and the
#' columns represent predictor variables.
#' @param y A vector of the response variable.
#' @param tau A quantile level, a number strictly between 0 and 1.
#' @param h A univariate bandwidth.  If not specified, the bandwidth is estimated
#'    using either "\code{rule}" or "\code{CV}".  See \code{method} below for
#'    details.
#' @param method A character string specifying the method to select the
#'    bandwidth, if it is missing.  Use "\code{rule}" for the rule-of-thumb
#'    bandwidth of Yu and Jones (1994) or "\code{CV}" for the method of
#'    cross-validation.
#' @param x0 A single observation for which to perform the estimation.  It needs
#'    to be a singular value, for a univariate predictor, or a vector, for a
#'    multivariate predictor.  If \code{x0} is missing, the estimation will be
#'    performed on the design matrix \code{x}.
#' @return \code{llqr} computes the local linear \eqn{\tau}th conditional
#'    quantile function of \code{y} given \code{x}, and returns:
#'    \itemize{
#'    \item{ll_est: }{The estimated function value at the design points \code{x}
#'    or, if specified, at the point \code{x0}.}
#'
#'  \item{h: }{The bandwidth for the local linear quantile regression fit.  If
#'  not specified by the user, \code{h} is estimated using either the
#'  rule-of-thumb given by Yu and Jones (1994) or the cross-validation
#'  criterion.} }
#'
#' @references Yu, K., and Jones, M.C. (1998), Local linear quantile regression.
#'  \emph{Journal of the American Statistical Association}, 93, 228-237.
#' @import stats
#' @include llqrcv.R
#' @examples
#'
#' # Example 1
#' # estimate the function for different quantile levels for simulated data
#' set.seed(1234)
#' n <- 100
#' x <- rnorm(n)
#' error <- rnorm(n)
#' y <- (x + 1)^3 + 0.1 * (x - 2)^3 + error
#' taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' par(mfrow=c(2, 3))
#' for (i in 1:length(taus)) {
#' plot(x, y, main = taus[i])
#' points(x, llqr(x, y, tau = taus[i])$ll_est, col = 'red', pch = 16)
#' }
#'
#' # Example 2
#' # estimate the function at a point x0
#' set.seed(1234)
#' n <- 100
#' x <- rnorm(n)
#' error <- rnorm(n)
#' y <- (x + 1)^3 + 0.1 * (x - 2)^3 + error
#' tau <- 0.5
#' x0 <- 1
#' llqr(x, y, tau = tau, x0 = x0)
#'
#' # Example 3
#' # estimate the function for different quantile levels
#' par(mfrow = c(1, 1))
#' data(mcycle, package = "MASS")
#' attach(mcycle)
#' plot(times, accel, xlab = "milliseconds", ylab = "acceleration")
#' taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' for(i in 1:length(taus)) {
#'  fit <- llqr(times, accel, tau = taus[i])$ll_est
#'  lines(times, fit, lty = i)
#' }
#' legend(45, -50, c("tau=0.1","tau=0.25","tau=0.5","tau=0.75", "tau=0.9"),
#'     lty=1:length(taus))
#'
#' # Example 4
#' # demonstrate a situation where the dimension of the predictor is large and
#' # the local linear fitting meets the 'curse of dimensionality' problem
#' set.seed(1234)
#' n <- 100
#' p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' error <- rnorm(n)
#' y <- 3 * x[, 1] + x[, 2] + error
#' tau <- 0.5
#'\dontrun{
#' fit1 <- llqr(x, y, tau = tau)
#' fit2 <- llqr(x, y, tau = tau, method = "CV")
#' }
#' fit.alt <- llqr(x, y, tau = tau, h=1)
#' fit.alt
#'
#' @export
llqr <- function(x, y, tau=0.5, h = NULL, method="rule", x0 = NULL) {

  x <- as.matrix(x)
  y <- as.matrix(y)

  # compatibility checks
  # checks if y is univariate
  if (dim(y)[2] > 1) {
      stop(paste("y needs to be a univariate response. y is a", dim(y)[2], "-dimensional response in this case."))
  }
  # checks if the number of observations for x and y agree
  if (length(y) != dim(x)[1]) {
    stop(paste("number of observations of y (", length(y), ") not equal to the number of rows of x (", dim(x)[1], ").", sep = ""))
  }
  # checks if the quantile level is one-dimensional
  if (length(tau) > 1) {
    stop(paste("quantile level needs to be one number."))
  }
  # checks if the quantile level is between 0 and 1 (strictly)
  if (tau >= 1 | tau <= 0) {
    stop(paste("quantile level needs to be a number strictly between 0 and 1."))
  }
  # checks for NAs
  if (sum(is.na(y)) > 0 | sum(is.na(x)) > 0) {
    stop(paste("Data include missing values."))
  }
  # checks if n>p
  if (length(y) <= dim(x)[2]) {
    stop(paste("number of observations of y (", length(y), ") should be greater than the number of columns of x (", dim(x)[2], ").", sep = ""))
  }

  n <- length(y)
  p <- dim(x)[2]

  # if bandwidth is missing, estimate it
  if (is.null(h)) {
    if (method == "CV") {
      h <- llqrcv(x, y, tau)
    }
    if (method == "rule") {
      h <- KernSmooth::dpill(x, y)
      h <- h * (tau * (1 - tau) / (dnorm(qnorm(tau)))^2)^.2
    }
  } else {
    h <- h
  }

  # if x0 is missing, perform estimation at the entire design matrix
  if (is.null(x0)) {
    ll_est <- list()
    # if the dimension of x is one, use univariate kernel, otherwise
    # use multivariate kernel
    if (p == 1) {
      # perform estimation at the design matrix x
      for (i in 1:dim(x)[1]) {
        z <- x - x[i, 1]
        w <- dnorm(z / h)
        q <- quantreg::rq(y ~ z, weights = w, tau = tau, ci = FALSE)
        ll_est[i] <- q$coef[1]
      }
    } else {
      for (i in 1:dim(x)[1]) {
        z <- list()
        z <- t(t(x) - x[i, ])
        w <- mvtnorm::dmvnorm(z / h)
        q <- quantreg::rq(y ~ z, weights = w, tau = tau, ci = FALSE)
        ll_est[i] <- q$coef[1]
      }
    }
  } else {
    # checks if the dimension of x0 is the same as p
    if (length(x0) != p) {
      stop(paste("x0 needs to be a p-dimensional vector, where p is the dimension of x (",dim(x)[2], ")."))
    }
    if (p == 1) {
      # perform estimation at the point x0
      z <- x - x0
      w <- dnorm(z / h)
      q <- quantreg::rq(y ~ z, weights = w, tau = tau, ci = FALSE)
      ll_est <- q$coef[1]
    } else {
      z <- list()
      z <- t(t(x) - x[i, ])
      w <- mvtnorm::dmvnorm(z / h)
      q <- quantreg::rq(y ~ z, weights = w, tau = tau, ci = FALSE)
      ll_est <- q$coef[1]
    }
  }
  list(ll_est = unlist(ll_est), h = h)
}
