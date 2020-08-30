#' Central quantile subspace
#'
#' \code{cqs} estimates the directions of the central quantile subspace.
#'
#' The function computes the directions that span the \eqn{\tau}th central
#' quantile subspace, i.e., the directions that define linear combinations of
#' the predictor \code{x} that contain all the information available on the
#' conditional quantile function.
#'
#' The function starts by estimating the initial vector, which is defined as the
#' least-squares estimator from regression the conditional quantile on \code{x}.
#' Then, if the dimension of the central quantile subspace is one, the algorithm
#' stops and reports that vector as the basis of the central quantile subspace.
#' Otherwise, the algorithm continues by creating more vectors and applying an
#' eigenvalue decomposition to extract linearly independent vectors.  If the
#' dimension of the central quantile subspace is unknown, it is estimated using
#' the modified-BIC type criterion of Zhu et al. (2010).
#'
#' @param x A design matrix.  The rows represent observations and the columns
#'   represent predictor variables.
#' @param y A vector of the response variable.
#' @param tau A quantile level, a number strictly between 0 and 1.
#' @param d The dimension of the central subspace.  If not specified, \code{d}
#'   is estimated using the modified-BIC type criterion of Zhu et al. (2010).
#' @param dtau The dimension of the central quantile subspace.  If not
#'   specified, \code{dtau} is estimated using the modified-BIC type criterion
#'   of Zhu et al. (2010)
#'
#' @return \code{cqs} computes the directions of the central quantile subspace,
#'   and returns:
#'   \itemize{
#'   \item{qvectors: }{The estimated directions of the
#'   \eqn{\tau} central quantile subspace, which, if \code{d_tau} is greater
#'   than 1, correspond to the eigenvalues of the matrix with column vectors the
#'   estimated vectors.}
#'
#'   \item{qvalues: }{The eigenvalues resulting from the eigenvalue decomposion
#'   of the matrix with column vectors the estimated vectors. If \code{d_tau} is
#'   one, the \code{evalues} output is not produced.}
#'
#'   \item{d: }{The dimension of the central subspace.  If not specified by the
#'   user, \code{d} is estimated using the modified-BIC type criterion of Zhu et
#'   al. (2010).}
#'
#'   \item{dtau: }{The dimension of the central quantile subspace.  If not
#'   specified by the user, \code{dtau_hat} is estimated using the modified-BIC
#'   type criterion of Zhu et al. (2010).}
#'   }
#' @references Zhu, L.-P., Zhu, L.-X., Feng, Z.-H. (2010) Dimension reduction in
#'   regression through cumulative slicing estimation. \emph{Journal of the
#'   American Statistical Association}, 105, 1455-1466.
#' @include llqr.R
#' @include bic_d.R
#'
#' @examples
#' # estimate the directions of a single-index model
#' set.seed(1234)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n * p), n, p); error <- rnorm(n)
#' y <- 3 * x[, 1] + x[, 2] + error
#' tau <- 0.5
#' out <- cqs(x, y, tau, d = 1, dtau = 1)
#' out
#' # without defining d and dtau
#' out <- cqs(x, y, tau)
#' out
#' out$qvectors[, 1:out$dtau]
#'
#' # estimate the directions of a multi-index model
#' set.seed(1234)
#' n <- 100; p <- 10
#' x <- matrix(rnorm(n * p), n, p); error <- rnorm(n)
#' y <- exp(x[, 1]) + x[, 2] + error
#' tau <- 0.5
#' cqs(x, y, tau)$qvectors[, 1:2]
#'
#' @export
cqs <- function(x, y, tau = 0.5, d, dtau) {

  # define the parameters
  n <- length(y); p <- dim(x)[2]
  # standardize the predictor variables
  xc <- scale(x, scale = FALSE)
  sig <- var(x)
  signrt <- MTS::msqrt(sig)$invsqrt
  #signrt <- matpower(sig, -1 / 2)
  xstand <- xc %*% signrt

  # use SIR for initial dimension reduction
  # apply BIC criterion if d is not specified
  if (missing(d)) {
    output <- dr::dr(y ~ xstand)
    lambdas <- output$evalues
    d_hat <- bic_d(lambdas, n, dim(x)[2])
    ahat <- cbind(output$evectors[, 1:d_hat])
    newx <- xstand %*% ahat
    d <- d_hat
  } else {
    ahat <- cbind(dr::dr(y ~ xstand)$evector[, 1:d])
    newx <- xstand %*% ahat
  }

  # define the bandwidth and estimate the conditional quantile
  h <- sd(y) * n^ (-1 / (d + 4))
  non_par <- llqr(newx, y, tau = tau, h = h)
  qhat <- non_par$ll_est

  # define the initial vector, i.e., the ordinary least squares estimator from
  # regressing qhat on x
  beta_hat <- lm(qhat ~ xstand)$coef[-1]

  # if dtau is missing, use the iterative procedure to produce all vectors
  # apply the BIC criterion to determine dtau
  if (missing(dtau)) {
    b <- matrix(0, p, p)
    b[, 1] <- beta_hat
    for (j in 2:(dim(x)[2])) {
      newx <- xstand %*% b[, j - 1]
      hatq <- llqr(newx, y, tau = tau, h = h)$ll_est
      mat <- matrix(0, dim(x)[1], dim(x)[2])
      for (i in 1:(dim(x)[1])) {
        mat[i, ] <- hatq[i] * xstand[i, ]
      }
      b[, j] <- apply(mat, 2, mean)
    }
    B <- b %*% t(b)
    eigenvalues <- eigen(B)$values
    out <- eigen(B)$vectors
    out <- signrt %*% out
    dtau <- bic_d(eigenvalues, n, dim(x)[2])
    list(qvectors = out, qvalues = eigenvalues, d = d, dtau = dtau)
  } else if (dtau > 1) {
    # if dtau is known to be greater than 1, then use the iterative procedure to
    # produce more vectors and select the eigenvectors associated with the dtau
    # largest eigenvalues
    b <- matrix(0, p, p)
    b[, 1] <- beta_hat
    for (j in 2:(dim(x)[2])) {
      newx <- xstand %*% b[, j - 1]
      hatq <- llqr(newx, y, tau = tau, h = h)$ll_est
      mat <- matrix(0, dim(x)[1], dim(x)[2])
      for (i in 1:(dim(x)[1])) {
        mat[i, ] <- hatq[i] * xstand[i, ]
      }
      b[, j] <- apply(mat, 2, mean)
    }
    B <- b %*% t(b)
    eigenvalues <- eigen(B)$values
    out <- eigen(B)$vectors
    out <- signrt %*% out
    list(qvectors = out, qvalues = eigenvalues, d = d, dtau = dtau)
  } else {
    # if dtau is known to be one, then the initial vector is sufficient
    out <- signrt %*% beta_hat
    out <- out / sqrt(sum(out^2))
    dtau <- dtau
    list(qvectors = out, d = d, dtau_hat = dtau)
  }
}
