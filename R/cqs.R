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
#' least-squares estimator from regressing the conditional quantile on the
#' predictor variable \code{x}. Then, if the dimension of the central quantile
#' subspace is one, the algorithm stops and reports that vector as the basis of
#' the central quantile subspace. Otherwise, the algorithm continues by creating
#' more vectors and applying an eigenvalue decomposition to extract linearly
#' independent vectors.
#'
#' @param x A design matrix (n x p).  The rows represent observations and the
#'   columns represent predictor variables.
#' @param y A vector of the response variable.
#' @param tau A quantile level, a number strictly between 0 and 1.
#' @param dtau An optional dimension of the central quantile subspace.  If
#'   specified, it should be an integer between 1 and p, the number of columns
#'   of the design matrix \code{x}.  In the context of the algorithm, if
#'   \code{dtau} is known to be one, i.e., the assumed model is a single-index
#'   model, then the algorithm stops after estimating the initial vector and
#'   saves computational time.  However, if \code{dtau} is greater than one or
#'   (more realistically) unknown, then the algorithm continues on creating more
#'   vectors.
#'
#' @return \code{cqs} computes the directions of the central quantile subspace
#'   and returns:
#'   \itemize{
#'   \item{qvectors: }{The estimated directions of the
#'   \eqn{\tau}th central quantile subspace.}
#'
#'   \item{qvalues: }{The eigenvalues resulting from the eigenvalue decomposition
#'   of the matrix with column vectors that span the central quantile subspace.
#'   If \code{dtau} is one, the \code{qvalues} output is not produced.}
#'
#'   \item{dtau: }{Suggested dimension of the central quantile subspace.  If
#'   \code{dtau} is specified by the user then the algorithm outputs the
#'   user-defined value.  If \code{dtau} is not specified by the user then the
#'   algorithm outputs a suggested dimension using the modified-BIC type
#'   criterion of Zhu et al. (2010).  Note that this is one suggested method to
#'   estimate the structural dimension and is not necessarily a perfect one. The
#'   user has the option to use the eigenvalues \code{qvalues} on other
#'   criteria, like cross-validation, and determine the estimated dimension of
#'   the subspace.}}
#'
#' @references Zhu, L.-P., Zhu, L.-X., Feng, Z.-H. (2010) Dimension reduction in
#'   regression through cumulative slicing estimation. \emph{Journal of the
#'   American Statistical Association}, 105, 1455-1466.
#' @include llqr.R
#' @include bic_d.R
#'
#' @examples
#' # estimate the directions of a single-index model
#' set.seed(1234)
#' n <- 100
#' p <- 10
#' x <- matrix(rnorm(n * p), n, p)
#' error <- rnorm(n)
#' y <- 3 * x[, 1] + x[, 2] + error
#' tau <- 0.5
#' out <- cqs(x, y, tau, dtau = 1)
#' out
#' # without specifying dtau
#' out <- cqs(x, y, tau)
#' out
#' out$qvectors[, 1:out$dtau]
#'
#' @export
cqs <- function(x, y, tau = 0.5, dtau = NULL) {

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
  if (length(x) == length(y)) {
    stop("x is one dimensional, no need for dimension reduction.")
  }

  # define the parameters
  n <- length(y)
  p <- dim(x)[2]
  # standardize the predictor variables
  xc <- scale(x, scale = FALSE)
  sig <- var(x)
  signrt <- MTS::msqrt(sig)$invsqrt
  xstand <- xc %*% signrt

  # use SIR for initial dimension reduction
  # use bic_d to estimate d, the dimension of the central subspace
  output <- dr::dr(y ~ xstand)
  lambdas <- output$evalues
  d_hat <- bic_d(lambdas, n)
  ahat <- cbind(output$evectors[, 1:d_hat])
  newx <- xstand %*% ahat
  d <- d_hat

  # define the bandwidth and estimate the conditional quantile
  h <- max(n^(-1 / (d + 4)), min(2, sd(y)) * n^(- 1 / (d + 4)))
  non_par <- llqr(newx, y, tau = tau, h = h)
  qhat <- non_par$ll_est

  # define the initial vector, i.e., the ordinary least squares estimator from
  # regressing qhat on x
  beta_hat <- (solve(crossprod(cbind(1, xstand))) %*% crossprod(cbind(1, xstand), qhat))[-1]

  # if dtau is missing, use the iterative procedure to produce all vectors
  # apply the BIC criterion to determine dtau
  if (is.null(dtau)) {
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
    B <- tcrossprod(b, b)
    eigenvalues <- eigen(B)$values
    out <- eigen(B)$vectors
    out <- signrt %*% out
    dtau <- bic_d(eigenvalues, n)
    list(qvectors = out, qvalues = eigenvalues, dtau = dtau)
  } else if (dtau > 1) {
    # if dtau is known to be greater than 1, then use the iterative procedure to
    # produce more vectors
    # check that dtau is integer
    if (dtau<1 | dtau != as.integer(dtau) | dtau > p){
      stop(paste("dtau needs to be an integer between 1 and p (", dim(x)[2], ")."))
    }
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
    B <- tcrossprod(b, b)
    eigenvalues <- eigen(B)$values
    out <- eigen(B)$vectors
    out <- signrt %*% out
    list(qvectors = out, qvalues = eigenvalues, dtau = dtau)
  } else {
    # check that dtau is integer
    if (dtau<1 | dtau != as.integer(dtau) | dtau > p){
      stop(paste("dtau needs to be an integer between 1 and p (", dim(x)[2], ")."))
    }
    # if dtau is known to be one, then the initial vector is sufficient
    out <- signrt %*% beta_hat
    out <- out / sqrt(sum(out^2))
    dtau <- dtau
    list(qvectors = out, dtau = dtau)
  }
}
