#' Central Quantile Subspace
#'
#' \code{cqs} ...
#'
#' .....
#'
#' @param
#'
#' @include llqr.R
#' @include misc.R
#' @export
cqs <- function(x, y, tau = 0.5, d, dtau, h, method = "rule") {
  # define the parameters
  n <- length(y); p <- dim(x)[2]
  # standardize the predictor variables
  xc <- center(x)
  sig <- var(x)
  signrt <- matpower(sig, -1 / 2)
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

  # define the bandwidth, if missing.
  if (missing(h)) {
    if (method == "CV") {
      non_par <- llqr(newx, y, tau = tau)
      qhat <- non_par$ll_est; h <- non_par$h
    }
    if (method == "rule") {
      h <- KernSmooth::dpill(newx, y);
      h <- h * (tau * (1 - tau) / (dnorm(qnorm(tau)))^2)^.2
      non_par <- llqr(newx, y, tau = tau, h = h)
      qhat <- non_par$ll_est; h <- h
    }
  } else {
    non_par <- llqr(newx, y, tau = tau, h = h)
    qhat <- non_par$ll_est; h <- h
  }

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
      hatq <- llqr(newx, y, tau = tau, h = h, )$ll_est
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
    list(out = out, eigenvalues = eigenvalues, d = d, dtau_hat = dtau, h = h)
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
    list(out = out, eigenvalues = eigenvalues, d = d, dtau_hat = dtau, h = h)
  } else {
    # if dtau is known to be one, then the initial vector is sufficient
    out <- signrt %*% beta_hat
    out <- out / sqrt(sum(out^2))
    dtau <- dtau
    list(out = out, d = d, dtau_hat = dtau, h = h)
  }
}
