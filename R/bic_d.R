#' Modified-BIC type criterion
#'
#' \code{bic_d} estimates the dimension of a subspace.
#'
#' The function estimates the dimension of a subspace using
#' the modified-BIC type criterion of Zhu et al. (2010).
#' The estimation is based on choosing \eqn{k} that maximizes
#' \eqn{G_n(k)=n \frac{\sum_{i=1}^{k} \hat{\lambda}_{i}^2}{\sum_{i=1}^{p}\hat{\lambda}_{i}^2}}
#'
#' @param lambdas The eigenvalues of the matrix that spans the subspace
#' @param n The number of observations.
#' @param p The

bic_d <- function(lambdas, n, p) {
  gn <- as.null(p)
  for (i in 1:length(lambdas)) {
    gn[i] <- n * sum((lambdas[1:i])^2) / sum((lambdas)^2) -
        2 * (n^ (3 / 4) / p) * i * (i + 1) / 2
  }
  return(which(gn == max(gn)))
}


