#' Modified-BIC type criterion for structural dimension selection
#'
#' \code{bic_d} estimates the dimension of a subspace.
#'
#' The function estimates the dimension of a subspace using the modified-BIC
#' type criterion of Zhu et al. (2010).
#'
#' @param lambdas The eigenvalues of the matrix that spans the subspace
#' @param n The number of observations.
#'
#' @return \code{bic_d} returns the dimension of the subspace.
#'
#' @references Zhu, L.-P., Zhu, L.-X., and Feng, Z.-H. (2010) Dimension
#'   reduction in regression through cumulative slicing estimation.
#'   \emph{Journal of the American Statistical Association}, 105, 1455-1466.
#'
#' @export
#' @examples
#' set.seed(1234)
#' n <- 100; p <- 10; tau <- 0.5
#' x <- matrix(rnorm(n * p), n, p); error <- rnorm(n)
#' y <- 3 * x[, 1] + x[, 2] + error
#' out <- cqs(x, y, tau)
#' bic_d(out$qvalues, n)
#'
bic_d <- function(lambdas, n) {
  lambdas <- sort(lambdas, decreasing = TRUE)
  p <- length(lambdas)
  gn <- as.null(p)
  for (i in 1:length(lambdas)) {
    gn[i] <- n * sum((lambdas[1:i])^2) / sum((lambdas)^2) -
        (n^ (3 / 4) / p) * i * (i + 1) / 2
  }
  return(which(gn == max(gn)))
}
