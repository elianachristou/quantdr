#' Value-at-Risk estimation using the central quantile subspace
#'
#' \code{ValAR} estimates the one-step ahead \eqn{\tau}th Value-at-Risk for a
#' vector of returns.
#'
#' The function calculates the \eqn{\tau}th Value-at-Risk of the next time occurence,
#' i.e., that number such that the probability that the returns fall below its
#' negative value is \eqn{\tau}.  The parameter \eqn{\tau} is typically chosen to be a
#' small number such as 0.01, 0.025, or 0.05.  By definition, the negative value of the
#' \eqn{\tau}th Value-at-Risk is the \eqn{\tau}th conditional quantile.  Therefore,
#' the estimation is performed using a local linear conditional quantile estimation.
#' However, prior to this nonparametric estimation, a dimension reduction technique
#' is performed to select linear combinations of the predictor variables.
#'
#' Specifically, the user provides a vector of returns \code{y} (usually log-returns)
#' and an integer \code{p} for the number of past observations to be used as the
#' predictor variables.  The function then forms the n x p design matrix x, where n is
#' either the number of all returns (if the user wants to use all observations)
#' or the number of returns defined by the moving window (if the user provides
#' an integer for the moving window). Value-at-Risk is then defined as the negative
#' value of the \eqn{\tau}th conditional quantile of y given x.  However, to aid the
#' nonparametric estimation of the \eqn{\tau}th conditional quantile, the
#' \code{cqs} function is applied to estimate the fewest linear combinations of the
#' predictor \code{x} that contain all the information available on the conditional
#' quantile function.  Finally, the \code{llqr} function is applied to estimate the
#' local linear conditional quantile of y using the extracted directions as the
#' predictor variables.
#'
#' For more details on the method and for an application to the Bitcoin data, see
#' Christou (2020).  Moreover, see Christou and Grabchak (2019) for a thorough
#' comparison between the proposed methodology and commonly used methods.
#'
#' @param y A vector of returns (1 x n).
#' @param p An integer for the number of past observations to be used as
#'     predictor variables.  This will form the n x p design matrix.
#' @param tau A quantile level, a number strictly between 0 and 1. Commonly
#'     used choices are 0.01, 0.025, and 0.05.
#' @param movwind An optional integer number for the moving window.
#'     If not specified, all n observations will be used to fit the model.
#'     If specified, it should be an integer between p and n.  Typical values
#'     for moving windows are 250 or 500 (one or two years of return values).
#' @param chronological A logical operator to indicate whether the returns are
#'      in standard chronological order (from oldest to newest).  The default
#'      value is TRUE.  If the returns are in reverse chronological order, the
#'      function will rearrange them.
#'
#' @return \code{ValAR} returns the one-step ahead \eqn{\tau}th Value-at-Risk.
#'
#' @references Christou, E. (2020) Central quantile subspace. \emph{Statistics and
#' Computing}, 30, 677–695.
#'
#' Christou, E., Grabchak, M. (2019) Estimation of value-at-risk using single index
#' quantile regression.  \emph{Journal of Applied Statistics}, 46(13), 2418–2433.
#'
#' @examples
#' # estimate the one-step ahead Value-at-Risk without a moving window
#' data(edhec, package = "PerformanceAnalytics")
#' y <- as.vector(edhec[, 1]) # Convertible Arbitrage
#' p <- 5 # use the 5 most recent observations as predictor variables
#' tau <- 0.05
#' ValAR(y, p, tau) # the data is already in standard chronological order
#' # compare it with the historical Value-at-Risk calculation
#' PerformanceAnalytics::VaR(y, 0.95, method = 'historical')
#'
#' @export
#'

ValAR <- function(y, p, tau, movwind = NULL, chronological = TRUE){

  # compatibility checks
  # check whether y is a vector
  if (is.vector(y) == FALSE) {
    stop(paste("y needs to be a vector."))
  }
  # checks for NAs
  if (sum(is.na(y)) > 0) {
    stop(paste("Data include missing values."))
  }
  #  checks if n >  p
  if (length(y) <= p){
    stop(paste("number of observations of y (", length(y), ") should be  greater than p."))
  }
  # checks if the quantile level is one-dimensional
  if (length(tau) > 1) {
    stop(paste("quantile level needs to be one number."))
  }
  # checks if the quantile level is between 0 and 1 (strictly)
  if (tau >= 1 | tau <= 0) {
    stop(paste("quantile level needs to be a number strictly between 0 and 1."))
  }
  # checks if p is an integer
  if (p != round(p)) {
    stop(paste("p needs to be an integer."))
  }

  # check that the returns is an increasing order, i.e., from past to present
  # if not, reverse the vector of returns
  if (chronological == FALSE) {
    y <- rev(y)
  }

  # delete the last p observations, since the full data can start from the p + 1
  # observation. This will give the new vector of returns where a matrix X will exist
  newy <- y[-c(1:p)]
  n <- length(newy)

  # define the design matrix X, which is defined as the previous p observations
  # for each return
  X <- matrix(0, n, p)
  for (i in (p + 1):(n + p)){
    X[i - p, ] <- y[(i - p):(i - 1)]
  }

  # If a moving window is provided, define the new reduced vector of returns
  # and design matrix X, i.e., take the last movwind observations.
  if (is.null(movwind) == FALSE){

    # first, checks
    # checks if moving window is an integer
    if (movwind != round(movwind)) {
      stop(paste('Moving window needs to be an integer.'))
    }
    # checks if moving window is a number between p and n
    if (movwind > (length(y) - p) | movwind < p){
      stop(paste("Moving window needs to be greater than the number of covariates
    (", p, ") and less than the number of used observations n - p (", length(y) - p, ")"))
    }

    newy <- newy[(n - movwind + 1):n]
    X <- X[(n - movwind + 1):n, ]
    n <- length(newy)
  }

  # one-step ahead prediction
  out <- cqs(x = X, y = newy, tau = tau)
  beta_hat <- as.matrix(out$qvectors[, 1:out$dtau])
  newx = X %*% beta_hat
  x0 = as.vector(newy[(n - p + 1):n] %*% beta_hat)
  qhat_est <- as.numeric(llqr(newx, newy, tau, x0 = x0)$ll_est)

  return(VaR = qhat_est)
}
