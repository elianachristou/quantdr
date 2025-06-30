#' Dimension Reduction Main Function
#'
#' \code{dr} performs dimension reduction for a multivariate response or scalar
#' response using a variety of dimension reduction techniques (e.g., SIR, SAVE,
#' etc.).  This function prepares the input data, handles missing values and
#' weights, and delegates the main computation to the internal \code{dr.compute}
#' function.
#'
#' @param formula A formula describing the model (e.g., \code{Y ~ X1 + X2 + ...}).
#' @param data A data frame containing the variables in the model.
#' @param subset An optional vector specifying a subset of observations to be used.
#' @param group An optional grouping variable used ini some dimension reduction
#'     methods.
#' @param na.action A function which indicates what should happen when the data
#'     contain NAs.  Default is \code{na.fail}.
#' @param weights Optional numeric vector of non-negaative weights.  If not specified,
#'     equal weights are used.
#' @param ... Additional arguments passed to \code{dr.compute}.
#'
#' @return A list of class \code{'dr'} containing the results of the dimension
#'     reduction procedure.  Components include:
#'     \item{call}{The matched call.}
#'   \item{y.name}{The name of the response variable.}
#'   \item{terms}{The model terms object.}
#'   \item{other elements}{As returned by the specific dimension reduction method
#'       in \code{dr.compute}.}
#' @noRd
dr <- function(formula, data, subset, group = NULL, na.action = na.fail,
               weights,...) {

  # Capture the function call and prepare model frame
  mf <- match.call(expand.dots = FALSE)
  mf$na.action <- na.action

  # Evaluate the group expression in the context of the data if it exists
  if(!is.null(mf$group)) {
   mf$group <- eval(parse(text = as.character(group)[2]), data,
                    environment(group)) }
  # Remoove unusued ... arguments for model.frame
  mf$... <- NULL
  mf[[1]] <- as.name("model.frame")

  # Evaluate the model frame with the current environment
  mf <- eval(mf, sys.frame(sys.parent()))

  # Extract model terms and response
  mt <- attr(mf, "terms")
  Y <- model.extract(mf, "response")
  X <- model.matrix(mt, mf)

  # Extract response name
  y.name <- if(is.matrix(Y)) colnames(Y) else as.character(attr(mt, "variables")[2])

  # Drop the intercept column from model matrix
  int <- match("(Intercept)", dimnames(X)[[2]], nomatch = 0)
  if (int > 0) X <- X[, -int, drop = FALSE]

  # Extract weights or assign equal weights
  weights <- mf$"(weights)"
  if(is.null(weights)) {
    weights <- rep(1, dim(X)[1])
    } else {
      if(any(weights < 0)) stop("Negative weights")
      pos.weights <- weights > 1.e-30
      weights <- weights[pos.weights]
      weights <- dim(X)[1] * weights / sum(weights) # normalize
      X <- X[pos.weights, ]
      Y <- if (is.matrix(Y)) Y[pos.weights, ] else Y[pos.weights]
    }

  # Run the main dimension reduction computation
  ans <- dr.compute(X, Y, weights = weights, group = mf$"(group)", ...)

  # Attach metadata and return
  ans$call <- match.call()
  ans$y.name <- y.name
  ans$terms <- mt
  ans
}

#' Internal computation function for dr
#'
#' Prepares and fits the initial dimension reduction object.  This function
#' is called internally by \code{dr()} and is not intended for direct use by end
#' users.
#'
#' @param x A numeric matrix of predictors (n x p).
#' @param y A response variable (vector or matrix), with the same number of rows
#'     as \code{x}.
#' @param weights A numeric vector of non-negative weights (length n).
#' @param group Optional grouping variable (e.g., for sliced inverse regression
#'     variants).
#' @param method A character string specifying the DR method to use (e.g., "sir",
#'     "save").
#' @param chi2approx A character string specifying the chi-square approximation
#'     to use (default: "bx").
#' @param ... Additional arguments passed to the method-specific fitting functions.
#'
#' @return An object of class \code{dr}, with additional method-specific subclass
#'     tags.
#' @noRd
dr.compute <- function(x, y, weights, group = NULL, method = "sir",
                       chi2approx = "bx", ...) {

  # Check that x and y have the same number of observations
  if (NROW(y) != nrow(x))
      stop("The response and predictors have different number of observations")

  # Warn if there are more predictors than observations
  if (NROW(y) < ncol(x))
    stop("The methods in dr require more observations than predictors")

  # Determine the class name of the object based on y and group
  classname <- if (is.matrix(y)) {
    c(paste("m", method, sep = ""), method)
    } else if (!is.null(group)) {
      c(paste("p", method, sep = ""), method)
      } else {
        method
      }

  genclassname <- "dr" # Generic class name for all objets from this framework

  # Compute square root of weights for use in weighted QR decomposition
  sweights <- sqrt(weights)

  # Apply weights to each column of x, then center (but do not scale) the matrix
  qrz <- qr(scale(apply(x, 2, function(a, sweights) a  * sweights, sweights),
                    center = TRUE, scale = FALSE))

  # Initialize the dr object with inputs and metadata
  ans <- list(x = x, y = y, weights = weights, method = method, cases = NROW(y),
                qr = qrz, group = group, chi2approx = chi2approx)

  # Assign class attributes (method-specific + general 'dr' class)
  class(ans) <-  c(classname, genclassname)

  # Fit the model using method-specific logic (delegted to dr.fit)
  ans <- dr.fit(object = ans, ...)

  # Reorder x columns according to QR pivoting and retain only the full-rank part
  ans$x <- ans$x[, qrz$pivot[1:qrz$rank]]

  # Return the compelted object
  ans
  }

#######################################################
#    accessor functions
#######################################################
#' Accessor for x componnent of a dr object
#' @noRd
dr.x <- function(object) {object$x}

#' Accessor for weights of a dr object
#' @noRd
dr.wts <- function(object) {object$weights}

#' Accessor for QR decomposition of a dr object
#' @noRd
dr.qr <- function(object) {object$qr}

#' Generic function for computing Q from a dr object
#' @noRd
#' @export
dr.Q <- function(object, ...){UseMethod("dr.Q")}

#' @method dr.Q default
#' @exportS3Method
dr.Q.default <- function(object) {
  qr.Q(dr.qr(object))[, 1:object$qr$rank]
  }

#' @noRd
#' @export
dr.R <- function(object) {UseMethod("dr.R")}

#' @method dr.R default
#' @exportS3Method
dr.R.default <- function(object) {
  qr.R(dr.qr(object))[1:object$qr$rank, 1:object$qr$rank]
}

#' Accessor for Z matrix from a dr object
#' @noRd
dr.z <- function(object) {
  sqrt(object$cases) * dr.Q(object)
}

#' Accessoor for response variable name from a dr object
#' @noRd
dr.yname <- function(object) {object$y.name}

#' Generic function for extracting basis vectoors from a dr object
#' @noRd
#' @export
dr.basis <- function(object, numdir, ...) {UseMethod("dr.basis")}

#' @method dr.basis default
#' @exportS3Method
dr.basis.default <- function(object, numdir = object$numdir) {
  object$evectors[, 1:numdir]
  }

#' @noRd
#' @export
dr.evalues <- function(object, ...) {UseMethod("dr.evalues")}

#' @method dr.evalues default
#' @exportS3Method
dr.evalues.default <- function(object, ...) object$evalues

#' Generic function for fitting a dimension reduction model
#'
#' This is the generic function for fitting a dimension reduction model.
#' Specific methods are dispatched based on the class of the object.
#' @noRd
#' @export
dr.fit <- function(object, numdir = 4, ...) UseMethod("dr.fit")

#' @method dr.fit default
#' @exportS3Method
dr.fit.default <- function(object, numdir = 4, ...) {

  # Compute the kernel matrix M (method-specific function)
  M <- dr.M(object, ...)

  # Compute eigenvalues/vectors based on shape of M
  D <- if (dim(M$M)[1] == dim(M$M)[2]) {
    eigen(M$M)
    } else {
    if (ncol(M$M) == 1) {
      eigen(M$M %*% t(M$M))
    } else {
        eigen(t(M$M) %*% M$M)
    }
    }

  # Order eigenvalues and corresponding vectors in decreasing magnitude
  or <- rev(order(abs(D$values)))
  evalues <- D$values[or]
  raw.evectors <- D$vectors[, or]

  # Transform eigenvectors using the inverse of weighted R matrix
  evectors <- backsolve(sqrt(object$cases) * dr.R(object), raw.evectors)
  evectors <- if (is.matrix(evectors)) evectors else matrix(evectors, ncol = 1)

  # Normalize each direction vector to unit norm
  evectors <- apply(evectors, 2, function(x) x / sqrt(sum(x^2)))

  # Name the rows and columns of the direction matrix
  names <- colnames(dr.x(object))[1:object$qr$rank]
  dimnames(evectors)<- list(names, paste("Dir", 1:NCOL(evectors), sep=""))

  # Construct the output object, including kernel matrix, directions, etc.
  aa <- c( object, list(evectors = evectors, evalues =evalues,
                      numdir = min(numdir, dim(evectors)[2], dim(M$M)[1]),
                      raw.evectors = raw.evectors), M)

  # Set the same class as input object and return
  class(aa) <- class(object)
  return(aa)
}

#####################################################################
# dr methods. Each method REQUIRES a dr.M function, and may also have
# a dr.y function and an function method.
#####################################################################

#####################################################################
# generic functions
#####################################################################
#' @noRd
#' @export
dr.M <- function(object, ...){UseMethod("dr.M")}

#' @noRd
#' @export
dr.y <- function(object) {UseMethod("dr.y")}

#' Default method for dr.y
#' @noRd
#' @method dr.y default
#' @exportS3Method
dr.y.default <- function(object, ...){object$y}

#' Generic function for performing dimension reduction test
#' @noRd
#' @export
dr.test <- function(object, numdir, ...){  UseMethod("dr.test")}

#' Default method for dr.test (returns NULL)
#' @noRd
#' @method dr.test default
#' @exportS3Method
dr.test.default <-function(object, numdir, ...) {NULL}

#' Generic function for coordinate-wise hypothesis tests
#' @noRd
#' @export
dr.coordinate.test <- function(object, hypothesis, d, chi2approx, ...) {
  UseMethod("dr.coordinate.test")
}

#' Default method for dr.coordinate.test (returns NULL)
#' @noRd
#' @method dr.coordinate.test default
#' @exportS3Method
dr.coordinate.test.default <- function(object, hypothesis, d, chi2approx, ...) {NULL}

#' Generic function for joint hypothesis testing
#' @noRd
#' @export
dr.joint.test <- function(object, ...){ UseMethod("dr.joint.test")}

#' Default method for dr.joint.test (returns NULL)
#' @noRd
#' @method dr.joint.test default
#' @exportS3Method
dr.joint.test.default <- function(object,...){NULL}

#####################################################################
#     OLS
#####################################################################
#' Computes the OLS kernel matri for dimension reduction
#'
#' This method computes the kernel matrix (M) for ordinary least squares (OLS),
#' which is a special case of dimension reduction with a univariate response.
#' @noRd
#' @method dr.M ols
#' @exportS3Method
dr.M.ols <- function(object, ...) {
  # Compute OLS projection of response onto predictors
  ols <- t(dr.z(object)) %*% (sqrt(dr.wts(object)) * dr.y(object))
  return(list(M = ols %*% t(ols), numdir = 1))}

#####################################################################
#     Sliced Inverse Regression (SIR) and Multivariate SIR Functions
#####################################################################
#' Compute kernel matrix M for SIR method
#'
#' This function calculates the kernel matrix (M) used in Sliced Inverse
#' Regression (SIR) by computing weighted means of the predictor transformations
#' within response slices.
#'
#' @noRd
#' @method dr.M sir
#' @exportS3Method
dr.M.sir <-function(object, nslices = NULL, slice.function = dr.slices,
                    sel = NULL, ...) {

  # Select subset of observations of all if sel is NULL
  sel <- if (is.null(sel)) 1:dim(dr.z(object))[1] else sel

  # Extract the predictor matrix and response vector/matrix
  z <- dr.z(object)[sel, , drop = FALSE]
  y <- dr.y(object)
  y <- if (is.matrix(y)) y[sel, , drop = FALSE] else y[sel]

  # Determine number of sllices if not specified
  h <- if (!is.null(nslices)) nslices else max(8, NCOL(z)+3)

  # Create slices based on the response
  slices <- slice.function(y, h)


  # Initialize matrix for weighted means and slice weights
  zmeans <- matrix(0, nrow = slices$nslices, ncol = NCOL(z))
  slice.weight <- rep(0, slices$nslices)

  # Extract observation weights
  wts <- dr.wts(object)

  # Weighted mean function
  wmean <- function (x, wts) sum(x * wts) / sum (wts)

  # Compute weighted means within each slice
  for (j in 1:slices$nslices) {
    sel_slice <- slices$slice.indicator==j
    zmeans[j, ]<- apply(z[sel_slice, , drop = FALSE], 2, wmean, wts = wts[sel_slice])
    slice.weight[j]<-sum(wts[sel_slice])
  }

  # Calculate kernel matrix M using weighted slice means
  M <- t(zmeans) %*% apply(zmeans, 2, "*", slice.weight) / sum(slice.weight)

  # Return list with kernel matrix and slice info
  return (list (M = M, slice.info = slices))
}

#' Alias of SIR kernel matrix function for multivariate SIR (MSIR)
#' @noRd
#' @method dr.M msir
#' @exportS3Method
dr.M.msir <-function(...) {dr.M.sir(...)}

#' Compute test statistics for SIR dimension reduction
#'
#' This function calculates the SIR test statistics for the first \code{numdir}
#' directions, including chi-square p-values.
#' @noRd
#' @method dr.test sir
#' @exportS3Method
dr.test.sir<-function(object, numdir = object$numdir, ...) {
  e <- sort(object$evalues) # sorted eigenvalues
  p <- length(object$evalues) # number of eigenvalues
  n <- object$cases # number of cases
  st <- df <- pv <- 0
  nt <- min(p,numdir)

  # Compute test statistics, degrees of freedom, and p-values for each dimension
  for (i in 0:(nt-1)) {
    st[i + 1] <- n * (p - i) * mean(e[seq(1, p - i)])
    df[i + 1] <- (p - i) * sum(object$slice.info$nslices - i - 1)
    pv[i + 1] <- 1 - pchisq(st[i + 1], df[i + 1])
  }

  # Create labeled data frame for output
  z <- data.frame(cbind(st, df, pv))
  rr <- paste(0:(nt - 1), "D vs >= ", 1:nt, "D", sep = "")
  dimnames(z) <- list(rr, c("Stat","df","p.value"))
  z
}

#  Written by Yongwu Shao, May 2006
#' Coordinate test for SIR dimension reduction
#'
#' Performs the coordinate test for dimension reduction using the SIR method.
#' It tests a hypothesis on the dimension reduction subspace specified by
#' \code{hypothesis}.
#'
#' @noRd
#' @method dr.coordinate.test sir
#' @exportS3Method
dr.coordinate.test.sir <- function(object, hypothesis, d = NULL,
                                 chi2approx = object$chi2approx,
                                 pval = "general", ...) {

  # Convert hypothesis to matrix if formula is given
  gamma <- if (inherits(hypothesis, "formula"))
    coord.hyp.basis(object, hypothesis)
  else as.matrix(hypothesis)

  p <- length(object$evalues)
  n <- object$cases
  z <- dr.z(object)
  ev <- object$evalues
  slices <- object$slice.info
  h <- slices$nslices

  # Number of directions of defaults to min(h, length of eigenvalues)
  d <- if(is.null(d)) min(h,length(ev)) else d

  M <- object$M
  r <- p - dim(gamma)[2]

  # Compute the matrix H and its QR decomposition (complete)
  H <- (dr.R(object)) %*% gamma
  H <- qr.Q(qr(H), complete = TRUE) # a p x p orthogonal matrix

  # Partition H into two parts
  QH <- H[, 1:(p - r), drop = FALSE] # first p-r columns
  H <- H[, (p - r + 1):p, drop = FALSE] # last r columns

  # Test statistic based on eigenvalues
  st <- n * sum(ev[1:d]) - n * sum(sort(eigen(t(QH) %*% M %*% QH)$values,
                                decreasing = TRUE)[1:min(d,p - r)])

  # Weights for test statistic
  wts <- 1 - ev[1:min(d, h - 1)]

  # Restricted test p-value
  testr <- dr.pvalue(rep(wts, r), st, chi2approx = chi2approx)

  # General test calculation
  epsilon <- array(0, c(n, h))
  zmeans <- array(0, c(p, h))

  for (i in 1:h) {
    sel <- (slices$slice.indicator == i)
    f_k <- sum(sel) / n
    zmeans[, i] <- apply(z[sel, , drop = FALSE], 2, mean)
    epsilon[, i] <- (sel - f_k - z %*% zmeans[, i] * f_k) / sqrt(f_k)
  }

  HZ <- z %*% H
  Phi <- svd(zmeans, nv = h)$v[, 1:min(d,h), drop = FALSE]

  epsilonHZ <- array(0, c(n, r * d))
  for (j in 1:n) epsilonHZ[j, ] <- t(Phi)% *%t (t(epsilon[j, ])) %*% t(HZ[j, ])

  wts <- eigen(((n - 1) / n) * cov(epsilonHZ))$values
  testg <- dr.pvalue(wts[wts > 0], st, chi2approx = chi2approx)

  # Select p-value according to user input
  pv <- if (pval == "restricted") testg$pval.adj else testr$pval.adj

  # Return test statistic and p-value as a data frame
  z <- data.frame(cbind(st, pv))
  dimnames(z) <- list("Test", c("Statistic", "P.value"))
  z
}

#####################################################################
#     Sliced Average Variance Estimation
#  Original by S. Weisberg; modified by Yongwu Shao 4/27/2006
#  to compute the A array needed for dimension tests
#####################################################################
#' SAVE kernel matrix computation
#'
#' Computes the kernel matrix \eqn{M} used in Sliced Average Variance Estimation
#' (SAVE), along with a 3D array \eqn{A} used for coordinate hypothesis testing.
#'
#' @param object A fitted \code{dr} object.
#' @param nslices Number of slices.  Default is \code{max(8, p + 3)} where \code{p}
#' is the number of predictors.
#' @param slice.function A function to slice the response.  Default is \code{dr.slices}.
#' @param sel Optional index vector selecting a subset of observations to use.
#' @param ... Additional arguments passed to the slicing function.
#'
#' @return A list with components:
#' \item{M}{The SAVE kerel matrix.}
#' \item{A}{A 3D array used for testing coordinate hypotheses.}
#' \item{slice.info}{Information about the slicing.}
#'
#' @noRd
#' @method dr.M save
#' @exportS3Method
dr.M.save <- function(object, nslices = NULL, slice.function = dr.slices,
                      sel=NULL, ...) {
  # Subset observations if specified
  sel <- if(is.null(sel)) 1:dim(dr.z(object))[1] else sel
  z <- dr.z(object)[sel, ]
  y <- dr.y(object)
  y <- if (is.matrix(y)) y[sel, ] else y[sel]
  wts <- dr.wts(object)[sel]

  # Determine number of slices
  h <- if (!is.null(nslices)) nslices else max(8, ncol(z) + 3)
  slices <- slice.function(y, h)

  # Initialize matrices
  M <- matrix(0, NCOL(z), NCOL(z)) # SAVE kernel matrix
  A <- array(0, c(slices$nslices, NCOL(z), NCOL(z))) # Slice-wise centered covariances
  ws <- rep(0, slices$nslices) # Slice weights

  # Helper: weighted covariance estimate
  wvar <- function(x, w) {
    (if (sum(w) > 1) {(length(w) - 1) / (sum(w) - 0)} else {0}) *
      var(sweep(x, 1, sqrt(w), "*"))}

  # Loop through slices and compute SAVE components
  for (j in 1:slices$nslices) {
    ind <- slices$slice.indicator == j
    IminusC <- diag(rep(1, NCOL(z))) - wvar(z[ind, ], wts[ind])
    ws[j] <- sum(wts[ind])
    A[j,,] <- sqrt(ws[j]) * IminusC
    M <- M + ws[j] * IminusC %*% IminusC
  }

  M <- M / sum(ws)
  A <- A / sqrt(sum(ws)) # Normalize A for use in tests

  return(list(M = M, A = A, slice.info = slices))
}

# Written by Yongwu Shao, 4/27/2006
#' Dimension test for SAVE
#'
#' Performs booth normal and general dimension tests using the SAVE method.
#'
#' @param object A fitted \code{dr} object with method \code{save}.
#' @param numdir Number of directions to test.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame containing test statistics and p-values under both the
#' normal and general assumptions.
#'
#' @noRd
#' @method dr.test save
#' @exportS3Method
dr.test.save <- function(object, numdir = object$numdir, ...) {
  p <- length(object$evalues)
  n <- object$cases
  h <- object$slice.info$nslices
  A <- object$A
  M <- object$M

  # Eigen decomposition of SAVE matrix M
  D <- eigen(M)
  or <- rev(order(abs(D$values)))
  evectors <- D$vectors[, or]

  # Initialize results
  st.normal <- df.normal <- st.general <- df.general <- 0
  pv.normal <- pv.general <- 0
  nt <- numdir

  for (i in 1:(nt -  1)) {
    theta <- evectors[, (i + 1):p, drop = FALSE]

    # Normal test statistic
    st.normal[i + 1] <- 0
    for (j in 1:h) {
      st.normal[i + 1] <- st.normal[i + 1] +
        sum((t(theta) %*% A[j, , ] %*% theta)^2) * n / 2
    }

    df.normal[i + 1] <- (h - 1) * (p - i) * (p - i + 1) / 2
    pv.normal[i + 1] <- 1 - pchisq(st.normal[i + 1], df.normal[i + 1])

    # General test statistic
    HZ <- dr.z(object) %*% theta
    ZZ <- array(0, c(n, (p - i) * (p - i)))
    for (j in 1:n) ZZ[j, ] <- t(t(HZ[j, ])) %*% t(HZ[j, ])
    Sigma <- cov(ZZ) / 2
    df.general[i + 1] <- sum(diag(Sigma))^2 / sum(Sigma^2) * (h - 1)
    st.general[i + 1] <- st.normal[i + 1] * sum(diag(Sigma)) / sum(Sigma^2)
    pv.general[i + 1] <- 1 - pchisq(st.general[i + 1], df.general[i + 1])
  }

  z <- data.frame(cbind(st.normal, df.normal, pv.normal, pv.general))
  rr <- paste(0:(nt - 1), "D vs >= ", 1:nt, "D", sep = "")
  dimnames(z) <- list(rr, c("Stat", "df(Nor)", "p.value(Nor)", "p.value(Gen)"))
  z
}

# Written by Yongwu Shao, 4/27/2006
#' Coordinate test for SAVE
#'
#' Performs a hypothesis test for coordinate subspaces using the SAVE method.
#'
#' @param object A fitted \code{dr} object using the SAVE method.
#' @param hypothesis A matrix or formula specifying the hypothesized coordinate
#'     subspace.
#' @param d Number of directions to test.  If \code{NIULL}, selected automatically.
#' @param chi2approx Method for chi-square approximation (used for general test).
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with test statistic, degrees of freedom, and p-values.
#'
#' @noRd
#' @method dr.coordinate.test save
#' @exportS3Method
dr.coordinate.test.save <- function (object, hypothesis, d = NULL,
                                     chi2approx = object$chi2approx, ...) {
  # Construct coordinate hypothesis basis
  gamma <- if (inherits(hypothesis, "formula"))
    coord.hyp.basis(object, hypothesis)
  else as.matrix(hypothesis)

  p <- length(object$evalues)
  n <- object$cases
  h <- object$slice.info$nslices
  A <- object$A

  st <- df <- pv <- 0
  gamma <- (dr.R(object)) %*% gamma
  r <- p - dim(gamma)[2]
  H <- as.matrix(qr.Q(qr(gamma), complete = TRUE)[, (p - r + 1):p])

  # Normal-theory test statistic
  st <- 0
  for (j in 1:h) {
    st <- st + sum((t(H) %*% A[j, , ] %*% H)^2) * n / 2
  }
  df.normal <- (h - 1) * r * (r + 1) / 2
  pv.normal  <- 1 - pchisq(st, df.normal)

  # General test statistic
  {
    HZ <- dr.z(object) %*% H
    ZZ <- array(0, c(n, r^2))
    for (j in 1:n) {
      ZZ[j, ] <- t(t(HZ[j, ])) %*% t(HZ[j, ])
    }
    wts <- rep(eigen(((n - 1) / n) * cov(ZZ) / 2)$values, h - 1)
    testg <- dr.pvalue(wts[wts > 0], st, a = chi2approx)
  }

  z <- data.frame(cbind(st, df.normal, pv.normal, testg$pval.adj))
  dimnames(z) <- list("Test", c("Statistic", "df(Nor)", "p.val(Nor)",
                                  "p.val(Gen)"))
    z
  }

#####################################################################
# Principal Hessian Directions (pHd), pHdy and pHdres Methods
#####################################################################
#' @method dr.M phdy
#' @exportS3Method
dr.M.phdy <- function(...) {dr.M.phd(...)}

#' @method dr.M mphd
#' @exportS3Method
dr.M.mphd <- function(...) stop("Multivariate pHd not implemented!")

#' @method dr.M phdres
#' @exportS3Method
dr.M.phdres <- function(...) {dr.M.phd(...)}

#' @method dr.M mphdres
#' @exportS3Method
dr.M.mphdres <- function(...) stop("Multivariate pHd not implemented!")

#' @method dr.M mphy
#' @exportS3Method
dr.M.mphy <- function(...) stop("Multivariate pHd not implemented!")

#' @method dr.M phd
#' @exportS3Method
dr.M.phd <-function(object,...) {
  # Computes the kernel matrix M for pHd method
  wts <- dr.wts(object) # Observation weights
  z <- dr.z(object) # Standardized predictors
  y <- dr.y(object) # Response
  # Kernell matrix as weighted covariance of y*z with z
  M<- (t(apply(z, 2, "*", wts * y)) %*% z) / sum(wts)
  return(list(M = M))
}

#' @method dr.y phdy
#' @exportS3Method
dr.y.phdy <- function(object) {
  # Centered response for phdy
  y <- object$y
  y - mean(y)
  }

#' @method dr.y phdres
#' @exportS3Method
dr.y.phdres <- function(object) {
  # Residuals from regressing y on the predictors (using QR)
  y <- object$y
  sw <- sqrt(object$weights)
  qr.resid(object$qr,sw*(y-mean(y)))
}

#' @method dr.y phd
#' @exportS3Method
dr.y.phd <- function(object) {dr.y.phdres(object)}

# Modified by Jorge de la Vega, February, 2001
#' @method dr.test phd
#' @exportS3Method
dr.test.phd<-function(object, numdir = object$numdir, ...) {
  # pHd asymptotic test based on OLS residuals
  e <- sort(abs(object$evalues)) # sorted eigenvalues
  p <- length(object$evalues)
  resi <- dr.y(object) # residuals
  varres <- 2 * var(resi) # variance under null
  n <- object$cases
  st <- df <- pv <- 0
  nt <- min(p, numdir)

  for (i in 0:(nt - 1)) {
  # the statistic is partial sums of squared absolute eigenvalues
  st[i + 1] <- n * (p - i) * mean(e[seq(1, p - i)]^2) / varres
  # compute the degrees of freedom
  df[i + 1] <- (p - i) * (p - i + 1) / 2
  # use asymptotic chi-square distribution for p.values.  Requires normality.
  pv[i + 1] <- 1 - pchisq(st[i + 1], df[i + 1])
  }

  # Additional tests
  indep <- dr.indep.test.phdres(object, st[1])
  lc <- dr.test2.phdres(object, st)

  # report results
  z <- data.frame(cbind(st, df, pv, c(indep[[2]], rep(NA, length(st) - 1)), lc[, 2]))
  rr <- paste(0:(nt - 1),"D vs >= ", 1:nt, "D", sep = "")
  cc <- c("Stat", "df", "Normal theory", "Indep. test", "General theory")
  dimnames(z) <- list(rr, cc)
  z
}

#' @method dr.test phdres
#' @exportS3Method
dr.test.phdres <- function(object, numdir, ...) {
  dr.test.phd(object, numdir)
  }

#' @method dr.test phdy
#' @exportS3Method
dr.test.phdy <- function(object, numdir, ...) {
  # pHd test using y directly (requires normal predictors)
  e <- sort(abs(object$evalues))
  p <- length(object$evalues)
  # get the response
  resi <- dr.y(object)
  varres <- 2 * var(resi)
  n <- object$cases
  st <- df <- pv <- 0
  nt <- min(p, numdir)

  for (i in 0:(nt - 1)) {
  # the statistic is partial sums of squared absolute eigenvalues
  st[i + 1] <- n * (p - i) * mean(e[seq(1, p - i)]^2) / varres
  # compute the degrees of freedom
  df[i + 1] <- (p - i) * (p - i + 1) / 2
  # use asymptotic chi-square distribution for p.values.  Requires normality.
  pv[i + 1] <- 1 - pchisq(st[i + 1], df[i + 1])
  }

  # report results
  z <- data.frame(cbind(st, df, pv))
  rr <- paste(0:(nt - 1), "D vs >= ", 1:nt, "D", sep = "")
  cc<-c("Stat","df","p.value")
  dimnames(z) <- list(rr, cc)
  z
  }

#####################################################################
# pHdq method for sufficient dimension reduction
# Reference: Li (1992, JASA)
# Quadratic form based on full second-order regression fit
# Corrected by Jorge de la Vega 7/10/01
#####################################################################
#' @method dr.M phdq
#' @exportS3Method
dr.M.phdq <- function(object, ...) {
  # Extracts parameters from a full quadratic regression
  pars <- fullquad.fit(object)$coef
  k <- length(pars)

  # Solve quadratic for p (number of predictors):
  # k = 1 + 2p + p(p-1)/2 --> quadratic equation
  p <- (-3 + sqrt(9 + 8 * (k - 1))) / 2
  p <- round(p) # Round in case of numeric precision error

  # Initialize matrix of second-order coefficients (diagonal)
  mymatrix <- diag(pars[(p + 2):(2 * p + 1)])
  pars <- pars[-(1:(2 * p + 1))]

  # Fill off-diagonal entries symetrically
  for (i in 1:(p - 1)) {
    mymatrix[i, (i + 1):p] <- pars[1:(p - i)] / 2
    mymatrix[(i + 1):p, i] <- pars[1:(p - i)] / 2
    pars <- pars[-(1:(p - i))]
  }
  return(list(M = mymatrix))
}

#' Full Quadratic Fit
#'
#' Fits a full quadratic model (including linear, squared, and cross-product terms)
#' to the response using weighted least squares.
#'
#' @param object A fitted \code{dr} object.
#' @return A weighted linear model object with quadratic terms.
#' @noRd
fullquad.fit <-function(object) {
  x <- dr.z(object) # Standardized predictors
  y <- object$y # Response
  w <- dr.wts(object) # Weights
  z <- cbind(x, x^2) # Linear and squared terms

  p <- NCOL(x)
  # Add cross-product terms
  for (j in 1:(p - 1)) {
    for (k in (j + 1):p) {
      z <- cbind(z, matrix(x[, j] * x[, k], ncol = 1))
    }
  }

  # Fit the full quadratic model
  lm(y ~ z, weights = w)
}

#' @method dr.y phdq
#' @exportS3Method
dr.y.phdq <- function(object) {
  # Returns the Pearson residuals from full quadratic fit
  residuals(fullquad.fit(object), type = "pearson")
}

#' @method dr.test phdq
#' @exportS3Method
dr.test.phdq <- function(object, numdir, ...){
  # Reuses the test procedure for pHd
  dr.test.phd(object, numdir)
}

#' @method dr.M mphdq
#' @exportS3Method
dr.M.mphdq <- function(...) stop("Multivariate pHd not implemented!")

#####################################################################
# Auxiliary Functions and Helpers
#####################################################################
#####################################################################
#  Auxiliary function to find matrix H used in marginal coordinate tests
#  We test Y indep X2 | X1, and H is a basis in R^p for the column
#  space of X2.
#####################################################################
#' Basis for coordinate hypothesis
#'
#' This function constructs a basis matrix H used in marginal coordinate tests,
#' identifying the subspace specified in 'spec'.
#'
#' @param object A fitted \code{dr} object.
#' @param spec A formula indicating the variables to test.
#' @param which Index multiplier (usually 1); allows for selection.
#'
#' @return A basis matrix for the coordinate hypothesis
#'
#' @noRd
coord.hyp.basis <- function(object, spec, which = 1) {
  mod2 <- update(object$terms, spec)
  Base <- attr(object$terms, "term.labels")
  New  <- attr(terms(mod2), "term.labels")
  cols <- na.omit(match(New, Base))
  if(length(cols) != length(New)) stop("Error---bad value of 'spec'")
  as.matrix(diag(rep(1, length(Base)))[, which * cols])
}

#####################################################################
# Direction Recovery Methods
#####################################################################
#' Generic dispatcher for direction recovery
#' @noRd
dr.directions <- function(object, which, x) {UseMethod("dr.direction")}

#' Generic method for direction recovery
#' @noRd
dr.direction  <- function(object, which, x) {UseMethod("dr.direction")}

#' Default method to compute directions
#'
#' This function standardizes predictor values and applied the estimated basis
#' vectors.
#' @noRd
#' @method dr.direction default
#' @exportS3Method
dr.direction.default <- function(object, which = NULL, x = dr.x(object)) {
    ans <- (apply(x, 2, function(x) {x - mean(x)}) %*% object$evectors)
    which <- if (is.null(which)) seq(dim(ans)[2]) else which
    ans <- apply(ans, 2, function(x) if(x[1] <= 0)  -x else x)[, which]
    if (length(which) > 1) {
      dimnames(ans) <- list ( attr(x,"dimnames")[[1]],
                              paste("Dir", which, sep = ""))
    }
    ans
}

#####################################################################
# Plotting methods
#####################################################################
#' Plot method for dr objects
#'
#' This function produces pairwise plots of the estimated dimension reduction
#' directions.
#' @noRd
#' @method plot dr
#' @exportS3Method
plot.dr <- function(x, which = 1:x$numdir, mark.by.y = FALSE,
                    plot.method = pairs, ...) {
  d <- dr.direction(x, which)
  if (mark.by.y == FALSE) {
    plot.method(cbind(dr.y(x), d), labels = c(dr.yname(x), colnames(d)), ...)
  } else {
    plot.method(d, labels = colnames(d), col = markby(dr.y(x)), ...)
    }
}

#' Mark observations by categorical response
#'
#' Helper function to assign colors or symbols based on response levels.
#' @noRd
markby <- function (z, use = "color", values = NULL, color.fn = rainbow,
            na.action = "na.use")
  {
    u <- unique(z)
    lu <- length(u)
    ans <- 0
    vals <- if (use == "color") {
      if (!is.null(values) && length(values) == lu)
        values
      else color.fn(lu)
    }
    else {
      if (!is.null(values) && length(values) == lu)
        values
      else 1:lu
    }
    for (j in 1:lu) if (is.na(u[j])) {
      ans[which(is.na(z))] <- if (na.action == "na.use")
        vals[j]
      else NA
    }
    else {
      ans[z == u[j]] <- vals[j]
    }
    ans
  }

###################################################################
# Print Method for dr Object
###################################################################
#' Print method for dr object
#'
#' This function prints estimated basis vectors and eigenvalues.
#' @noRd
#' @method print dr
#' @exportS3Method
print.dr <- function(x, digits = max(3, getOption("digits") - 3), width = 50,
           numdir = x$numdir, ...) {
    fout <- deparse(x$call, width.cutoff = width)
    for (f in fout) cat("\n",f)
    cat("\n")
    cat("Estimated Basis Vectors for Central Subspace:\n")
    evectors <- x$evectors
    print.default(evectors[, 1:min(numdir, dim(evectors)[2])])
    cat("Eigenvalues:\n")
    print.default(x$evalues[1:min(numdir, dim(evectors)[2])])
    cat("\n")
    invisible(x)
  }

###################################################################
# Summary Method for dr Object
###################################################################
#' Summary method for dr object
#'
#' This function summarizes reduction results including eigenvalues and angles.
#' @noRd
#' @method summary dr
#' @exportS3Method
summary.dr <- function (object, ...) {
  z <- object
  ans <- z[c("call")]
  nd <- min(z$numdir, length(which(abs(z$evalues) > 1.e-15)))
  ans$evectors <- z$evectors[, 1:nd]
  ans$method <- z$method
  ans$nslices <- z$slice.info$nslices
  ans$sizes <- z$slice.info$slice.sizes
  ans$weights <- dr.wts(z)
  sw <- sqrt(ans$weights)
  y <- z$y
  ans$n <- z$cases
  ols.fit <- qr.fitted(object$qr, sw * (y - mean(y)))
  angles <- cosangle(dr.direction(object), ols.fit)
  angles <- if (is.matrix(angles)) angles[, 1:nd] else angles[1:nd]
  if (is.matrix(angles)) dimnames(angles)[[1]] <- z$y.name
  angle.names <- if (!is.matrix(angles)) "R^2(OLS|dr)" else
  paste("R^2(OLS for ",dimnames(angles)[[1]], "|dr)", sep = "")
  ans$evalues <-rbind (z$evalues[1:nd],angles)
  dimnames(ans$evalues)<- list(c("Eigenvalues", angle.names),
       paste("Dir", 1:NCOL(ans$evalues), sep = ""))
ans$test <- dr.test(object,nd)
class(ans) <- "summary.dr"
ans
}

###################################################################
# Print Method for summary.dr
###################################################################
#' Print method for summary.dr object
#'
#' This function prints summary of dimension reduction analysis incuding test
#' results.
#' @noRd
#' @method print summary.dr
#' @exportS3Method
print.summary.dr <- function (x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nCall:\n")#S: ' ' instead of '\n'
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Method:\n")#S: ' ' instead of '\n'
    if(is.null(x$nslices)) {
      cat(paste(x$method, ", n = ", x$n,sep=""))
      if(diff(range(x$weights)) > 0)
        cat(", using weights.\n") else cat(".\n")}
    else {
      cat(paste(x$method," with ",x$nslices, " slices, n = ",
                x$n,sep=""))
      if(diff(range(x$weights)) > 0)
        cat(", using weights.\n") else cat(".\n")
      cat("\nSlice Sizes:\n")#S: ' ' instead of '\n'
      cat(x$sizes, "\n")}
    cat("\nEstimated Basis Vectors for Central Subspace:\n")
    print(x$evectors, digits = digits)
    cat("\n")
    print(x$evalues,digits=digits)
    if (length(x$omitted) > 0){
      cat("\nModel matrix is not full rank.  Deleted columns:\n")
      cat(x$omitted, "\n")}
    if (!is.null(x$test)){
      cat("\nLarge-sample Marginal Dimension Tests:\n")
      print(as.matrix(x$test), digits = digits)}
    invisible(x)
  }

##################################################################
## Translations of pHd test methods from Arc to R
##  Original Lisp functions by R. D. Cook
##  Translated to R by S. Weisberg, February, 2001
##################################################################
#' Covariance matrix for eW sttistics
#'
#' This function computes the matrices W and eW described in Sec. 12.3.1 of
#' Cook (1998).
#'
#' @param object A fitted \code{dr} object.
#' @param scaled Logical; whether to center y by weighted mean.
#'
#' @return A covariance matrix used in pHd tests.
#' @noRd
cov.ew.matrix <- function(object, scaled = FALSE) {
  mat.normalize <- function(a) apply(a, 2, function(x) {x / (sqrt(sum(x^2)))})
  n <- dim(dr.x(object))[1]
  TEMPwts <- object$weights
  sTEMPwts <- sqrt(TEMPwts)
  v <- sqrt(n) * mat.normalize(
    apply(scale(dr.x(object), center = TRUE, scale = FALSE), 2, "*", sTEMPwts) %*%
      object$evectors)
  y <- dr.y(object) # get the response
  y <- if (scaled) y - mean(sTEMPwts * y) else 1 # a multiplier in the matrix
  p <- dim(v)[2]
  ew0 <- NULL
  for (i in 1:p) {
    for (j in i:p) {
      ew0 <- cbind(ew0, if (i == j) y * (v[, j]^2 - 1) else y * sqrt(2) * v[, i] * v[, j])
    }
  }
  wmean <- function (x, w) sum(x * w) / sum (w)
  tmp <- apply(ew0, 2, function(x, w, wmean) sqrt(w) * (x - wmean(x, w)),
               TEMPwts, wmean)
  ans <- (1 / sum(TEMPwts)) * t(tmp) %*% tmp
  ans
}

#' General p-values for phd using eW matrix
#'
#' This function computes adjusted p-values using eigenvalues of the eW covariance
#' matrix.
#'
#' @param object A fitted \code{dr} object.
#' @param stats Vector of test statistics.
#'
#' @return A data frame with statistics and general p-values.
#' @noRd
#' @export
dr.test2.phdres <- function(object, stats) {
  covew <- cov.ew.matrix(object, scaled = TRUE)
  C <- 0.5 / var(dr.y(object))
  p <- length(stats)
  pval <- NULL
  d2 <- dim(dr.x(object))[2]
  start <- -d2
  end <- dim(covew)[2]
  for (i in 1:p) {
    start <- start + d2 - i + 2
    evals <- eigen(as.numeric(C) * covew[start:end, start:end],
                   only.values = TRUE)$values
    pval <- c(pval, dr.pvalue(evals, stats[i],
                              chi2approx = object$chi2approx)$pval.adj)
  }

  # report results
  z <- data.frame(cbind(stats, pval))
  rr <- paste(0:(p - 1), "D vs >= ", 1:p, "D", sep = "")
  dimnames(z) <- list(rr, c("Stat", "p.value"))
  z
}

#' Independence test for phd residuals
#'
#' This function tests independence using the uncentered eW matrix.
#'
#' @param object A fitted \code{dr} object.
#' @param stat The test statistic.
#'
#' @return A data frame with test result.
#' @noRd
#' @export
dr.indep.test.phdres <- function(object, stat) {
  eval <- eigen(cov.ew.matrix(object, scaled = FALSE), only.values = TRUE)
  pval <- dr.pvalue(0.5 * eval$values, stat,
                  chi2approx = object$chi2approx)$pval.adj

  # report results
  z <- data.frame(cbind(stat, pval))
  dimnames(z) <- list(c("Test of independence"), c("Stat", "p.value"))
  z
}

#' Compute p-value using Bentler-Xie or Wood's approximation
#'
#' This wrapper function selects the approximation method.
#'
#' @param coef Vector of coefficients (weights).
#' @param f Observed test statistic.
#' @param chi2approx Approximation method: 'bx' or 'wood'
#'
#' @return A data frame with adjusted test and p-value.
#' @noRd
dr.pvalue <- function(coef, f, chi2approx = c("bx", "wood"), ...) {
  method <- match.arg(chi2approx)
  if (method == "bx") {
    bentlerxie.pvalue(coef, f)
    } else {
      wood.pvalue(coef, f, ...)
  }
}

#' Bentler-Xie p-value approximation
#'
#' This function approximates P(coef'X > f), where X ~ chi^2(1), coef is a vector
#' of weights, and f is the observed value.
#'
#' @param coef Vector of weights.
#' @param f Observed statistic.
#'
#' @return A data frame with adjusted p-value and degrees of freedom.
#' @noRd
bentlerxie.pvalue <- function(coef, f) {
  trace1 <- sum(coef)
  trace2 <- sum(coef^2)
  df.adj <- trace1^2 / trace2
  stat.adj <- f *  df.adj / trace1
  bxadjusted <- 1 - pchisq(stat.adj, df.adj)
  ans <- data.frame(test = f, test.adj = stat.adj, df.adj = df.adj,
                    pval.adj = bxadjusted)
  ans
}

#' Wood's p-value approximaation
#'
#' This function implements Wood (1989) approximation for mixture chi-squared.
#'
#' @param coef Vector of weights (eigenvalues).
#' @param f Observed test statistic.
#' @param tol Tolerance for approximation.
#' @param print Logical; print method used.
#'
#' @return A data frame with test result and p-value.
#' @noRd
wood.pvalue <- function (coef, f, tol = 0.0, print = FALSE) {
  if (min(coef) < 0) stop("negative eigenvalues")
  if (length(coef) == 1) {
    pval <- 1 - pchisq(f / coef, 1)
    } else {
      k1 <- sum(coef)
      k2 <- 2 * sum(coef^2)
      k3 <- 8 * sum(coef^3)
      t1 <- 4 * k1 * k2^2 + k3 * (k2 - k1^2)
      t2 <- k1 * k3 - 2 * k2^2
      if ((t2 <= tol) && (tol < t2)) {
        a1 <- k1^2 / k2
        b  <- k1 / k2
        pval <- 1 - pgamma(b * f, a1)
        if (print) print(paste("Satterthwaite-Welsh Approximation =", pval))
        } else if ((t1 <= tol) && (tol < t2)) {
        a1 <- 2 + (k1 / k2)^2
        b  <- (k1 * (k1^2 + k2)) / k2
        pval <- if (f < tol) 1.0 else 1 - pgamma(b / f, a1)
        if (print) print(paste("Inverse gamma approximation =", pval))
        } else if ((t1 > tol) && (t2 > tol)) {
        a1 <- (2 * k1 * (k1 * k3 + k2 * k1^2 - k2^2)) / t1
        b <- t1 / t2
        a2 <- 3 + 2 * k2 * (k2 + k1^2) / t2
        pval <- 1 - pbeta(f / (f + b), a1, a2)
        if (print) print(paste("Three parameter F(Pearson Type VI) approximation =",
                               pval))
        } else {
        pval <- NULL
        if (print) print("Wood's Approximation failed")
        }
      }
  data.frame(test = f, test.adj = NA, df.adj = NA, pval.adj = pval)
  }

#########################################################################
# permutation tests for dimension reduction
#########################################################################

#' Permutation Test for Dimension Reduction
#'
#' This function performs permutation tests for the dimension of the central
#' subspace.  This test compares observed test statistics to their distribution
#' under permutation of predictors (or weights).
#'
#' @param object A fitted \code{dr} object.
#' @param npermute Number of permutations.  Default is 50.
#' @param numdir Maximum number of directions to test.  Default to
#' \code{object$numdir}.
#'
#' @return An object of class \code{dr.permutation.test} with components:
#' \item{summary}{A data frame with test statistics and p-values.}
#' \item{npermute}{The number of permutations used.}
#'
#' @noRd
#' @export
dr.permutation.test <- function(object, npermute = 50, numdir = object$numdir) {
  if (inherits(object, "ire"))
    stop("Permutation tests not implemented for ire")

    permute.weights <- TRUE
    call <- object$call
    call[[1]] <- as.name("dr.compute")
    call$formula <- call$data <- call$subset <- call$na.action <- NULL
    x <- dr.directions(object) # project data on estimated directions
    call$y <- object$y
    weights <- object$weights

    # Number of dimensions to test
    nd <- min(numdir, length(which(abs(object$evalues) > 1.e-8)) - 1)
    nt <- nd + 1

    # Observed test statistics
    obstest<-dr.permutation.test.statistic(object,  nt)
    count <- rep(0, nt)
    val <- rep(0, nt)

    # Main permutation loop
    for (j in 1:npermute) {
      perm <- sample(1:object$cases)
      call$weights <- if (permute.weights) weights[perm] else weights

      # Inner loop for each test (0D vs 1D, 1D vs 2D, ..., (nd-1)D vs ndD)
      for (col in 0:nd) {
        # Permute all but the first col components
        call$x <- if (col==0) x[perm, ] else cbind(x[, 1:col], x[perm, -(1:col)])
        iperm <- eval(call)
        val[col + 1] <- dr.permutation.test.statistic(iperm, col + 1)[col + 1]
      }

      # Update counter if permuted start >= observed
      count[val > obstest] <- count[val > obstest] + 1
    }

    # Compute permutation-based p-values
    pval <- (count) / (npermute + 1)
    ans1 <- data.frame(cbind(obstest, pval))
    dimnames(ans1) <- list(paste(0:(nt - 1), "D vs >= ", 1:nt, "D", sep = ""),
                          c("Stat", "p.value"))
    ans <- list(summary = ans1, npermute = npermute)
    class(ans) <- "dr.permutation.test"
    ans
}

#' @method print dr.permutation.test
#' @export
print.dr.permutation.test <- function(x, digits = max(3, getOption("digits") - 3),
                                      ...) {
    cat("\nPermutation tests\nNumber of permutations:\n")
    print.default(x$npermute)
    cat("\nTest results:\n")
    print(x$summary, digits = digits)
    invisible(x)
 }

#' @method summary dr.permutation.test
#' @export
summary.dr.permutation.test <- function(...) {
  print.dr.permutation.test(...)
}


#########################################################################
# Permutation Test Statistic Methods for Dimension Reduction
#
# These functions compute the test statistics used in permutation tests
# for various dimension reduction methods.
#########################################################################
#' Generic method for computing permutation test statistics.
#'
#' This function computes the test statistics for dimennsion reduction
#' permutation tests.
#'
#' @param object A fitted \code{dr} object.
#' @param numdir Number of directions t otest.
#'
#' @return A numeric vector of test statistics for each hypothesized dimension.
#'
#' @noRd
#' @export
dr.permutation.test.statistic <- function(object, numdir) {
  UseMethod("dr.permutation.test.statistic")
  }

#' Default method for permutation test statistics
#'
#' This function uses the sum of eigenvalues for computing test statistics.
#'
#' @param object A fitted \code{dr} object.
#' @param numdir Number of directions to test.
#'
#' @return A numeric vector of test statistics.
#'
#' @noRd
#' @method dr.permutation.test.statistic default
#' @export
dr.permutation.test.statistic.default <- function(object, numdir) {
  object$cases * rev(cumsum(rev(object$evalues)))[1:numdir]
}

#' Permutation test statistic for method = "phdy"
#'
#' This function reuses the implementation for \code{phd}.
#'
#' @param object A \code{dr} object using the \code{phdy} method.
#' @param numdir Number of directions to test.
#'
#' @noRd
#' @method dr.permutation.test.statistic phdy
#' @exportS3Method
dr.permutation.test.statistic.phdy <- function(object, numdir) {
  dr.permutation.test.statistic.phd(object,numdir)
}

#' Permutation test statistic for method = "phdres"
#'
#' This function reuses the implementation for \code{phd}.
#'
#' @param object A \code{dr} object using the \code{phdres} method.
#' @param numdir Number of directions to test.
#'
#' @noRd
#' @method dr.permutation.test.statistic phdres
#' @exportS3Method
dr.permutation.test.statistic.phdres <- function(object, numdir) {
  dr.permutation.test.statistic.phd(object,numdir)
}

#' Permutation test statistic for method = "phd"
#'
#' This function computes test statistics based on squared eigenvalues
#' and response variance.
#'
#' @param object A \code{dr} object using the \code{phd} method.
#' @param numdir Number of directions to test.
#'
#' @noRd
#' @method dr.permutation.test.statistic phd
#' @exportS3Method
dr.permutation.test.statistic.phd <- function(object,numdir) {
  (0.5 * object$cases * rev(cumsum(rev(object$evalues^2))) /
     var(dr.y(object)))[1:numdir]
}

#####################################################################
# dr.slices returns non-overlapping slices based on y
#####################################################################

#' Slice the response variable for dimension reduction
#'
#' Divides a univariate or multivariate response \code{y} into non-overlapping
#' slices for use in sufficient dimension reduction methods like SIR or pHd.
#'
#' @param y A numeric vector or matrix representing the response.
#' @param nslices Either the total number of slices (if \code{y} is univariate),
#'   or a vector giving the number of slices per dimension (if \code{y} is
#'   multivariate).
#'
#' @return A list with:
#'   \item{slice.indicator}{An integer vector indicating slice membership.}
#'   \item{nslices}{Total number of slices.}
#'   \item{slice.sizes}{A vector of counts per slice.}
#'
#' @noRd
#' @export
dr.slices <- function(y, nslices) {
  # Slice a univeriate response into h levels
  dr.slice.1d <- function(y, h) {
    z <- unique(y)
    if (length(z) > h) dr.slice2(y, h) else dr.slice1(y, sort(z))
    }

  # Exact slicing by unique values
  dr.slice1 <- function(y, u) {
    z <- sizes <- 0
    for (j in 1:length(u)) {
      temp <- which(y == u[j])
      z[temp] <- j
      sizes[j] <- length(temp)
      }
    list(slice.indicator = z, nslices = length(u), slice.sizes = sizes)
  }

  # Approximate equal-size slicing
  dr.slice2 <- function(y, h) {
    myfind <- function(x, cty) {
      ans <- which(x <= cty)
      if (length(ans) == 0) length(cty) else ans[1]
      }
    or <- order(y)
    cty <- cumsum(table(y))
    names(cty) <- NULL # drop class names
    n <- length(y)
    m <- floor(n / h) # nominal number of obs per slice
    sp <- end <- 0
    j <- 0            # slice counter will end up <= h
    ans <- rep(1, n)

    # find slice boundaries: all slices have at least 2 obs
    while(end < n - 2) {
      end <- end + m
      j <- j + 1
      sp[j] <- myfind(end, cty)
      end <- cty[sp[j]]
      }
    sp[j] <- length(cty)

    for (j in 2:length(sp)) {
      firstobs <- cty[sp[j - 1]] + 1
      lastobs <- cty[sp[j]]
      ans[or[firstobs:lastobs]] <- j
      }

    list(slice.indicator = ans, nslices = length(sp),
         slice.sizes = c(cty[sp[1]], diff(cty[sp])))
  }

  # Handle multivariate y
  p <- if (is.matrix(y)) dim(y)[2] else 1
  h <- if (length(nslices) == p) nslices else rep(ceiling(nslices^(1 / p)), p)
  a <- dr.slice.1d(if(is.matrix(y)) y[, 1] else y, h[1])

  if (p > 1) {
    for (col in 2:p) {
      ns <- 0
      for (j in unique(a$slice.indicator)) {
        b <- dr.slice.1d(y[a$slice.indicator == j, col], h[col])
        a$slice.indicator[a$slice.indicator == j] <-
          a$slice.indicator[a$slice.indicator == j] + 10^(col - 1) * b$slice.indicator
        ns <- ns + b$nslices
        }
      a$nslices <- ns
    }

    # Recode unique values to sequential slice numbers
    v <- unique(a$slice.indicator)
    L <- slice.indicator <- NULL
    for (i in 1:length(v)) {
      sel <- a$slice.indicator == v[i]
      slice.indicator[sel] <- i
      L <- c(L, length(a$slice.indicator[sel]))
      }
    a$slice.indicator <- slice.indicator
    a$slice.sizes <- L
    }
  a
}

#' Slice univariate y as in Arc (original S-PLUS/Arc way)
#'
#' Matches the slicing procedure used in the Arc software by Cook et al.
#'
#' @param y A univariate numeric vector (response).
#' @param nslices Number of desired slices.
#'
#' @return A list with slice indicators and sizes.
#'
#' @noRd
#' @export
dr.slices.arc <- function(y, nslices) {
  if(is.matrix(y))  stop("dr.slices.arc is used for univariate y only.
                         Use dr.slices")
  h <- nslices
  or <- order(y)
  n <- length(y)
  m <- floor(n / h)
  r <- n - m * h
  start <- sp <- ans <- 0
  j <- 1

  # Determine slice boundaries, handling tied values carefully
  while((start + m) < n) {
    if (r==0) {
    start<-start
  } else {
    start <- start + 1
    r <- r - 1
  }
  while (y[or][start + m] == y[or][start + m + 1]) {
    start <- start + 1
  }
  sp[j] <- start + m
  start <- sp[j]
  j <- j + 1
  }

  # Ensure the last slice has at least 2 observations
  if (sp[j - 1] == n - 1) j <- j - 1
  sp[j] <- n

  # Assign slice labels
  ans[or[1:sp[1]]] <- 1
  for (k in 2:j) {
    ans[or[(sp[k - 1] + 1):sp[k]]] <- k
    }

  list(slice.indicator = ans, nslices = j, slice.sizes = c(sp[1], diff(sp)))
}

#####################################################################
#     Misc. Auxillary functions: cosine of the angle between two
#     vectors and between a vector and a subspace.
#####################################################################

#' Cosine of angle between vectors and a subspace
#'
#' This function computes the cosine of the angle(s) between one or more
#' vectors and the subspace spanned by the columns of a matrix. This is
#' useful in dimension reduction diagnostics to evaluate alignment of directions.
#'
#' @param mat A numeric matrix whose columns define the subspace.
#' @param vecs A numeric vector or matrix. If a matrix, angles are computed for each column.
#'
#' @return A numeric vector or matrix of squared cosine values.
#'
#' @noRd
#' @export
cosangle <- function(mat, vecs) {
  ans <-NULL
  # If vecs is a single vector
  if (!is.matrix(vecs)) {
    ans <- cosangle1(mat, vecs)
    } else {
      # Compute angle for each column vector of vecs
      for (i in 1:dim(vecs)[2]) {
        ans <- rbind(ans, cosangle1(mat, vecs[, i]))
      }
      dimnames(ans) <- list(colnames(vecs), NULL)
      }
  ans
}

#' Internal function for angle computation
#'
#' This function computes the squared cosine of the angle between a vector
#' and a subspace.
#'
#' @param mat A matrix defining the subspace.
#' @param vec A numeric vector.
#'
#' @return A numeric vector of cumulative squared cosines.
#' @noRd
cosangle1 <- function(mat, vec) {
  # Project vec onto orthogonal basis of mat and compute squared projection
  # lengths
  ans <- cumsum((t(qr.Q(qr(mat))) %*% scale(vec) / sqrt(length(vec) - 1))^2)

  # Fix for compatibility with S-PLUS
  if (version$language == "R") ans else t(ans)
}

#####################################################################
# R Functions for reweighting for elliptical symmetry
# modified from reweight.lsp for Arc
# Sanford Weisberg, sandy@stat.umn.edu
# March, 2001, rev June 2004

# Here is an outline of the function:
#   1.  Estimates of the mean m and covariance matrix S are obtained.  The
#       function cov.rob in the MASS package is used for this purpose.
#   2.  The matrix X is replaced by Z = (X - 1t(m))S^{-1/2}.  If the columns
#       of X were from a multivariate normal N(m,S), then the rows of Z are
#       multivariate normal N(0, I).
#   3.  Randomly sample from N(0, sigma*I), where sigma is a tuning
#       parameter.  Find the row of Z closest to the random data, and increase
#       its weight by 1
#   4.  Return the weights divided by nsamples and multiplied by n.
#
#     dr.weights
#
#####################################################################
#' Compute Reweighting for Elliptical Symmetry
#'
#' This function computes reweighted observations to approximate elliptical
#' symmetry using the method described in Cook and Weisberg (1999). Based on
#' random sampling in standardized space.
#'
#' @param formula A model formula.
#' @param data A data frame containing variables in the model.
#' @param subset Optional subset of the data.
#' @param na.action A function which indicates what should happen when the data
#'    contain NAs.
#' @param sigma A tuning parameter controlling dispersion of normal samples
#'    (default 1).
#' @param nsamples Number of Monte Carlo samples. Defaults to \code{10 * n}.
#' @param ... Additional arguments passed to \code{cov.rob}.
#'
#' @return A numeric vector of weights.
#' @noRd
#' @export
dr.weights <- function (formula, data = list(), subset, na.action=na.fail,
            sigma = 1, nsamples = NULL, ...) {
    # Build model frame and extract matrix
    mf1 <- match.call(expand.dots = FALSE)
    mf1$... <- NULL
    mf1$covmethod <- mf1$nsamples <- NULL
    mf1[[1]] <- as.name("model.frame")
    mf <- eval(mf1, sys.frame(sys.parent()))
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    int <- match("(Intercept)", dimnames(x)[[2]], nomatch = 0)
    if (int > 0) x <- x[, -int, drop = FALSE] # remove intercept

    # Estimate robust center and covariance
    ans <- cov.rob(x, ...)
    m <- ans$center
    s <- svd(ans$cov)

    # Standardize design matrix
    z <- sweep(x, 2, m) %*% s$u %*% diag(1 / sqrt(s$d))
    n <- dim(z)[1]
    p <- dim(z)[2]
    ns <- if (is.null(nsamples)) 10 * n else nsamples

    dist <- wts <- rep(0, n)
    for (i in 1:ns) {
      point <- rnorm(p) * sigma
      dist <- apply(z, 1, function(x, point) {sum((point - x)^2)}, point)
      #Dist to each point
      sel <- dist == min(dist)               # Find closest point(s)
      wts[sel] <- wts[sel] + 1 / length(wts[sel])
    }

    w <- n * wts / ns
    if (missing(subset)) return(w)
    if (is.null(subset)) return(w) else {
      # find dimension of mf without subset specified
      mf1$subset <- NULL
      w1 <- rep(NA, length = dim(eval(mf1))[1])
      w1[subset] <- w
      return(w1)
    }
}

#' Drop1 Method for Dimension Reduction
#'
#' This function implements drop-one-variable testing based on coordinate
#' hypotheses.
#'
#' @param object A \code{dr} object.
#' @param scope Terms not to be dropped.
#' @param update Whether to update the model or return a table only.
#' @param test Type of test: "general" or "marginal".
#' @param trace If > 0, print progress.
#' @param ... Additional arguments passed to \code{dr.coordinate.test}.
#'
#' @return An updated \code{dr} object or a data frame of test results.
#'
#' @noRd
#' @export
drop1.dr <- function (object, scope = NULL, update = TRUE, test = "general",
            trace=1, ...) {
    keep <- if (is.null(scope))
      NULL
    else attr(terms(update(object$terms, scope)), "term.labels")
    all <- attr(object$terms, "term.labels")
    candidates <- setdiff(all, keep)
    if (length(candidates) == 0) stop("Error---nothing to drop")

    ans <- NULL
    for (label in candidates) {
      ans <- rbind(ans, dr.coordinate.test(object, as.formula(paste("~.-",
                                label, sep = "")), ...))
    }
    row.names(ans) <- paste("-", candidates)

    ncols <- ncol(ans)
    or <- order(-ans[, if (test == "general")
      ncols
      else (ncols - 1)])

    form <- formula(object)
    attributes(form) <- NULL
    fout <- deparse(form, width.cutoff = 50)
    if(trace > 0) {
      for (f in fout) cat("\n", f)
      cat("\n")
      print(ans[or, ]) }
    if (is.null(object$stop)) {
      object$stop <- 0
    }
    stopp <- if (ans[or[1], if (test == "general")
      ncols
      else (ncols - 1)] < object$stop)
      TRUE
    else FALSE
    if (stopp == TRUE) {
      if(trace > 0) cat("\nStopping Criterion Met\n")
      object$stop <- TRUE
      object
    }
    else if (update == TRUE) {
      update(object, as.formula(paste("~.", row.names(ans)[or][1],
                                      sep = "")))
    }
    else invisible(ans)
  }

#' Stepwise Dimension Reduction
#'
#' Applies recursive drop-one-variable testing until stopping criterion is met.
#'
#' @param object A \code{dr} object.
#' @param scope Optional scope for variables to retain.
#' @param d Currently unused.
#' @param minsize Minimum number of variables to retain (default 2).
#' @param stop Threshold test statistic for stopping.
#' @param trace If > 0, print progress.
#' @param ... Additional arguments passed to \code{drop1.dr}.
#'
#' @return A reduced \code{dr} object.
#' @noRd
#' @export
dr.step <- function (object, scope = NULL, d = NULL, minsize = 2, stop = 0,
            trace=1,...) {
    if (is.null(object$stop)) {
      object$stop <- stop
    }
    if (object$stop == TRUE)
      object
    else {
      minsize <- max(2, minsize)
      keep <- if (is.null(scope))
        NULL
      else attr(terms(update(object$terms, scope)), "term.labels")
      all <- attr(object$terms, "term.labels")
      if (length(keep) >= length(all)) {
        if(trace > 0) cat("\nNo more variables to remove\n")
        object
      }
      else if (length(all) <= minsize) {
        if (trace > 0) cat("\nMinimum size reached\n")
        object$numdir <- minsize
        object
      }
      else {
        if (dim(object$x)[2] <= minsize) {
          if (trace > 0) cat("\nMinimum size reached\n")
          object}
        else {
          obj1 <- drop1(object, scope = scope, d = d, trace=trace, ...)
          dr.step(obj1, scope = scope, d = d, stop = stop,
                  trace=trace, ...)
        }
      }
    }
  }

#####################################################################
# Utility Functions (R-compatible replacements for S-Plus)
#####################################################################
#' Compatibility Utilities for S-Plus in R
#' This function provides utility functions (`is.empty.model`, `NROW`,
#' `NCOL`, `colnames`, `getOption`) to ensure compatibility between R
#' and legacy S-Plus code.  These definitions are only added when the langugage
#' is not R (i.e., running in S-Plus).
#' @noRd
if (version$language != "R") {
  "is.empty.model" <- function (x)
  {
    tt <- terms(x)
    (length(attr(tt, "factors")) == 0) & (attr(tt, "intercept") == 0)
  }
  "NROW" <-
    function(x) if(is.array(x)||is.data.frame(x)) nrow(x) else length(x)
  "NCOL" <-
    function(x) if(is.array(x)||is.data.frame(x)) ncol(x) else as.integer(1)
  "colnames" <-
    function(x, do.NULL = TRUE, prefix = "col")
    {
      dn <- dimnames(x)
      if(!is.null(dn[[2]]))
        dn[[2]]
      else {
        if(do.NULL) NULL else paste(prefix, seq(length=NCOL(x)), sep="")
      }
    }

  "getOption" <- function(x) options(x)[[1]]
  # end of special functions
}

#' Fit Inverse Regression Estimator (IRE)
#'
#' This function estimates the central subspace using inverse regression
#' methodology as described in Cook and Ni (2005).  It implements a minimum
#' discrepancy approach based on slicing the response and computing slice
#' means in the inverse regression space.
#'
#' @param object An object of class \code{dr}.
#' @param numdir Integer. Number of directions (dimension) to estimate.
#'     Default is 4.
#' @param nslices Integer. Number of slices to use in the inverse regression.
#'     If \code{NULL}, a default value based on the number of predictors is used.
#' @param slice.function A function for computing slices. Default is \code{dr.slices}.
#' @param tests Logical. If \code{TRUE}, performs hypothesis testing for each dimension
#'     up to \code{numdir}.
#' @param ... Additional arguments passed to lower-level functions like \code{dr.test}.
#'
#' @return An updated object of class \code{dr} containing fitted IRE results,
#'     including estimated directions, slice means, and test statistics.
#' @noRd
dr.fit.ire <- function(object, numdir = 4, nslices = NULL,
                       slice.function = dr.slices, tests = TRUE, ...) {

  # Extract the centered and scaled predictor matrix
  z <- dr.z(object)

  # Determine number of slices (default: max(8, p+3))
  h <- if (!is.null(nslices)) nslices else max(8, NCOL(z) + 3)

  # Slice the response variable using the specified slicing function
  slices <- slice.function(dr.y(object), h)
  object$slice.info <- slices # Store slice information in object

  # Compute slice proportions
  f <- slices$slice.sizes / sum(slices$slice.sizes)

  n <- object$cases # Total number of cases
  weights <- object$weights # Observation weights
  p <- dim(z)[2] # Number of predictors

  # Ensure number of directions does not exceed p-1
  numdir <- min(numdir, p - 1)
  h <- slices$nslices # Actual number of slices (may differ slightly from requested)

  # Initialize matrix of slice means in predictor space
  xi <- matrix(0, nrow = p, ncol = slices$nslices)
  for (j in 1:h) {
    sel <- slices$slice.indicator == j # Indicator for observations in slice j
    xi[, j] <- apply(z[sel, ], 2, function(a, w) sum(a * w) / sum(w), weights[sel])
  }

  # Compute eigenvectors of the weighted covariance of slice means
  object$sir.raw.evectors <- eigen(xi %*% diag(f) %*% t(xi))$vectors

  # Save slice means to object
  object$slice.means <- xi

  # Create orthonormal basis for contrast space (Helmert matrix)
  An <- qr.Q(qr(contr.helmert(slices$nslices)))

  # Compute zeta matrix (Eq. (1), Cook & Ni 2005) in Q-transformed coordinates
  object$zeta <- xi %*% diag(f) %*% An
  rownames(object$zeta) <- paste("Q",1:dim(z)[2], sep = "")

  ans <- NULL # Initialize result list

  if (tests == TRUE) {
    # Perform independence test for null dimension
    object$indep.test <- dr.test(object, numdir = 0, ...)

    # Compute Gamma_zeta (Gz) matrix once for efficiency
    Gz <- Gzcomp(object, numdir)

    # Perform tests for dimensions d=1 to numdir
    for (d in 1:numdir) {
      ans[[d]] <- dr.test(object, numdir = d, Gz, ...)
      colnames(ans[[d]]$B) <- paste("Dir", 1:d, sep = "")
    }

    # Compute final Gz for highest fitted dimension using the final estimated B
    object$Gz <- Gzcomp(object, d, span = ans[[numdir]]$B)
  }

  # Return updated dr object with resullts
  aa <- c(object, list(result = ans, numdir = numdir, n = n, f = f))
  class(aa) <- class(object)
  return(aa)
}

#' Large-Sample Test for IRE
#'
#' This function performs marginal dimension hypothesis tests for IRE directions.
#'
#' @param object An object returned by \code{dr.fit.ire}.
#' @param numdir Number of directions to test.
#' @param Gz Precomputed Gz matrix.
#' @param steps Number of re-estimation steps.
#' @param ... Additional arguments.
#'
#' @return A list with direction matrix \code{B} and test summary.
#' @noRd
dr.test.ire <- function(object, numdir, Gz = Gzcomp(object,numdir),
                        steps = 1, ...) {

  # Initial estimate of B using raw eigenvectors from SIR step
  ans <- dr.iteration(object, Gz, d = numdir,
                      B = object$sir.raw.evectors[, 1:numdir, drop = FALSE], ...)

  # Optionally re-estimate B using iterative updates for stability/refinement
  if (steps > 0 & numdir > 0) {
    for (st in 1:steps) {
      Gz <- Gzcomp(object, numdir, ans$B[, 1:numdir]) # Update Gz using current span
      ans <- dr.iteration(object, Gz, d = numdir,
                          B = ans$B[, 1:numdir, drop = FALSE],...)
    }
  }

  # Sequential reordering of directions according to importance (Cook & Ni, p. 414)
  if (numdir > 1) {
    ans0 <- dr.iteration(object, Gz, d = 1, T = ans$B) # First direction
    sumry <- ans0$summary # Save test summary
    B0 <- matrix(ans0$B / sqrt(sum((ans0$B)^2)), ncol = 1) # Normalize
    C <- ans$B # Copy of current B for modification

    # Loop over remaining directions
    for (d in 2:numdir) {
      common <- B0[, d - 1] # Last added direction

      # Orthogonalize remaining columns of C against the current B0
      for (j in 1:dim(C)[2]) {
        C[, j] <- C[, j] - sum(common * (C[, j])) * B0[, d - 1]
        }

      # QR decomposition to re-orthogonalize remaining directions
      C <- qr.Q(qr(C))[, -dim(C)[2], drop = FALSE]

      # Find the next best direction orthoogonal to previous ones
      ans0 <- dr.iteration(object, Gz, d = 1, T = C)

      # Normalize and append to B0
      B0 <- cbind(B0, ans0$B / sqrt(sum((ans0$B)^2)))

      # Store the sequential test result
      sumry <- rbind(sumry, dr.iteration(object, Gz, d = d, T = B0)$summary)
    }

    # Ensure directions have unit norm and first component positive
    ans$B <- apply(B0, 2, function(x) {
      b <- x / sqrt(sum(x^2))
      if (b[1] < 0) - b else b
      })

    ans$sequential <- sumry # Store test summaries
  }

  # Add column and row names for readability
  if (numdir > 0) colnames(ans$B) <- paste("Dir", 1:numdir, sep = "")
  if (numdir > 1) rownames(ans$sequential) <- paste(1:numdir)

  return(ans)
}

#' @export
Gzcomp <- function(object, numdir, span) {
  UseMethod("Gzcomp")
}

#' Compute Gz Matrix for Inverse Regression Estimator (IRE)
#'
#' This function computes the inner product matrix \eqn{\Gamma_z} (or its
#' Cholesky factor) for the inverse regression estimator (IRE), following Cook
#' and Ni (2005).
#'
#' @param object A fitted \code{dr} object with method \code{"ire"}.
#' @param numdir Integer. Number of directions to project onto.
#' @param span Optional. A projection matrix (e.g., from fitted directions); if NULL, uses identity.
#'
#' @return A Cholesky factor of the \eqn{\Gamma_z} matrix.
#' @noRd
Gzcomp.ire <- function(object, numdir, span = NULL) {
  slices <- object$slice.info # slicig structure
  An <- qr.Q(qr(contr.helmert(slices$nslices))) # Orthonormal contrast matrix
  n <- object$cases # Sample size
  weights <- object$weights # Weights for each observation
  z <- dr.z(object)  # Standardized covariate matrix (centered & weighted)
  p <- dim(object$slice.means)[1] # Number of predictors
  h <- slices$nslices # Number of slices
  f <- slices$slice.sizes/sum(slices$slice.sizes) # Slice proportions

  # Use provided projection span or default to identify matrix
  span <- if (!is.null(span)) span else diag(rep(1, p))

  # Project the slice means to the subspace spanned by span, or use 0 if numdir = 0
  xi <- if (numdir > 0) {
    qr.fitted(qr(span), object$slice.means) # Projected slice means
    } else {
    matrix(0, nrow = p, ncol = h)
    }

  # Initialize the Gram matrix Gz, which will store vec(z_i * epsilon_i') rows
  Gz <- array(0, c(p * h, p * h))
  Gmat <- NULL

  # Construct each vec(z_i * epsilon_i') row for all i in 1, ..., n
  for (i in 1:n) {
    # Compute epsilon_i = e_j - f - z_i %*% xi %*% diag(f)
    epsilon_i <- -f - z[i, , drop = FALSE] %*% xi %*% diag(f)
    epsilon_i[slices$slice.indicator[i]] <- 1 + epsilon_i[slices$slice.indicator[i]]

    # Compute vec(z_i * epsilon_i') as a row vector
    vec <- as.vector(z[i, ] %*% epsilon_i)

    # Stack vecs row-wise into Gmat
    Gmat <- rbind(Gmat, vec)
  }

  # Estimate Gamma_z: coompute sample covariance of Gmat, scale by (n-1)/n
  # Then transform it with kronecker products involving the An matrix
  # Only the Cholesky decomposition is computed sinnce it is sufficient for later use
  Gz <- chol(kronecker(t(An), diag(rep(1, p))) %*% (((n - 1) / n) * cov(Gmat)) %*%
               kronecker(An, diag(rep(1, p))))

  return(Gz)
}

#####################################################################
## Iterative Estimation for IRE (Inverse Regression Estimator)
##
## Purpose:
##   This function performs iterative estimation of the central
##   subspace using the algorithm proposed by Cook & Ni (2005).
##
## Arguments:
##   - B: A p × d matrix whose columns provide an initial estimate
##        of the basis vectors for the central subspace.
##        Defaults to the first d columns of the identity matrix I_p.
##
##   - R: A restriction matrix (optional). If specified, the estimated
##        spanning directions are returned as R %*% B, rather than B itself.
##
## Return:
##   - The returned matrix B is the product R %*% B, if R is specified.
##
## Algorithm Notes:
##   - fn:        The objective function to minimize (Eq. (5), Cook & Ni, 2005).
##
##   - updateC:   Step 2 of the algorithm (p. 414): updates the intermediate
##                coefficient matrix C based on the current estimate of B.
##
##   - updateB:   Step 3 of the algorithm (p. 414): updates B by solving an
##                optimization problem using the most recent C.
##
##   - PBk:       In the paper, PBk denotes the projection onto the orthogonal
##                complement of the space spanned by QBk. In the code, we compute
##                the QR decomposition of the **projection**, not the complement.
##                Hence, we use `qr.resid()` instead of `qr.fitted()` to compute
##                residuals from the projection.
#####################################################################

#' @export
dr.iteration <- function(object, Gz, d = 2, B, T, eps, itmax, verbose) {
  UseMethod("dr.iteration")
}

#' Internal Iteration Function for IRE Optimization
#'
#' Performs iterative optimization of the inverse regression estimator (IRE)
#' using the algorithm described in Cook & Ni (2005). It alternates between
#' updating the basis matrix \code{B} and coefficient matrix \code{C} to minimize
#' the objective function.
#'
#' @param object A \code{dr} object containing required components such as
#'     \code{zeta}.
#' @param Gz A precomputed Cholesky decomposition of the \eqn{\Gamma_\zeta}
#'    matrix.
#' @param d Integer. Target number of sufficient directions.
#' @param B Optional. Initial orthonormal basis matrix for the subspace (p × d).
#' @param T Optional. Transformation matrix (typically identity or restriction).
#' @param eps Convergence tolerance. Default is \code{1e-6}.
#' @param itmax Maximum number of iterations. Default is 200.
#' @param verbose Logical. If \code{TRUE}, prints progress at each iteration.
#'
#' @return A list with:
#'   \item{\code{B}}{Estimated basis matrix of the central subspace.}
#'   \item{\code{summary}}{A data frame with test statistic, degrees of
#'       freedom, p-value, and number of iterations.}
#' @noRd
#' @exportS3Method dr.iteration ire
dr.iteration.ire <- function(object, Gz, d = 2, B = NULL, T = NULL, eps = 1.e-6,
                             itmax = 200, verbose = FALSE) {
  n <- object$cases
  zeta <- object$zeta
  p <- dim(zeta)[1]
  h1 <- dim(zeta)[2]  # h1 is h-1 for ire, sum(h_i - 1) for pire

  # If d=0, return the test statistic (no iteration needed)
  if (d == 0) {
    err <- n * sum(forwardsolve(t(Gz), as.vector(zeta))^2)
    return(data.frame(Test = err, df = (p - d) * (h1 - d),
               p.value = pchisq(err, (p - d) * (h1 - d), lower.tail = FALSE),
               iter = 0))
    } else {
      # Initialize transformation matrix T if not provided
      T <- if(is.null(T)) diag(rep(1, p)) else T

      # Initialize starting B matrix (p x d) if not provided
      B <- if(is.null(B)) diag(rep(1, ncol(T)))[, 1:d, drop = FALSE] else B

      # Define objective function: squared norm of transformed residuals
      fn <- function(B, C) {
      n * sum(forwardsolve(t(Gz), as.vector(zeta) - as.vector(T %*% B %*% C))^2 )
      }

      # Step 2 of the algorithm (update C given B)
      updateC <- function() {
        design <- forwardsolve(t(Gz), kronecker(diag(rep(1, h1)), T %*% B))
        response <- forwardsolve(t(Gz), as.vector(zeta))
        matrix(qr.coef(qr(design), response), nrow = d)
      }

      # Step 3 of the algorithm (update B given C)
      updateB <- function() {
        for (k in 1:d) {
          # Compute partial residuals excluding direction k
          alphak <- as.vector(zeta - T %*% B[, -k, drop = FALSE] %*% C[-k, ])

          # Project onto orthogonal complement of current B excluding column k
          PBk <- qr(B[, -k])

          # Regress residuals on current column of C to update B[, k]
          bk <- qr.coef(
            qr(forwardsolve(t(Gz), t(qr.resid(PBk, t(kronecker(C[k, ], T)))))),
            forwardsolve(t(Gz), as.vector(alphak)))

          # If NA (due to singularities), replace with 0
          bk[is.na(bk)] <- 0

          # Orthogonalize relative to previous columns
          bk <- qr.resid(PBk, bk)

          # Normalize to unit length
          B[, k] <- bk / sqrt(sum(bk^2))
        }
        B
      }

      # Initialize C and error
      C <- updateC()
      err <- fn(B, C)
      iter <- 0

      # Main optimization loop
      repeat{
        iter <- iter + 1
        B <- updateB()
        C <- updateC()
        errold <- err
        err <- fn(B, C)

        # Optional progress output
        if (verbose == TRUE)
          print(paste("Iter =", iter, "Fn =", err), quote = FALSE)

        # Check convergence or iteration limit
        if (abs(err - errold) / errold < eps || iter > itmax ) break
      }

      # Return final result: transformed B and summary statistics
      B <- T %*% B
      rownames(B) <- rownames(zeta)

      list(B = B, summary = data.frame(
        Test = err, df = (p - d) * (h1 - d),
        p.value = pchisq(err, (p - d) * (h1 - d), lower.tail = FALSE),
        iter = iter))
    }
}

#' Coordinate Hypothesis Test for IRE
#'
#' Performs a coordinate hypothesis test for the Inverse Regression Estimator (IRE),
#' based on equations (12) and (14) in Cook & Ni (2005). The test evaluates whether
#' a specified linear combination (subspace) of the predictors contributes to the
#' model.
#'
#' @param object A \code{dr} object fitted using the IRE method.
#' @param hypothesis Either a matrix or a \code{formula} specifying the
#'     coordinate hypothesis to test.
#' @param d Integer (optional). Dimension of the sufficient subspace under the
#'     alternative hypothesis.  If \code{NULL}, the full model is tested.
#' @param ... Additional arguments passed to \code{dr.iteration} or
#'     \code{dr.joint.test}.
#'
#' @return A data frame containing:
#' \item{\code{Test}}{Chi-squared test statistic.}
#' \item{\code{df}}{Degrees of freedom.}
#' \item{\code{p.value}}{Associated p-value.}
#'
#' @exportS3Method dr.coordinate.test ire
#' @noRd
dr.coordinate.test.ire <- function(object, hypothesis, d = NULL, ...) {
  # Convert hypothesis to matrix (if formula, extract design matrix basis)
  gamma <- if (inherits(hypothesis, "formula")) {
    coord.hyp.basis(object, hypothesis)
  } else {
    as.matrix(hypothesis)
  }

  # Project the hypothesis into the Q-coordinate system
  gamma <- dr.R(object) %*% gamma  # Rotate too orthonomal coordinate system

  p <- dim(object$x)[2] # Number of predictors
  r <- p - dim(gamma)[2] # Rank (dimension) of the null hypothesis
  maxdir <- length(object$result) # Maximum dimension estimated in the model

  if (is.null(d)) {
    # Case: Full model hypothesis test (Equation 12 in Cook & Ni, p. 415)
    h1 <- dim(object$zeta)[2] # h-1 for ire
    H <- qr.Q(qr(gamma), complete = TRUE)[, -(1:(p - r)), drop=FALSE]
    n <- object$cases
    Gz <- object$Gz
    zeta <- object$zeta

    # Transform H and zeta for test statistic computation
    m1 <- Gz %*% kronecker(diag(rep(1, h1)), H)
    m1 <- chol(t(m1) %*% m1)
    T_H <- n * sum (forwardsolve(t(m1), as.vector(t(H) %*% zeta))^2)

    df <- r * h1 # Degrees of freedom
    z <- data.frame(Test = T_H, df = df,
                    p.value = pchisq(T_H, df, lower.tail = FALSE))
    z
    } else {
      # Case: Coordinate test for dimension d (Equation 14 in Cook & Ni, p. 415)

      # Null modell test statistic (F0)
      F0 <-if(maxdir >= d) {
        object$result[[d]]$summary$Test
        } else {
          dr.iteration(object, object$Gz, d = d, ...)$summary$Test
        }

      # Alternative model test statistic (F1) with hypothesis restriction
      F1 <- dr.joint.test(object, hypothesis, d = d, ...)$summary$Test


      data.frame(Test = F1 - F0, df = r * d,
                 p.value = pchisq(F1 - F0, r * d, lower.tail = FALSE))
    }
}

#' Joint Hypothesis Test for IRE
#'
#' Performs a joint hypothesis test in the Inverse Regression Estimator (IRE)
#' framework, based on the unnumbered equation in the second column of p. 415
#' in Cook & Ni (2005). This test assesses whether the sufficient directions lie
#' in the subspace spanned by a hypothesized transformation matrix.
#'
#' @param object A \code{dr} object fitted using the IRE method.
#' @param hypothesis A matrix or \code{formula} specifying the null hypothesis
#'    subspace.
#' @param d Integer (optional). Number of sufficient directions under the
#'    alternative hypothesis.
#' @param ... Additional arguments passed to \code{dr.iteration}.
#'
#' @return Either the result of \code{\link{dr.coordinate.test}} (if
#'    \code{d = NULL}), or a result from \code{\link{dr.iteration}} using the
#'    hypothesized transformation matrix.
#'
#' @exportS3Method dr.joint.test ire
#' @noRd
dr.joint.test.ire <- function(object, hypothesis, d = NULL, ...) {
  # If no specific dimension is given, defer to coordinate test
  if (is.null(d)) {
    dr.coordinate.test(object, hypothesis, ...)
    } else {
      # Convert hypothesis to matrix
      gamma <- if (inherits(hypothesis, "formula")) {
        coord.hyp.basis(object, hypothesis)
      } else {
        as.matrix(hypothesis)
      }

      # Project the hypothesis into the rotated Q-coordinates used by IRE
      gamma <- dr.R(object) %*% gamma

      # Run IRE estimation constrained to lie in the span of gamma
      dr.iteration(object, object$Gz, d = d, T = gamma)
    }
}

#' Print Method for IRE Objects
#'
#' Displays summary output from an IRE model object, including the fitted call
#' and large-sample marginal dimension tests for each estimated direction.
#'
#' @param x An object of class \code{"ire"}, typically produced by a function
#'     such as \code{dr.fit.ire}.
#' @param width Integer. Maximum line width when printing the call. Default is 50.
#' @param ... Further arguments passed to other print methods (ignored here).
#'
#' @method print ire
#' @exportS3Method
#' @noRd
print.ire <- function(x, width = 50, ...) {
  # Print the call used to fit the model, wrapping lines at 'width' characters
  fout <- deparse(x$call, width.cutoff = width)
  for (f in fout) cat("\n", f)
  cat("\n")

  # Number of sufficient directions estimated
  numdir <- length(x$result)

  # Initialize tests with the independence test (for 0 directions)
  tests <- x$indep.test

  # Append test results for each direction
  for (d in 1:numdir) {
    tests <- rbind(tests, x$result[[d]]$summary)
  }

  # Label rows as comparisons of 0:D, 1:D, ..., D:D
  rownames(tests) <- paste(0:numdir, "D vs", " > ", 0:numdir, "D", sep = "")

  # Display the marginal dimension tests
  cat("Large-sample Marginal Dimension Tests:\n")
  print(tests)
  cat("\n")

  # Return the object invisibly (stadard for print methods)
  invisible(x)
}

#' Summary Method for IRE Objects
#'
#' Summarizes the results from an inverse regression estimator (IRE) fit.
#' Includes information about the call, slicing, weights, estimated directions,
#' and large-sample dimension tests.
#'
#' @param object An object of class \code{"ire"} returned by \code{dr.fit.ire}.
#' @param ... Currently ignored. Included for S3 method consistency.
#'
#' @return An object of class \code{"summary.ire"} containing the fitted call,
#' method, slice information, weights, estimated directions, and test results.
#' @exportS3Method
#' @noRd
summary.ire <- function (object, ...) {
  # Start with storing the original call
  ans <- object[c("call")]

  # Extract estimated directions and compute how many were fitted
  result <- object$result
  numdir <- length(result)

  # Combine independence test (0D) and tests for each fitted direction
  tests <- object$indep.test
  for (d in 1:numdir) {
  tests <- rbind(tests, result[[d]]$summary)
  }

  # Label rows as comparisons of 0D, 1D, ..., up to numdir D
  rownames(tests) <- paste(0:numdir, "D vs", " > ", 0:numdir, "D", sep = "")

  # Add additional model components to the summary output
  ans$method <- object$method
  ans$nslices <- object$slice.info$nslices
  ans$sizes <- object$slice.info$slice.sizes
  ans$weights <- dr.wts(object)

  # Add result list with direction matrices replaced by actual bases
  ans$result <- object$result
  for (j in 1:length(ans$result)) {
    ans$result[[j]]$B <- dr.basis(object,j)
    }

  # Add number of observations and test summary
  ans$n <- object$cases
  ans$test <- tests

  # Set the summary class for future method dispatch
  class(ans) <- "summary.ire"
  ans
}

#' Print Method for Summary of IRE Objects
#'
#' Nicely formats and prints the summary information for an inverse regression
#' estimator (IRE) object, including call, method details, slice sizes,
#' large-sample dimension tests, and estimated directions.
#'
#' @param x An object of class \code{"summary.ire"}, typically from
#'     \code{summary.ire}.
#' @param digits Number of significant digits to print. Defaults to
#'     \code{max(3, getOption("digits") - 3)}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The input object \code{x}, invisibly.
#' @exportS3Method
#' @noRd
print.summary.ire <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  # Print the original function call that created the IRE object
    cat("\nCall:\n") #S: ' ' instead of '\n'
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    # Print method name and number of slices & observations
    cat("Method:\n") #S: ' ' instead of '\n'
    cat(x$method, "with", x$nslices, "slices, n =", x$n)

    # Indicate if weights were used
    if (diff(range(x$weights)) > 0) cat(", using weights.\n") else cat(".\n")

    # Show the sizes of the slices
    cat("\nSlice Sizes:\n") #S: ' ' instead of '\n'
    cat(x$sizes,"\n")

    # Display large-sample marginal dimension tests as a matrix
    cat("\nLarge-sample Marginal Dimension Tests:\n")
    print(as.matrix(x$test), digits = digits)

    # Print estimated directions (basic vectors) for each dimension
    cat("\n")
    cat("\nSolutions for various dimensions:\n")
    for (d in 1:length(x$result)) {
      cat("\n", paste("Dimension =", d), "\n")
      print(x$result[[d]]$B, digits = digits)
    }

    cat("\n")
    invisible(x)
}

#' Extract Orthonormal Basis for Estimated Directions (IRE)
#'
#' Returns an orthonormal basis for the estimated sufficient dimension
#' reduction directions in the original coordinate system, for inverse
#' regression estimator (IRE) objects.
#'
#' @param object An object of class \code{ire}, typically the result of a
#'    fit using IRE.
#' @param numdir Integer. Number of directions to extract. Defaults to the
#'   full number of estimated directions in the object.
#'
#' @return A matrix whose columns form an orthonormal basis for the estimated
#'   central subspace, with appropriate column and row names.
#' @noRd
dr.basis.ire <- function(object, numdir = length(object$result)) {
  # Normalize and flip sign of vectors so that the first elementn is always positive
  fl <- function(z) apply(z, 2, function(x) {
    b <- x / sqrt(sum(x^2)) # Normalize to unit length
    if (b[1] < 0) -1*b else b
    })

  # Transform from Q-space to original coordinate system via backsolve
  ans <- fl(backsolve(dr.R(object), object$result[[numdir]]$B))
  dimnames(ans) <- list(colnames(dr.x(object)),
                        paste("Dir", 1:dim(ans)[2], sep = ""))
  ans
}

#' Project Data onto Estimated IRE Directions
#'
#' Computes the projections of the original predictors onto the estimated
#' inverse regression (IRE) directions.
#'
#' @param object An object of class \code{ire}, typically the result of a fit
#'    using IRE.
#' @param which Integer vector. Specifies which directions to return.
#'   Defaults to all estimated directions.
#' @param x Optional matrix of predictor variables. Defaults to \code{dr.x(object)}.
#'
#' @return A numeric matrix containing the projected data onto the selected
#'    directions.

#' @noRd
dr.direction.ire <- function(object, which = 1:length(object$result),
                             x = dr.x(object)){
  # Determine how many directions are requested
  d <- max(which)

  # Project the centered preditors onto the estimated directions
  scale(x, center = TRUE, scale = FALSE) %*% dr.basis(object, d)
  }

#############################################################################
# partial ire --- see Wen and Cook (in press), Optimal sufficient dimension
# reduction in regressions with categorical predictors. Journal of Statistical
# planning and inference.
#############################################################################

dr.fit.pire <-function(object,numdir=4,nslices=NULL,slice.function=dr.slices,...){
  y <- dr.y(object)
  z <- dr.z(object)
  p <- dim(z)[2]
  if(is.null(object$group)) object$group <- rep(1,dim(z)[1])
  group.names <- unique(as.factor(object$group))
  nw <- table(object$group)
  if (any(nw < p) ) stop("At least one group has too few cases")
  h <- if (!is.null(nslices)) nslices else NCOL(z)
  group.stats <- NULL
  for (j in 1:length(group.names)){
    name <- group.names[j]
    group.sel <- object$group==name
    ans <- dr.compute(z[group.sel,],y[group.sel],method="ire",
                      nslices=h,slice.function=slice.function,
                      tests=FALSE,weights=object$weights[group.sel])
    object$zeta[[j]] <- ans$zeta
    group.stats[[j]] <- ans
  }
  object$group.stats <- group.stats
  numdir <- min(numdir,p-1)
  object$sir.raw.evectors <- dr.compute(z,y,nslices=h,
                                        slice.function=slice.function,
                                        weights=object$weights)$raw.evectors
  class(object) <- c("pire","ire","dr")
  object$indep.test <- dr.test(object,numdir=0,...)
  Gz <- Gzcomp(object,numdir)  # This is the same for all numdir > 0
  ans <- NULL
  for (d in 1:numdir){
    ans[[d]] <- dr.test(object,numdir=d,Gz,...)
    colnames(ans[[d]]$B) <- paste("Dir",1:d,sep="")
  }
  object$Gz <- Gzcomp(object,d,span=ans[[numdir]]$B)
  aa<-c(object, list(result=ans,numdir=numdir))
  class(aa) <- class(object)
  return(aa)
}

dr.iteration.pire <- function(object,Gz,d=2,B=NULL,T=NULL,eps=1.e-6,itmax=200,
                              verbose=FALSE){
  gsolve <- function(a1,a2){  #modelled after ginv in MASS
    Asvd <- svd(a1)
    Positive <- Asvd$d > max(sqrt(.Machine$double.eps) * Asvd$d[1], 0)
    if(all(Positive))
      Asvd$v %*% (1/Asvd$d * t(Asvd$u)) %*% a2
    else Asvd$v[, Positive, drop = FALSE] %*% ((1/Asvd$d[Positive]) *
                                                 t(Asvd$u[, Positive, drop = FALSE])) %*% a2}
  n <- object$cases
  zeta <- object$zeta
  n.groups <- length(zeta)
  p <- dim(zeta[[1]])[1]
  h1 <- 0
  h2 <- NULL
  for (j in 1:n.groups){
    h2[j] <- dim(zeta[[j]])[2]
    h1 <- h1 + h2[j]}
  if (d == 0){
    err <- 0
    for (j in 1:n.groups){
      err <- err + n*sum(forwardsolve(t(Gz[[j]]),as.vector(zeta[[j]]))^2)}
    data.frame(Test=err,df=(p-d)*(h1-d),
               p.value=pchisq(err,(p-d)*(h1-d),lower.tail=FALSE),iter=0)} else {
                 T <- if(is.null(T)) diag(rep(1,p)) else T
                 B <- if(is.null(B)) diag(rep(1,ncol(T)))[,1:d,drop=FALSE] else B
                 fn <- function(B,C){
                   ans <- 0
                   for (j in 1:n.groups){ans <- ans +
                     n * sum( forwardsolve(t(Gz[[j]]),as.vector(zeta[[j]])-
                                             as.vector(T%*%B%*%C[[j]]))^2) }
                   ans}
                 updateC <- function() {
                   C <- NULL
                   for (j in 1:n.groups){ C[[j]]<-
                     matrix( qr.coef(qr(forwardsolve(t(Gz[[j]]),kronecker(diag(rep(1,h2[j])),T%*%B))),
                                     forwardsolve(t(Gz[[j]]),as.vector(zeta[[j]]))), nrow=d)}
                   C}
                 updateB <- function() {
                   for (k in 1:d) {
                     PBk <- qr(B[,-k])
                     a1 <- a2 <- 0
                     for (j in 1:n.groups){
                       alphak <- as.vector(zeta[[j]]-T%*%B[,-k,drop=FALSE]%*%C[[j]][-k,])
                       m1 <-  forwardsolve(t(Gz[[j]]),
                                           t(qr.resid(PBk,t(kronecker(C[[j]][k,],T)))))
                       m2 <- forwardsolve(t(Gz[[j]]),alphak)
                       a1 <- a1 + t(m1) %*% m1
                       a2 <- a2 + t(m1) %*% m2}
                     bk <- qr.resid(PBk, gsolve(a1,a2))
                     bk[is.na(bk)] <- 0  # can use any OLS estimate; eg, set NA to 0
                     bk <- qr.resid(PBk,bk)
                     B[,k] <- bk/sqrt(sum(bk^2))}
                   B}
                 C <- updateC()  # starting values
                 err <- fn(B,C)
                 iter <- 0
                 repeat{
                   iter <- iter+1
                   B <- updateB()
                   C <- updateC()
                   errold <- err
                   err <- fn(B,C)
                   if(verbose==TRUE) print(paste("Iter =",iter,"Fn =",err),quote=FALSE)
                   if ( abs(err-errold)/errold < eps || iter > itmax ) break
                 }
                 B <- T %*% B
                 rownames(B) <- rownames(zeta)
                 list(B=B,summary=data.frame(Test=err,df=(p-d)*(h1-d),
                                             p.value=pchisq(err,(p-d)*(h1-d),lower.tail=FALSE),iter=iter))
               }}

dr.coordinate.test.pire<-function(object,hypothesis,d=NULL,...){
  gamma <- if (inherits(hypothesis, "formula"))
    coord.hyp.basis(object, hypothesis)
  else as.matrix(hypothesis)
  gamma <- dr.R(object)%*%gamma  # Rotate to Q-coordinates:
  p <- object$qr$rank
  r <- p-dim(gamma)[2]
  maxdir <- length(object$result)
  n.groups <- length(object$group.stats)
  h1 <- 0
  h2 <- NULL
  zeta <- object$zeta
  for (j in 1:n.groups){
    h2[j] <- dim(zeta[[j]])[2]
    h1 <- h1 + h2[j]}
  if(is.null(d)){
    H <- qr.Q(qr(gamma),complete=TRUE)[,-(1:(p-r)),drop=FALSE]
    n<-object$cases
    Gz <- object$Gz
    T_H <- 0
    for (j in 1:n.groups){
      m1 <- Gz[[j]] %*% kronecker(diag(rep(1,h2[j])),H)
      m1 <- chol(t(m1) %*% m1)
      T_H <- T_H + n * sum (forwardsolve(t(m1),as.vector(t(H)%*%zeta[[j]]))^2)}
    df <- r*(h1)
    z <- data.frame(Test=T_H,df=df,p.value=pchisq(T_H,df,lower.tail=FALSE))
    z}
  else {
    F0 <-if(maxdir >= d) object$result[[d]]$summary$Test else
      dr.iteration(object,object$Gz,d=d,...)$summary$Test
    F1 <- dr.joint.test(object,hypothesis,d=d,...)$summary$Test
    data.frame(Test=F1-F0,df=r*d,p.value=pchisq(F1-F0,r*d,lower.tail=FALSE))
  }}

Gzcomp.pire <- function(object,numdir,span=NULL){
  Gz <- NULL
  n.groups <- length(object$group.stats)
  pw <- sum(object$group.stats[[1]]$weights)/sum(object$weights)
  Gz[[1]] <- Gzcomp(object$group.stats[[1]],numdir=numdir,span=span)/sqrt(pw)
  if (n.groups > 1){
    for (j in 2:n.groups){
      pw <- sum(object$group.stats[[j]]$weights)/sum(object$weights)
      Gz[[j]] <-
        Gzcomp(object$group.stats[[j]],numdir=numdir,span=span)/sqrt(pw)}}
  Gz
}

"summary.pire" <- function (object, ...)
{   ans <- object[c("call")]
result <- object$result
numdir <- length(result)
tests <- object$indep.test
for (d in 1:numdir) {
  tests <- rbind(tests,result[[d]]$summary)}
rownames(tests) <- paste(0:numdir,"D vs"," > ",0:numdir,"D",sep="")
ans$method <- object$method
ans$nslices <- ans.sizes <- NULL
ans$n <-NULL
for (stats in object$group.stats){
  ans$n <- c(ans$n,stats$n)
  ans$nslices <- c(ans$nslices,stats$slice.info$nslices)
  ans$sizes <- c(ans$sizes,stats$slice.info$slice.sizes)
}
ans$result <- object$result
for (j in 1:length(ans$result)) {ans$result[[j]]$B <- dr.basis(object,j)}
ans$weights <- dr.wts(object)
ans$test <- tests
class(ans) <- "summary.ire"
ans
}

