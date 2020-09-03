
<!-- README.md is generated from README.Rmd. Please edit that file -->

# quantdr

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/elianachristou/quantdr.svg?branch=master)](https://travis-ci.com/elianachristou/quantdr)
<!-- badges: end -->

# quantdr: Dimension Reduction Techniques for Conditional Quantiles

The R package ‘quantdr’ performs dimension reduction techniques for
conditional quantiles by estimating the fewest linear combinations of X
that contain all the information on that function. For details of the
methodology, see [Christou, E. (2020) Central quantile subspace.
*Statistics and Computing*, 30,
677–695](https://link.springer.com/article/10.1007/s11222-019-09915-8).

The main function of the package is `cqs`, which estimates the
directions of the central quantile subspace. Once the directions are
determined, one can form the new sufficient predictors and estimate the
conditional quantile function using `llqr`.

## Installation

You can install the released version of quantdr from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("quantdr") 
```

and the development version from [GitHub](https://github.com/)
with:

<!-- You can install the development version of quantdr from [GitHub](https://github.com/) with:  -->

``` r
# install.packages("devtools")
devtools::install_github("elianachristou/quantdr")
```

## Example

This is a basic example which shows you how to solve the problem.

``` r
library(quantdr)

## basic example code - a homoscedastic single-index model

# Setting
set.seed(1234)
n <- 100; p <- 10; taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
x <- matrix(rnorm(n * p), n, p); error <- rnorm(n)
y <- 3 * x[, 1] + x[, 2] + error

# true direction that spans each central quantile subspace
beta_true <- c(3, 1, rep(0, p - 2))
beta_true / sqrt(sum(beta_true^2))
#>  [1] 0.9486833 0.3162278 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#>  [8] 0.0000000 0.0000000 0.0000000

# sufficient direction
dir1 <- x %*% beta_true

# Estimate the directions of each central quantile subspace
# Since dtau is known to be one, the algorithm will produce only one vector
out <- matrix(0, p, length(taus))
for (i in 1:length(taus)) {
  out[, i] <- cqs(x, y, tau = taus[i], dtau = 1)$qvectors
}
out
#>               [,1]         [,2]          [,3]         [,4]         [,5]
#>  [1,]  0.949525419  0.951345485  0.9481249240  0.949893722  0.951285567
#>  [2,]  0.297675031  0.290211374  0.3036884256  0.295630714  0.280349811
#>  [3,]  0.066064403  0.063422461  0.0539648971  0.042117726  0.047074410
#>  [4,]  0.026048171  0.020410164  0.0255546948  0.049771597  0.040678843
#>  [5,]  0.002763456  0.008308760  0.0043920967  0.008982865  0.003477829
#>  [6,] -0.050276673 -0.061732986 -0.0534880659 -0.046628222 -0.057344442
#>  [7,]  0.032939936  0.035874469  0.0284903647  0.016886469  0.015840191
#>  [8,]  0.016294499  0.016271051  0.0228861044  0.039207685  0.061479219
#>  [9,] -0.001967384 -0.003376335 -0.0002880524  0.012076001  0.026121381
#> [10,]  0.029295759  0.028931331  0.0324252190  0.042780536  0.067642983

# compare each estimated direction with the true one using the angle between the two subspaces
library(pracma)
for (i in 1:length(taus)) {
print(subspace(out[, i], beta_true) / (pi / 2)) # the angle is measured in radians, so divide by pi/2
  }
#> [1] 0.06412042
#> [1] 0.06801387
#> [1] 0.06038418
#> [1] 0.06597454
#> [1] 0.08488814

# Estimate and plot the conditional quantile function using the new sufficient predictors
newx <- x %*% out
par(mfrow=c(2,3))
for (i in 1:length(taus)) {
  plot(dir1, y, xlab = "sufficient direction", ylab = "y", main = taus[i], pch = 16)
  qhat <- llqr(newx[, i], y, tau = taus[i])$ll_est
  points(dir1, qhat, pch = 16, col = "red")
}
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Another example using the Boston housing data from the `MASS` library in
`R`.

``` r
library(MASS)
attach(Boston)
# read the data 
y <- medv
x <- cbind(rm, log(tax), ptratio, log(lstat))
n <- length(y); p <- dim(x)[2]

# plot the estimated coefficient of each predictor variable for multiple quantiles
tau <- seq(0.1, 0.9, by = 0.005)
beta_hat <- matrix(0, p, length(tau))

for (k in 1:length(tau)) {
  out <- cqs(x, y, tau = tau[k])
  beta_hat[, k] <- out$qvectors[, 1:out$dtau] # the suggested dimension of the central quantile subspace is 1
}
   
par(mfrow=c(2,2))
plot(tau, beta_hat[1, ], type = 'l', xlab = 'Quantile', main = 'RM', ylab = 'Coefficient')
plot(tau, beta_hat[2, ], type = 'l', xlab = 'Quantile', main = 'log(TAX)', ylab = 'Coefficient')
plot(tau, beta_hat[3, ], type = 'l', xlab = 'Quantile', main = 'PTRATIO', ylab = 'Coefficient')
plot(tau, beta_hat[4, ], type = 'l', xlab = 'Quantile', main = 'log(LSTAT)', ylab = 'Coefficient')
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />
