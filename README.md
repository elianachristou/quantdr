
<!-- README.md is generated from README.Rmd. Please edit that file -->

# quantdr

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/elianachristou/quantdr.svg?branch=master)](https://travis-ci.com/elianachristou/quantdr)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/elianachristou/quantdr?branch=master&svg=true)](https://ci.appveyor.com/project/elianachristou/quantdr)

[![](https://cranlogs.r-pkg.org/badges/quantdr)](https://cran.rstudio.com/web/packages/quantdr/index.html)
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
n <- 100
p <- 10
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
x <- matrix(rnorm(n * p), n, p)
error <- rnorm(n)
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
#>               [,1]         [,2]         [,3]         [,4]         [,5]
#>  [1,]  0.954434245  0.954157118  0.954284686  0.954102998  0.955280813
#>  [2,]  0.284442247  0.285286950  0.284568946  0.285361444  0.281802100
#>  [3,]  0.041573731  0.041511173  0.041973256  0.042316333  0.040029455
#>  [4,]  0.037601576  0.037516750  0.038225173  0.037488274  0.037592677
#>  [5,] -0.002847888 -0.002493953 -0.003067594 -0.002846182 -0.002591189
#>  [6,] -0.050483133 -0.050151910 -0.049459112 -0.050186911 -0.049889434
#>  [7,]  0.028661501  0.028834868  0.028800844  0.028717673  0.030183956
#>  [8,]  0.017102923  0.017507555  0.017927522  0.017078933  0.017124945
#>  [9,]  0.012950444  0.014204362  0.016656443  0.013684454  0.014670655
#> [10,]  0.034165243  0.034694450  0.035417490  0.035068327  0.033096014

# compare each estimated direction with the true one using the angle between the two subspaces
library(pracma)
for (i in 1:length(taus)) {
print(subspace(out[, i], beta_true) / (pi / 2)) # the angle is measured in radians, so divide by pi/2
  }
#> [1] 0.06105685
#> [1] 0.06102884
#> [1] 0.06173004
#> [1] 0.06121337
#> [1] 0.06126765

# Estimate and plot the conditional quantile function using the new sufficient predictors
library(ggplot2)
newx <- x %*% out
qhat <- as.null()
for (i in 1:length(taus)) {
  qhat <- c(qhat, llqr(newx[, i], y, tau = taus[i])$ll_est)
}  

data1 <- data.frame(rep(dir1, n), rep(y, n), c(newx), rep(taus, each = n), qhat)
names(data1) <- c("dir1", "y", "newx", "quantiles", 'qhat')
ggplot(data1, aes(x = dir1, y = y)) + geom_point(size = 1) + 
  geom_point(aes(x = dir1, qhat), colour = 'red', size = 1) +
  facet_wrap(~quantiles, ncol = 3) + xlab('sufficient direction')
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" style="display: block; margin: auto;" />

Another example using the Boston housing data from the `MASS` library in
`R`.

``` r
library(MASS)
attach(Boston)
# read the data 
y <- medv
x <- cbind(rm, log(tax), ptratio, log(lstat))
n <- length(y)
p <- dim(x)[2]

# plot the estimated coefficient of each predictor variable for multiple quantiles
tau <- seq(0.1, 0.9, by = 0.005)
beta_hat <- matrix(0, p, length(tau))

for (k in 1:length(tau)) {
  out <- cqs(x, y, tau = tau[k])
  beta_hat[, k] <- out$qvectors[, 1:out$dtau] # the suggested dimension of the central quantile subspace is 1
}

data2 <- data.frame(c(t(beta_hat)), rep(tau, p), rep(c('RM', 'log(TAX)', 'PTRATIO', 'log(LSTAT)'), each = length(tau)))
names(data2) <- c('beta_hat', 'tau', 'coefficient')
ggplot(data2, aes(x = tau, y = beta_hat)) + geom_line() + 
  facet_wrap(~coefficient, ncol = 2, scales = "free_y") + 
  ylab('Coefficient') + xlab('Quantile')
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" style="display: block; margin: auto;" />
