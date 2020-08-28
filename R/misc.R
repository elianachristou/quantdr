bic_d <- function(lambdas, n, p) {
  gn <- as.null(p)
  for (i in seq_len(lambdas)) {
    gn[i] <- n * sum((lambdas[1:i])^2) / sum((lambdas)^2)
                     - 2 * (n^ (3 / 4) / p) * i * (i + 1) / 2
  }
  return(which(gn == max(gn)))
}

center <- function(x) {
  return(t(t(x) - apply(x, 2, mean)))
  }

matpower <- function(a, alpha) {
  a <- (a + t(a)) / 2
  tmp <- eigen(a)
  return(tmp$vectors %*% diag((tmp$values)^alpha) %*%
           t(tmp$vectors))
  }
