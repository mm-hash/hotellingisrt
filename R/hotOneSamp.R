hotOneSamp <- function(x, mu = rep(0, ncol(x))) {
  n <- nrow(x)
  p <- ncol(x)
  bar.x <- colMeans(x)
  Sigma <- var(x)
  Sigma.inv <- solve(Sigma)
  T2 <- n * t(bar.x - mu) %*% Sigma.inv %*% (bar.x - mu)
  F_stat = (n - p) / (p * (n - 1)) * T2
  pval = pf(F_stat, p, n - p, lower.tail = F)
  return(list(T2=T2, F_statistic = F_stat, pvalue = pval))
}

#' hotelling one sample T^2 test
#' 
#' takes multivariate data, with one independent vector. 
#' @param x A matrix
#' @returns T^2 statistic, F statistic and the p-value.

