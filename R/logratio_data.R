#' transforms compositional data for a biplot. First, each row is centered on its geometrical mean, followed by log transformation. This is equivalent to the logratio of every single pair of variables. Then the data is scaled.
#' @param data Matrix of compositions. Each row sums to 1
#' @return data_lr Matrix of logratios, column centered. Suitable for a biplot.
logratio_data <- function(data) {
data_lr <- t(apply(data, 1, function(x) x/sum(x)))
gm <- apply(data_lr, 1, function(x) exp(mean(log(x))))
data_lr <- sweep(data_lr, 1, gm, '/')
data_lr <- scale(log(data_lr))
return(data_lr)
}

#' replaces zero with a small value
#' @param X matrix of compositions
#' @return matrix of compositions, with zero replaced by a small value and rescaled
replace_zero <- function(X) {
  X[X == 0] <- min(X[X > 0])/100
  X <- t(apply(X, 1, function(x) x/sum(x)))
  return(X)
}


