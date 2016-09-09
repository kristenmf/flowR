#' trim the data before starting
#' @param Data A list of data matrices
#' @param lim1, lim2 The upper and lower limits for trimming on the interval (0, 1) (quantiles)
#' @return The list of matrices is combined into one single matrix. Trimming is performed for each column. Rows are filtered out if at least one value falls outside the bounds set by lim1, lim2. The filted data matrix is returned.  
outlier_pretrim <- function(Data, lim1, lim2) {
  Q <- apply(do.call(rbind, Data), 2, function(x) quantile(x, c(lim1, lim2)))
  trim <- mapply(function(x) rowSums(sapply(1:ncol(x), function(i) x[, i] > Q[1, i] & x[, i] < Q[2, i])), Data, SIMPLIFY = FALSE)
  trim <- mapply(function(x) x == ncol(Q), trim, SIMPLIFY = FALSE)
  d_t <- mapply(function(x, y) x[y, ], Data, trim, SIMPLIFY = FALSE)
  return(d_t)
}
