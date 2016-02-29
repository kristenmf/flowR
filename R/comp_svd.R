#' Spectral Map Analysis
#' Input V must be linear
#' @param V a matrix on the linear scale
#' @return The singular value decomposition of the double centered matrix
comp_svd <- function(V) {
  V1 <- row_col_center(V)
  s <- svd(V1)
  return(s)
}
