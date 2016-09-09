#' pre-processing a data matrix for spectral map analysis
#' @param V A matrix with rows as observations, columns as features
#' @return V The log-transformed matrix is row and column centered, where the row and column means have been weighted by linear column and row sums respectively.
#' @export
row_col_center <- function(V) { # V is linear
  Vn <- V/sum(V)
  r <- rowSums(Vn)
  c <- colSums(Vn)
  V <- log(V, 10)
  mr <- rowSums(sweep(V, 2, c, `*`))
  mc <- colSums(sweep(V, 1, r, `*`))
  m <- sum(sweep(sweep(V, 1, r, `*`), 2, c, `*`))
  V1 <- sweep(V, 1, mr, `-`)
  V1 <- sweep(V1, 2, mc, `-`)
  V1 <- V1 + m
  V1 <- sweep(V1, 1, sqrt(r), `*`)
  V1 <- sweep(V1, 2, sqrt(c), `*`)
  return(V1)
}
