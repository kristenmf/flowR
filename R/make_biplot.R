#' make a biplot from a singular value decomposition
#' @param s A singular value decomposition, returned from svd
#' @param alpha either 0 or 1. Commonly alpha = 1, which favours display of the individuals (form biplot). alpha = 0 (covariance biplot) favours display of the variables
#' @param dim The number of dimensions to retain
#' @return F Cartesian coordinates of individuals. Columns correspond to components
#' @return G Cartesian coordinates of variables. Columns correspond to components. Generally represented as rays emanating from origin.
#' @export
make_biplot <- function(s, alpha, dim) {
  F <- s$u[, 1:dim] * matrix(rep(s$d[1:dim], nrow(s$u)), byrow = TRUE, ncol = dim)^alpha
  G <- s$v[, 1:dim] * matrix(rep(s$d[1:dim], nrow(s$v)), byrow = TRUE, ncol = dim)^(1-alpha)
  return(list(F = F, G =G))
}


# biplot <- function(data, alpha) {
#   S <- svd(data)
#   #alpha <- 1
#   f <- sweep(S$u, 2, S$d^alpha, `*`)
#   g <- sweep(S$v, 2, S$d^(1 - alpha), `*`)
#   return(list(f = f, g = g))
# }
