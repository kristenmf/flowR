#' This takes a combined data matrix (generally dimension has already been reduced) which is gridded by spliting each component into l bins of equal width (defined by minimum and maximum of each component. Thus it is generally a good idea to trim outliers). Helper function for make_grid, not generally used on its own.
#' @param d combined data matrix d, generally dimension reduced
#' @param l The number of bins along each parameter
#' @return a factors with which to split d
#' @export
make_grid_index <- function(d, l) {
  d <- data.frame(d)
  grid_index <- apply(d, 2, function(x) x - min(x))
  grid_index <- apply(grid_index, 2, function(x) floor(x/max(x)*l))
  grid_index <- apply(grid_index, 2, function(x) {x[x == max(x)] <- (max(x) - 1); return(x)})
  grid_index <- if (ncol(grid_index) == 1) {factor(grid_index)} else {factor(apply(grid_index, 1, function(x) paste(x, sep = '', collapse = '-')))}
  #grid_index_split <- split(grid_index, sample_fac)
  return(grid_index)
}

