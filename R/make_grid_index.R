#' This takes a combined data matrix (generally dimension has already been reduced) which is gridded by spliting each component into l bins of equal width (defined by minimum and maximum of each component. Thus it is generally a good idea to trim outliers). Helper function to be used downstream
#' @param d combined data matrix d, generally dimension reduced
#' @param l The number of bins along each parameter
#' @return a factors with which to split d
make_grid_index <- function(d, l) { # here d is all the datasets together, 'average face', make one single reference grid, then return each individual grid. 
  #d <- D_transform[, 1:2]
  #l <- 10
  d <- data.frame(d)
  grid_index <- apply(d, 2, function(x) x - min(x))
  grid_index <- apply(grid_index, 2, function(x) floor(x/max(x)*l))
  grid_index <- apply(grid_index, 2, function(x) {x[x == max(x)] <- (max(x) - 1); return(x)})
  grid_index <- if (ncol(grid_index) == 1) {factor(grid_index)} else {factor(apply(grid_index, 1, function(x) paste(x, sep = '', collapse = '-')))}
  #grid_index_split <- split(grid_index, sample_fac)
  return(grid_index)
}

