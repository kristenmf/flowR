#' bins a dimension-reduced pooled dataset. Can be with PCA or SMA
#' @param DONORS_reduced pooled dataset decomposed by PCA or SMA. Generally keep all components, because this is chosen within this function
#' @param sample_fac index of which sample each event originated from, length(sample_fac) == nrow(D), as given by sample_data
#' @param k  the number of components to include in analysis
#' @param l  number of bins along each component
#' @return grid a list of factors, one for each donor, corresponding to the binned events
#' @return xtab a table with number of events in each bin, one for each donor.
#' @export
make_bins <- function(DONORS_reduced, sample_fac, k, l) {
  DONORS_grid <- make_grid_index(DONORS_reduced[, 1:k],  l)
  DONORS_grid <- factor(DONORS_grid, levels = unique(DONORS_grid))
  DONORS_grid <- split(DONORS_grid, factor(sample_fac))
  DONORS_crosstab <- mapply(table, DONORS_grid, SIMPLIFY = FALSE)
  return(list(grid = DONORS_grid, xtab = DONORS_crosstab))

}
