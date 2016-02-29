#' After binning the data, make a table with the nearly empty bins combined into a 'leftover bin'. often dramatically reduces size of the table without much information loss. 
#' @param data_binned output of make_bins
#' @param w the 'vacancy rate', a number between 0 and 1. The final table will have no column (bin) with a greater proportion than w of donors with zero events in that bin. 
#' @param C The row names - donor sample type
#' @return The reduced data matrix. Each row (donor) sums to one. The first column is the number of 'leftover' events, i.e. the events that don't fit into the bins that are common to most samples. 
get_xtab <- function(data_binned, w, C) {
  #data_binned <- bins_k2[[1]]
  #w <- 0.2
  #C <- sample_names1
  tab <- data_binned$xtab
  tab <- do.call(rbind, tab)
  empty <- apply(tab, 2, function(x) length(which(x == 0)))
  which_empty <- which(empty < (length(C) * w))
  
  donor_coords <- tab
  donor_coords <- donor_coords[, which_empty]
  if (ncol(donor_coords) == (ncol(tab) - 1)) {
    donor_coords <- cbind(tab[, -which_empty], donor_coords)
  }
  if (ncol(donor_coords) < (ncol(tab) - 1)) {
    donor_coords <- cbind(rowSums(tab[, -which_empty]), donor_coords)
  }
  donor_coords <- t(apply(donor_coords, 1, function(x) x/sum(x)))
  rownames(donor_coords) <- C
  return(donor_coords)
}
