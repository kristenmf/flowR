#' take samples from data. Always use with set.seed for reproducibility
#' @param Data is a list of uncompensated data matrices
#' @param D_comp is a list of compensated data matrices
#' @param S is the number of samples to be taken
#' @return samples A list of sample indices applied to each matrix
#' @return sample_fac a factor to split combined data matrix into individuals
#' @return D combined sampled uncompensated matrix
#' @return D_comp combined sampled compensated matrix
#' @export
sample_data <- function(Data, D_comp, S) {
  samples <- mapply(function(x) if (nrow(x) < S) {1:nrow(x)} else {sample(nrow(x), S)}, Data, SIMPLIFY = FALSE)
  sample_fac <- rep(1:length(Data), mapply(length, samples))
  D <- do.call(rbind, mapply(function(x, y) x[y, ], Data, samples, SIMPLIFY = FALSE))
  Dc <- do.call(rbind, mapply(function(x, y) x[y, ], D_comp, samples, SIMPLIFY = FALSE))
  return(list(samples = samples, sample_fac = sample_fac, D = D, D_comp = Dc))
}
