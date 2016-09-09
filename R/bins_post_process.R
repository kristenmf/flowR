#' splits a datamatrix D into bins with respect to the output of make_bins.
#' @param data_binned  output of make_bins
#' @param w  proportion of zero counts (over samples) allowed in a bin
#' @param  D  the matrix to split. Probably compensated, so you can see what the markers are doing
#' @param  sample_fac  factor indicating which sample type each donor originated from, length(sample_fac) == nrow(D)
#' @return bins data matrices for each bin, combined over all donors
#' @return binnames The bin names
#' @export
get_bins_lineplot <- function(data_binned, w, D, sample_fac) {
#   data_binned <- model_rebin[[5]][[1]]
#   w <- 0.8
#   sample_fac <- D$sample_fac
#   D <- D$D_comp[, markers_index]
  tab <- data_binned$xtab
  tab <- do.call(rbind, tab)
  binnames <- colnames(tab)
  empty <- apply(tab, 2, function(x) length(which(x == 0)))
  which_empty <- which(empty < (nrow(tab) * w))

  grid <- data_binned$grid
  x <- split(data.frame(D), factor(sample_fac))
  data_bins <- mapply(function(x, y) split(x, y), x, grid, SIMPLIFY = FALSE)
  data_bins <- mapply(function(x) x[which_empty], data_bins, SIMPLIFY = FALSE)

  combine_bins <- lapply(1:length(which_empty), function(j) do.call(rbind, lapply(1:length(data_bins), function(i) data_bins[[i]][[j]])))
  return(list(bins = combine_bins, binnames =  binnames[which_empty]))

}


#' Threshold a correlation matrix such that all values less than q are zero.
#' @param mtx A correlation matrix
#' @param q The threshold. Values < q are set to zero.
#' @return The thresholded correlation matrix
mtx_thres <- function(mtx, q) {
  ind <- which(mtx < q, arr.ind = TRUE)
  mtx1 <- array(1, dim = dim(mtx))
  mtx1[ind] <- 0
  return(mtx1)
}


#' Calculates eigenvalue spectrum at different values of q
#' @param C A correlation matrix, not yet thresholded
#' @param seq A vector of thresholds q
#' @return e The eigenvalues of matrix Q
#' @return CQ The ratio of consecutive eigenvalues
CQ <- function(C, seq) {
  L <- lapply(seq, function(.q) {
    mtx_q <- mtx_thres(C, .q)
    diag(mtx_q) <- 0
    d <- colSums(mtx_q)
    w <- which(d > 0)
    D <- diag(1/d[w])
    P <- D %*% mtx_q[w, w]
    Q <- diag(nrow(P)) - P
    e <- eigen(Q)$values
    e <- e[length(e):1]
    CQ <- e[2:(length(e) - 1)]/e[3:length(e)]
    return(list(e = e, CQ = CQ, q = .q))
  }
  )
  return(L)
}

# P_mtx <- function(mtx, q) {
#   mtx_q <- mtx_thres(mtx, q)
#   diag(mtx) <- 0
#   d <- colSums(mtx_q)
#   w <- which(d > 0)
#   D <- diag(1/d[w])
#   P <- D %*% mtx_q[w, w]
#   return(list(d = d, P = P))
# }





























