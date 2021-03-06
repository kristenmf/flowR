#' Calculates common breaks of all datasets for each channel
#' @param X list of matrices, rows are events, columns are parameters
#' @return list of intervals for each parameter
#' @export
make_breaks <- function(X) {
  lapply(1:ncol(X[[1]]), function(i) {
    x <- unlist(mapply(function(x) x[, i], X))
    c(min(x), min(x) + ((max(x) - min(x))/2), max(x))
  }
  )
}

#' Converts each matrix into a burt matrix, consisting of marker combination and frequency it occurs in that sample
#' @param X list of data frames, rows are events, columns are parameters. MUST be dataframes.
#' @param breaks, output from make_breaks or custom if desired.
#' @param markers vector of character strings indicating marker names, in same order as they appear in the columns of the matrices in X
#' @return list of dataframes, with markers coded as 1 for 'positive' and 0 for 'negative', and additional column 'counts' indicating the proportion of events with that marker combination. Only combinations with non-zero proportions are included. Column 'which' gives a unique index to each sample.
#' @export
make_burt <- function(X, breaks, markers) {
  tmp1 <- lapply(1:length(X), function(i) {
    x <- mapply(function(x, y) cut(x, y, FALSE, TRUE), X[[i]], breaks)
    tmp <- ddply(as.data.frame(x), markers, nrow)
    colnames(tmp)[ncol(X[[1]]) + 1] <- 'counts'
    tmp$counts <- tmp$counts/sum(tmp$counts)
    tmp$which <- i
    tmp
  }
  )
  return(tmp1)
}

#' Finds the ratio of marker frequency in a sample to average marker frequency. This is counting 'total positive' for each marker, e.g. the bulk behaviour. Events may be double-counted, e.g. cells which are IL17-IL22 double positive are counted twice in IL17 total counts and IL22 total counts.
#' @param burt , output of make_burt
#' @return A matrix, where each row is a vector of the ratio of a sample's total marker frequency to the average total marker frequency (averaged over all samples). Columns are markers, additionally the final column is the number of unstained cells. When a marker combination doesn't appear in the sample, a value of 1 is given.
#' @export
# make_sample_coexp <- function(burt) {
#   burt <- mapply(function(x) rbind(x, c(rep(1, ncol(x) - 2), 0, x$which[1])), burt, SIMPLIFY = FALSE)
#   burt_all <- do.call(rbind, burt)
#   tot_pos <- apply(burt_all[, 1:(ncol(burt_all) - 2)], 2, function(y) by(burt_all$counts, y, sum)[2]/length(burt))
#   unst <- lapply(burt, function(.x) .x[rowSums(.x[, 1:(ncol(burt_all) - 2)]) == (ncol(burt_all) - 2), 'counts'])
#   unst <- mapply(sum, unst)
#   unst[unst == 0] <- 1
#
#   expression_donor <-
#     lapply(burt, function(.x) {
#       apply(.x[, 1:(ncol(burt_all) - 2)], 2, function(y) by(.x$counts, y, sum)[2])/tot_pos
#     })
#   expression_donor <- do.call(rbind, expression_donor)
#   expression_donor[is.na(expression_donor)] <- 1
#   expression_donor <- cbind(expression_donor, unst)
#   return(expression_donor)
# }

make_sample_coexp <-
  function(burt) {
    burt <- mapply(function(x) rbind(x, c(rep(1, ncol(x) - 2), 0, x$which[1])), burt, SIMPLIFY = FALSE)
    burt_all <- do.call(rbind, burt)
    tot_pos <- apply(burt_all[, 1:(ncol(burt_all) - 2)], 2, function(y) by(burt_all$counts, y, sum)[2]/length(burt))
    unst <- lapply(burt, function(.x) .x[rowSums(.x[, 1:(ncol(burt_all) - 2)]) == (ncol(burt_all) - 2), 'counts'])
    unst <- mapply(sum, unst)

    expression_donor <-
      lapply(burt, function(.x) {
        apply(.x[, 1:(ncol(burt_all) - 2)], 2, function(y) by(.x$counts, y, sum)[2])
      })
    expression_donor <- do.call(rbind, expression_donor)
    expression_donor[is.na(expression_donor)] <- 0
    expression_donor <- cbind(expression_donor, unst)
    return(expression_donor)
  }


#' For every possible pair of markers, find the ratio of the frequency in each sample to the average frequency over all samples.
#' @param burt output from make_burt
#' @param markers vector of character strings, marker names in same order as they appear in columns of burt
#' @return A matrix, where each row is a vector of the ratio of a sample's  frequency to the average  frequency (averaged over all samples). Columns are marker pairs. When a marker combination doesn't appear in the sample, a value of 1 is given.
#' @export
make_cell_coexp <-
  function(burt, markers) {

    burt_all <- do.call(rbind, burt)
    N <- ncol(burt_all) - 2

    double_pos <- t(sapply(1:length(burt), function(k) {
      unlist(lapply(1:(N - 1), function(i) sapply((i+1):N, function(j) {
        sum(burt[[k]][burt[[k]][, i] == 2 & burt[[k]][, j] == 2, 'counts'])
      }
      )
      )
      )
    }
    )
    )
    cn <- unlist(lapply(1:(N - 1), function(i) sapply((i+1):N, function(j) paste(markers[j], markers[i], sep = '_'))))
    colnames(double_pos) <- cn
    return(double_pos)

  }



#' For every possible pair of markers, find the ratio of the frequency in each sample to expected frequency given no association between the markers.
#' @param burt output from make_burt
#' @param markers vector of character strings, marker names in same order as they appear in columns of burt
#' @return A matrix, where each row is a vector of the ratio of a sample's  frequency to the expected frequency given no association between the markers. Columns are marker pairs.
#' @export
make_twoD_marginal <- function(burt, markers) {

  N <- ncol(burt[[1]]) - 2
  P <- t(sapply(1:length(burt), function(k) {
    unlist(lapply(1:(N - 1), function(i) sapply((i+1):N, function(j) {
      sum(burt[[k]][burt[[k]][, i] == 2 & burt[[k]][, j] == 2, 'counts'])/(sum(burt[[k]][burt[[k]][, i] == 2, 'counts']) * sum(burt[[k]][burt[[k]][, j] == 2, 'counts']))
    }
    )
    )
    )
  }
  )
  )

  cn <- unlist(lapply(1:(N - 1), function(i) sapply((i+1):N, function(j) paste(markers[j], markers[i], sep = '_'))))

  colnames(P) <- cn
  P[is.na(P)] <- 1
  P[P == 0] <- 1
  return(P)

}

#' Looks at the relationship between correlation of total marker frequency and cell-level coexpression
#' @param sample_coexp output of make_sample_coexp
#' @param cell_coexp output of make_cell_coexp
#' @param markers vector of character strings, marker names in same order as they appear in columns of burt
#' @return A dataframe. For each pair of markers: R squared value of multiple regression of cell coexpression frequency on both the total frequencies. Correlation between both total frequencies.
#' @export
make_rsq <- function(sample_coexp, cell_coexp, markers) { # can be donor-centered versions if necessary
  N <- ncol(L) - 1
  rsq <- lapply(1:(N - 1), function(i) lapply((i+1):N, function(j) {
    df <- data.frame(x1 = sample_coexp[, i], x2 = sample_coexp[, j], y = cell_coexp[, which(paste(colnames(sample_coexp)[j], colnames(sample_coexp)[i], sep= '_') == colnames(cell_coexp))])
    if (ncol(df) == 3) {
      lmodel <- lm(y ~ x1 + x2, df)
      summary(lmodel)$r.squared} else {0}
  }
  )
  )

  rsq <- unlist(rsq)

  cn <- unlist(lapply(1:(N - 1), function(i) sapply((i+1):N, function(j) paste(markers[j], markers[i], sep = '_'))))
  C <- cor(do.call(rbind, (by(L, factor(sample_fac$stim), colMeans))))
  cor_markers <- unlist(sapply(1:(N-1), function(i) C[(i+1):N, i]))
  names(cor_markers) <- cn

  return(list(rsq = rsq, cor = cor_markers))
}

#' Finds the cross-correlation between total frequencies and cell coexpression.
#' @param sample_coexp output from make_sample_coexp
#' @param cell_coexp output from make_cell_coexp
#' @return Matrix of cross-correlation values. Rows correspond to total frequencies, columns correspond to cell coexpression frequencies.
#' @export
make_cross_correlation <- function(sample_coexp, cell_coexp) {
  cor(cbind(sample_coexp, cell_coexp))[1:ncol(sample_coexp), (ncol(sample_coexp)+1):(ncol(sample_coexp) + ncol(cell_coexp))]
}

#' Finds the matrix residual after singular value decomposition components are removed.
#' @param X matrix
#' @param k Vector of components to be removed
#' @return Residual vector
#' @export
svd_residuals <- function(X, k) {
  s <- svd(X)
  u <- s$u
  u[, k] <- 0
  v <- s$v
  v[, k] <- 0
  d <- s$d
  d[k] <- 0
  s1 <- u %*% diag(d) %*% t(v)
  colnames(s1) <- colnames(X)
  rownames(s1) <- rownames(X)
  return(s1)

}


#' summarises the frequency of events that are at least single positive, double positive etc for each marker. Pooled over all donors
#' @param burt  A list of burt matrices, output from make_burt
#' @param N N = 1 for at least single positive, N = 2 for at least double positive, N = 3 for at least triple positive etc. For example, any event positive for IL17 is at least single positive; any event positive for IL17 and at least one other marker is double positive, any event positive for IL17 and at least two other markers is triple positive etc
#' @export
make_summary <- function(burt, N) {
  burt_all <- do.call(rbind, burt)
  N1 <- ncol(burt_all) - 2
  N2 <- ncol(burt_all) - 2 + N - 1
  apply(burt_all[, 1:N1], 2, function(x) sum(burt_all$counts[rowSums(burt_all[, 1:N1]) > N2 & x == 2])/sum(burt_all$counts))
}


#' Makes a square matrix indicating overlapping marker pairs
#' @param single a vector of character strings, containing the set of markers from which the pairs are drawn.
#' @param double a vector of character strings, containing the set of pairs.
#' @param sep the separating character. For example, single = c('A', 'B'), double = c('A_B') and sep = '_'
#' @return an square matrix with dimension equal to length of double. Matrix element equals 1 if the two pairs share a marker, zero otherwise.
#' @export
make_mask <- function(single, double, sep) {
  cn <- single
  indicator <- sapply(double, function(.m) as.numeric(sapply(cn, function(.c) .m %in% strsplit(.c, sep)[[1]])))
  indicator_square <- indicator %*% t(indicator)
  diag(indicator_square) <- 1
  colnames(indicator_square) <- cn
  rownames(indicator_square) <- cn
  return(indicator_square)

}


#' takes a burt table for one donor and a data frame describing the double positive pairs. returns a vector specifying whether a row in the burt table belongs to the subset described by the dataframe.
#' @param donor a burt table for a single donor, output of make_burt
#' @param pairs A dataframe or matrix summarising marker pairs of interest. Dataframe has two columns, one for each element of the pair.
#' @return A vector with length equal to number of rows in donor, encoding whether marker combination is in subset or not. Vector element equals 2 when in the subset, and 1 when not.
#' @export
doublepos_events <- function(donor, pairs) {
  as.integer(Reduce('+', lapply(1:nrow(pairs), function(i) donor[[pairs[i, 1]]] == 2 & donor[[pairs[i, 2]]] == 2)) > 0) + 1
}

#' takes a burt table for one donor and a data frame describing the double positive pairs. returns a vector specifying whether a row in the burt table belongs to the subset described by the dataframe.
#' @param donor a burt table for a single donor, output of make_burt
#' @param pairs A dataframe or matrix summarising marker pairs of interest. Dataframe has two columns, one for each element of the pair.
#' @return A vector with length equal to number of rows in donor, encoding whether marker combination is in subset or not. Vector element equals 2 when in the subset, and 1 when not.
#' @export
singlepos_events <- function(donor, marker, bulkmarkers) { # takes a burt table for one donor and a marker name. returns a vector specifying whether a row in the burt table is single positive for that marker
  as.integer(rowSums(donor[bulkmarkers]) == 8 & donor[[marker]] == 2) + 1
}




