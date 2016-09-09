#' Calculates the Aitchison (compositional) distance between samples.
#' @param table of counts in each bin for each donor (row). Rows must sum to one. 
#' @return a distance matrix between donors. 
aDist_mtx <- function(table) { # table > 0
  gm <- apply(table, 1, function(x) exp(mean(log(x))))
  table_transformed <- log(sweep(table, 1, gm, `/`))
  ad <- dist(table_transformed)
  return(ad)
}
