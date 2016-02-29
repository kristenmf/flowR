#' puts a data matrix (compensated, uncompensated) into a format for lineplots
#' ncol(X) == length(markers)
#' @param X data matrix with observations/events as rows, features/channels as columns
#' @param markers The column/channel names, e.g. the biomarkers used
#' @return a data frame for making lineplots/parallel coordinate plots
make_lineplot_data <- function(X, markers) {
  spec = t(X)[1:(nrow(X) * ncol(X))]
  ch = factor(rep(markers, nrow(X)), levels = markers, ordered = TRUE)
  bin <- factor(rep(1:nrow(X), rep(ncol(X), nrow(X))))
  df <- data.frame(spec = spec, ch = ch, bin = bin)
  return(df)
}
