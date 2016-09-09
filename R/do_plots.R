#' makes lineplots
#' @param bins  a list of data matrices, corresponding to bins
#' @param dir  the directory where the plots should be made
#' @param markers  vector of strings indicating the markers
#' @param threshold  the minimum sized bin that should be plotted
#' @param y1, y2  vertical plot limits
#' @return a plot in the desired directory (png format)
#' @export
do_plots <- function(bins, dir, markers, threshold, y1, y2) {
  for (i in 1:length(bins[[1]])) {
    tmp <- bins[[1]][[i]]
    colnames(tmp) <- markers
    if (nrow(tmp) > threshold) {
      if (nrow(tmp) > 10000) {line_data <- make_lineplot_data(tmp[sample(nrow(tmp), 10000), ], markers)} else {line_data <- make_lineplot_data(tmp, markers)}
      size <- 10/log(nrow(tmp)/length(markers), 1.05)
      fn <- paste(dir, 'bin_', bins[[2]][[i]], '.png', sep = '')
      png(fn)
      print(ggplot(line_data, aes(x = ch, y = spec, group = bin)) + geom_path(colour = rgb(0.2, 0.2, 0.2, 0.4), size = size) + ylim(y1, y2) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)))
      dev.off()
    }
  }
}

