# library(flowCore)
# library(ggplot2)
# library(pheatmap)
# library(robCompositions)
# library(e1071)
# library(vegan)
# library(energy)
# 
# 
# ###################
# # data preprocess
# 
# #' trim the data before starting
# #' @param Data A list of data matrices
# #' @param lim1, lim2 The upper and lower limits for trimming, expressed as percentiles. Computed column-wise
# #' @return The list of matrices is combined into one single matrix. Then rows are filtered out if at least one value falls outside the bounds set by lim1, lim2. The filted data matrix is returned.  
# outlier_pretrim <- function(Data, lim1, lim2) {
#   Q <- apply(do.call(rbind, Data), 2, function(x) quantile(x, c(lim1, lim2)))
#   trim <- mapply(function(x) rowSums(sapply(1:ncol(x), function(i) x[, i] > Q[1, i] & x[, i] < Q[2, i])), Data, SIMPLIFY = FALSE)
#   trim <- mapply(function(x) x == ncol(Q), trim, SIMPLIFY = FALSE)
#   d_t <- mapply(function(x, y) x[y, ], Data, trim, SIMPLIFY = FALSE)
#   return(d_t)
# }
# 
# #' take samples from data. Always use with set.seed for reproducibility
# #' @param Data is a list of uncompensated data matrices
# #' @param D_comp is a list of compensated data matrices
# #' @param S is the number of samples to be taken
# #' @return samples A list of sample indices applied to each matrix
# #' @return sample_fac a factor to split combined data matrix into individuals
# #' @return D combined sampled uncompensated matrix
# #' @return D_comp combined sampled compensated matrix
# sample_data <- function(Data, D_comp, S) {
#   samples <- mapply(function(x) if (nrow(x) < S) {1:nrow(x)} else {sample(nrow(x), S)}, Data, SIMPLIFY = FALSE)
#   sample_fac <- rep(1:length(Data), mapply(length, samples))
#   D <- do.call(rbind, mapply(function(x, y) x[y, ], Data, samples, SIMPLIFY = FALSE))
#   Dc <- do.call(rbind, mapply(function(x, y) x[y, ], D_comp, samples, SIMPLIFY = FALSE))
#   return(list(samples = samples, sample_fac = sample_fac, D = D, D_comp = Dc))
# }
# 
# # # trim the top and bottom N outliers of each channel after sampling (for when there is only a low cell count e.g. t-cells)
# # # update sample_fac accordingly
# # outlier_trim <- function(D, D_comp, sample_fac, N) {
# #   outlier <- apply(D, 2, rank)
# #   outlier <- rowSums(apply(outlier, 2, function(x) x > N & x < (length(x) - N))) == ncol(D)
# #   D <- D[outlier, ]
# #   D_comp <- D_comp[outlier, ]
# #   sample_fac <- sample_fac[outlier]
# #   return(list(outlier = outlier, D = D, sample_fac = sample_fac, D_comp = D_comp))
# # }
# 
# ###############
# # some functions
# 
# #' pre-processing a data matrix for spectral map analysis
# #' @param V A matrix with rows as observations, columns as features
# #' @return V The log-transformed matrix is row and column centered, where the row and column means have been weighted by linear column and row sums respectively. 
# row_col_center <- function(V) { # V is linear
#   Vn <- V/sum(V)
#   r <- rowSums(Vn)
#   c <- colSums(Vn)
#   V <- log(V, 10)
#   mr <- rowSums(sweep(V, 2, c, `*`))
#   mc <- colSums(sweep(V, 1, r, `*`))
#   m <- sum(sweep(sweep(V, 1, r, `*`), 2, c, `*`))
#   V1 <- sweep(V, 1, mr, `-`)
#   V1 <- sweep(V1, 2, mc, `-`)
#   V1 <- V1 + m
#   V1 <- sweep(V1, 1, sqrt(r), `*`)
#   V1 <- sweep(V1, 2, sqrt(c), `*`)
#   return(V1)
# }
# 
# #' Spectral Map Analysis
# #' Input V must be linear
# #' @param V a matrix on the linear scale
# #' @return The singular value decomposition of the double centered matrix
# comp_svd <- function(V) {
#   V1 <- row_col_center(V)
#   s <- svd(V1)
#   return(s)
# }
# 
# #' make a biplot from a singular value decomposition
# #' @param s A singular value decomposition, returned from svd
# #' @param alpha either 0 or 1. Commonly alpha = 1, which favours display of the individuals (form biplot). alpha = 0 (covariance biplot) favours display of the variables
# #' @param dim The number of dimensions to retain
# #' @return F Cartesian coordinates of individuals. Columns correspond to components
# #' @return G Cartesian coordinates of variables. Columns correspond to components. Generally represented as rays emanating from origin. 
# make_biplot <- function(s, alpha, dim) {
#   F <- s$u[, 1:dim] * matrix(rep(s$d[1:dim], nrow(s$u)), byrow = TRUE, ncol = dim)^alpha
#   G <- s$v[, 1:dim] * matrix(rep(s$d[1:dim], nrow(s$v)), byrow = TRUE, ncol = dim)^(1-alpha)
#   return(list(F = F, G =G))
# }
# 
# ###############
# 
# #' This takes a combined data matrix (generally dimension has already been reduced) which is gridded by spliting each component into l bins of equal width (defined by minimum and maximum of each component. Thus it is generally a good idea to trim outliers). Helper function to be used downstream
# #' @param d combined data matrix d, generally dimension reduced
# #' @param l The number of bins along each parameter
# #' @return a factors with which to split d
# make_grid_index <- function(d, l) { # here d is all the datasets together, 'average face', make one single reference grid, then return each individual grid. 
#   #d <- D_transform[, 1:2]
#   #l <- 10
#   d <- data.frame(d)
#   grid_index <- apply(d, 2, function(x) x - min(x))
#   grid_index <- apply(grid_index, 2, function(x) floor(x/max(x)*l))
#   grid_index <- apply(grid_index, 2, function(x) {x[x == max(x)] <- (max(x) - 1); return(x)})
#   grid_index <- if (ncol(grid_index) == 1) {factor(grid_index)} else {factor(apply(grid_index, 1, function(x) paste(x, sep = '', collapse = '-')))}
#   #grid_index_split <- split(grid_index, sample_fac)
#   return(grid_index)
# }
# 
# 
# #' bins a dimension-reduced pooled dataset. Can be with PCA or SMA
# #' @param DONORS_reduced pooled dataset decomposed by PCA or SMA. Generally keep all components, because this is chosen within this function
# #' @param sample_fac index of which sample each event originated from, length(sample_fac) == nrow(D), as given by sample_data
# #' @param k  the number of components to include in analysis
# #' @param l  number of bins along each component
# #' @return grid a list of factors, one for each donor, corresponding to the binned events
# #' @return xtab a table with number of events in each bin, one for each donor. 
# make_bins <- function(DONORS_reduced, sample_fac, k, l) {
#   DONORS_grid <- make_grid_index(DONORS_reduced[, 1:k],  l)
#   DONORS_grid <- factor(DONORS_grid, levels = unique(DONORS_grid))
#   DONORS_grid <- split(DONORS_grid, factor(sample_fac))
#   DONORS_crosstab <- mapply(table, DONORS_grid, SIMPLIFY = FALSE)
#   return(list(grid = DONORS_grid, xtab = DONORS_crosstab))
#   
# }
# 
# # ##############
# # not currently needed
# # ##############
# # bins a dimension-reduced pooled dataset. Can be with PCA or SMA. Then rebins using SMA.
# # DONORS_reduced : pooled dataset transformed by PCA or SMA
# # D : raw dataset, linear
# # sample_fac : index of which sample each event originated from, length(sample_fac) == nrow(D)
# # k1 : the number of components to include in analysis
# # l1 : number of bins along each component
# # k2 : the number of components for re-binning
# # l2_seq : either a vector of number of bins, or a single value 
# # make_bins_rebin <- function(DONORS_reduced, D, sample_fac, k1, l1, k2, l2_seq) {
# #   #k2 <- 1;
# #   #k1<-k1;k2<-1;l1<-l1;l2_seq<-l2_seq;
# #   #DONORS_reduced <- D_transform
# #   grid <- make_grid_index(DONORS_reduced[, 1:k1], l1)
# #   o <- order(grid)
# #   tmp <- D[o, ]
# #   ABC <- 
# #     lapply(l2_seq, function(.L2) {
# #       comp1 <- by(data.frame(tmp), grid[o], function(x) if (nrow(x) > 50) {comp_svd(x)$u[, 1:k2]} else {nrow(x)})
# #       nr <- by(data.frame(tmp), grid[o], nrow)
# #       rebin <- mapply(function(x, y) if (y > 50) {make_grid_index(x, .L2)} else (factor(rep(1, x))), comp1, nr, SIMPLIFY = FALSE)
# #       bins <- data.frame(bin1 = factor(grid[o], levels = unique(grid)), bin2 = factor(unlist(rebin), levels = unique(unlist(rebin))))
# #       c <- apply(bins, 1, function(x) paste(x, sep = '', collapse = '!'))
# #       c <- factor(c, levels = unique(c))
# #       bins$combined <- c
# #       bins <- bins[order(o), ]
# #       bins <- split(bins, sample_fac)
# #       xtab <- mapply(function(x) c(table(x$combined)), bins, SIMPLIFY = FALSE)
# #       return(list(grid = bins, xtab = xtab))
# #     }
# #     )
# #   
# #   return(list( k1 = k1, l1 = l1, k2 = k2, l2 = l2_seq, output = ABC))  
# # }
# 
# # bins a dimension-reduced pooled dataset with PCA. Then rebins using PCA.
# # DONORS_reduced :compensated data
# # sample_fac : index of which sample each event originated from, length(sample_fac) == nrow(D)
# # k1 : the number of components to include in analysis
# # l1 : number of bins along each component
# # k2 : the number of components for re-binning
# # l2_seq : either a vector of number of bins, or a single value 
# 
# 
# # ############
# # not currently needed
# # ############
# # make_bins_rebin_pca <- function(D, sample_fac, k1, l1, k2, l2_seq) {
# #   DONORS_reduced <- D %*% eigen(cov(D))$vectors
# #   grid <- make_grid_index(DONORS_reduced[, 1:k1], l1)
# #   o <- order(grid)
# #   tmp <- D[o, ]
# #   ABC <- 
# #     lapply(l2_seq, function(.L2) {
# #       comp1 <- by(data.frame(tmp), grid[o], function(x) if (nrow(x) > 50) {(as.matrix(x) %*% eigen(cov(x))$vectors)[, 1:k2]} else {nrow(x)})
# #       nr <- by(data.frame(tmp), grid[o], nrow)
# #       rebin <- mapply(function(x, y) if (y > 50) {make_grid_index(x, .L2)} else (factor(rep(1, x))), comp1, nr, SIMPLIFY = FALSE)
# #       bins <- data.frame(bin1 = factor(grid[o], levels = unique(grid)), bin2 = factor(unlist(rebin), levels = unique(unlist(rebin))))
# #       c <- apply(bins, 1, function(x) paste(x, sep = '', collapse = '!'))
# #       c <- factor(c, levels = unique(c))
# #       bins$combined <- c
# #       bins <- bins[order(o), ]
# #       bins <- split(bins, sample_fac)
# #       xtab <- mapply(function(x) c(table(x$combined)), bins, SIMPLIFY = FALSE)
# #       return(list(grid = bins, xtab = xtab))
# #     }
# #     )
# #   
# #   return(list( k1 = k1, l1 = l1, k2 = k2, l2 = l2_seq, output = ABC))  
# # }
# # 
# # ################
# 
# # ############
# # not currently needed
# # ############
# # model selection using SVM
# 
# 
# # C : class labels, one for each sample (fcs file). THis should be a vector of characters, not a factor
# # xtab : output of make_grid. list of count vectors, one for each sample
# # r : number of test/training rounds
# # w : proportion of zero counts (over samples) allowed in a bin
# 
# # make_svm_donor <- function(xtab, r, w, C) {
# #   categories <- unique(C)
# #   categories_index <- lapply(categories, function(.c) which(C == .c))
# #   
# #   R <- lapply(1:r, function(.dummy) {
# #     #xtab <- XT
# #     samples <- lapply(categories_index, function(.c) { s <- sample(length(.c), floor(length(.c)/2))})
# #     train <- unlist(mapply(function(x, y) x[y], categories_index, samples))
# #     test <- unlist(mapply(function(x, y) x[-y], categories_index, samples))
# #     tab1 <- do.call(rbind, xtab[train])
# #     tab2 <- do.call(rbind, xtab[test])
# #     empty <- apply(tab1, 2, function(x) length(which(x == 0)))
# #     tab1 <- tab1 + rpois(length(tab1), 1)
# #     tab2 <- tab2 + rpois(length(tab2), 1)
# #     tab <- table(C[test]) 
# #     nc <- ncol(tab1)
# #     g <- c(0.1/nc, 0.5/nc, 1/nc, 2/nc, 10/nc)
# #     loop <- lapply(w, function(.w) {
# #       lapply(g, function(.g) {
# #         ssvm <- svm(tab1[, empty < (nrow(tab1) * .w)], factor(C[train]), C = 1, gamma = .g)
# #         pred <- predict(ssvm, tab2[, empty < (nrow(tab1) * .w)])
# #         match <- cbind(levels(pred)[pred], C[test])
# #         match <- apply(match, 1, function(x) x[1] == x[2])
# #         match <- as.vector(by(match, C[test], sum)/table(C[test]))
# #         mean(match)
# #       }
# #       )
# #     }
# #     )
# #     loop1 <-  mapply(max, mapply(unlist, loop, SIMPLIFY = FALSE))
# #     return(loop1) 
# #   } 
# #   )
# #   return(R)
# # }
# 
# ###############
# 
# # ##########
# # old version
# # ##########
# # after binning, find the Aitchison distance (compositional distance) between all samples. Find a non-parametric F-statistic with respect to this inter-sample distance matrix and the categories of the sample
# # xtab : output of make_grid. list of count vectors, one for each sample
# # C : class labels, one for each sample (fcs file)
# # w : proportion of zero counts (over samples) allowed in a bin
# 
# # make_adist <- function(xtab, w, C) {
# #   #tab <- bp_model$xtab
# #   tab <- do.call(rbind, xtab)
# #   empty <- apply(tab, 2, function(x) length(which(x == 0)))
# #   which_empty <- which(empty < (length(nrow(tab)) * w))
# #   
# #   donor_coords <- tab
# #   donor_coords <- cbind(rowSums(donor_coords[, -which_empty]), donor_coords[, which_empty])
# #   donor_coords <- t(apply(donor_coords, 1, function(x) x/sum(x)))
# #   nr <- nrow(donor_coords)
# #   ad <- donor_coords + 0.0001
# #   ad <- apply(ad, 2, function(x) x/sum(x))
# #   ad <- lapply(1:(nr - 1), function(i) sapply((i+1):nr, function(j) aDist(ad[i, ], ad[j, ])))
# #   tmp <- matrix(0, nrow = nr - 1, ncol = nr - 1)
# #   for (i in 1:(nr - 1))
# #     tmp[i:(nr - 1), i] <- ad[[i]]
# #   ad <- rbind(0, cbind(tmp, 0))
# #   ad <- ad + t(ad)
# #   Ftest <- permutest(betadisper(as.dist(ad), bias.adjust = TRUE, group = factor(C)))$tab[1, c(4, 6)]
# #   return(list(distmtx = ad, Ftest = Ftest))
# # }
# 
# 
# 
# 
# ############
# 
# #' After binning the data, make a table with the nearly empty bins combined into a 'leftover bin'. often dramatically reduces size of the table without much information loss. 
# #' @param data_binned output of make_bins
# #' @param w the 'vacancy rate', a number between 0 and 1. The final table will have no column (bin) with a greater proportion than w of donors with zero events in that bin. 
# #' @param C The row names - donor sample type
# #' @return The reduced data matrix. Each row (donor) sums to one. The first column is the number of 'leftover' events, i.e. the events that don't fit into the bins that are common to most samples. 
# get_xtab <- function(data_binned, w, C) {
#   #data_binned <- bins_k2[[1]]
#   #w <- 0.2
#   #C <- sample_names1
#   tab <- data_binned$xtab
#   tab <- do.call(rbind, tab)
#   empty <- apply(tab, 2, function(x) length(which(x == 0)))
#   which_empty <- which(empty < (length(C) * w))
#   
#   donor_coords <- tab
#   donor_coords <- donor_coords[, which_empty]
#   if (ncol(donor_coords) == (ncol(tab) - 1)) {
#     donor_coords <- cbind(tab[, -which_empty], donor_coords)
#   }
#   if (ncol(donor_coords) < (ncol(tab) - 1)) {
#     donor_coords <- cbind(rowSums(tab[, -which_empty]), donor_coords)
#   }
#   donor_coords <- t(apply(donor_coords, 1, function(x) x/sum(x)))
#   rownames(donor_coords) <- C
#   return(donor_coords)
# }
# 
# #' Calculates the Aitchison (compositional) distance between samples.
# #' @param table of counts in each bin for each donor (row). Rows must sum to one. 
# #' @return a distance matrix between donors. 
# aDist_mtx <- function(table) { # table > 0
#   gm <- apply(table, 1, function(x) exp(mean(log(x))))
#   table_transformed <- log(sweep(table, 1, gm, `/`))
#   ad <- dist(table_transformed)
# }
# 
# 
# 
# 
# #############
# 
# #' puts a data matrix (compensated, uncompensated) into a format for lineplots
# #' ncol(X) == length(markers)
# #' @param X data matrix with observations/events as rows, features/channels as columns
# #' @param markers The column/channel names, e.g. the biomarkers used
# #' @return a data frame for making lineplots/parallel coordinate plots
# make_lineplot_data <- function(X, markers) {
#   spec = t(X)[1:(nrow(X) * ncol(X))]
#   ch = factor(rep(markers, nrow(X)), levels = markers, ordered = TRUE)
#   bin <- factor(rep(1:nrow(X), rep(ncol(X), nrow(X))))
#   df <- data.frame(spec = spec, ch = ch, bin = bin)
#   return(df)
# }
# 
# #' makes lineplots
# #' @param bins  a list of data matrices, corresponding to bins
# #' @param dir  the directory where the plots should be made
# #' @param markers  vector of strings indicating the markers
# #' @param threshold  the minimum sized bin that should be plotted
# #' @param y1, y2  vertical plot limits
# #' @return a plot in the desired direction (png format)
# do_plots <- function(bins, dir, markers, threshold, y1, y2) {
#   for (i in 1:length(bins[[1]])) {
#     tmp <- bins[[1]][[i]]
#     colnames(tmp) <- markers
#     if (nrow(tmp) > threshold) {
#       if (nrow(tmp) > 10000) {line_data <- make_lineplot_data(tmp[sample(nrow(tmp), 10000), ], markers)} else {line_data <- make_lineplot_data(tmp, markers)}
#       size <- 10/log(nrow(tmp)/length(markers), 1.05) 
#       fn <- paste(dir, 'bin_', bins[[2]][[i]], '.png', sep = '')
#       png(fn)
#       print(ggplot(line_data, aes(x = ch, y = spec, group = bin)) + geom_path(colour = rgb(0.2, 0.2, 0.2, 0.4), size = size) + ylim(y1, y2) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)))
#       dev.off()
#     }
#   }
# }
# 
# # ###########
# # not currently needed
# # ###########
# # makes lineplots
# # comp1 : a list of 1st components of each bin. can be either pca or SMA
# # bins : a list of data matrices, corresponding to bins
# # dir : the directory where the plots should be made
# # markers : vector of strings indicating the markers
# # threshold : the minimum sized bin that should be plotted
# # y1, y2 : vertical plot limits
# # no_panels : number of re-bins
# # do_plots_rebin <-
# #   function(comp1, bins, dir, markers, threshold, y1, y2, no_panels) {
# #     
# #     for (i in 1:length(bins[[1]])) {
# #       tmp <- bins[[1]][[i]]     
# #       if (nrow(tmp) > threshold) {
# #         if (nrow(tmp) > 10000) {
# #           s <- sample(nrow(tmp), 10000)
# #           line_data <- make_lineplot_data(tmp[s, ], markers)
# #           group <- comp1[[i]][s]
# #           group <- group - min(group)
# #           group <- group/max(group)
# #           group <- floor(group * no_panels)
# #           group[group == max(group)] <- max(group - 1)
# #           group <- rep(group, rep(length(markers), length(group)))
# #           line_data$group <- factor(group)
# #         } else {
# #           line_data <- make_lineplot_data(tmp, markers)
# #           group <- comp1[[i]]
# #           group <- group - min(group)
# #           group <- group/max(group)
# #           group <- floor(group * no_panels)
# #           group[group == max(group)] <- max(group - 1)
# #           group <- rep(group, rep(length(markers), length(group)))
# #           line_data$group <- factor(group)
# #         }
# #         size <- 10/log(nrow(tmp)/length(markers), 1.05) 
# #         fn <- paste(dir, 'bin_', bins[[2]][[i]], '.png', sep = '')
# #         png(fn, height = 700, width = 700)
# #         print(ggplot(line_data, aes(x = ch, y = spec, group = bin)) + geom_path(colour = mygrey, size = size) + facet_wrap(~ group) + ylim(y1, y2) + ggtitle(paste('Bin', i)) + theme(plot.title = element_text(size=50), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))) 
# #         dev.off()
# #       }
# #     }
# #   }
# # 
# # 
# # # splits a datamatrix D into bins with respect to the output of make_bins. 
# # # data_binned : output of make_bins
# # # w : proportion of zero counts (over samples) allowed in a bin
# # # D : the matrix to split. Probably compensated, so you can see what the markers are doing
# # # sample_fac : index of which sample each event originated from, length(sample_fac) == nrow(D)
# # get_bins_lineplot <- function(data_binned, w, D, sample_fac, rebin) {
# # #   data_binned <- model_rebin[[5]][[1]]
# # #   w <- 0.8
# # #   sample_fac <- D$sample_fac
# # #   D <- D$D_comp[, markers_index]
# #   tab <- data_binned$xtab
# #   tab <- do.call(rbind, tab)
# #   binnames <- colnames(tab)
# #   empty <- apply(tab, 2, function(x) length(which(x == 0)))
# #   which_empty <- which(empty < (nrow(tab) * w))
# #   
# #   grid <- data_binned$grid
# #   x <- split(data.frame(D), factor(sample_fac))
# #   if (!rebin) {
# #     data_bins <- mapply(function(x, y) split(x, y), x, grid, SIMPLIFY = FALSE)
# #   } else {
# #     data_bins <- mapply(function(x, y) split(x, y$combined), x, grid, SIMPLIFY = FALSE) 
# #   }
# #   
# #   data_bins <- mapply(function(x) x[which_empty], data_bins, SIMPLIFY = FALSE)
# #   
# #   combine_bins <- lapply(1:length(which_empty), function(j) do.call(rbind, lapply(1:length(data_bins), function(i) data_bins[[i]][[j]])))
# #   return(list(combine_bins, binnames[which_empty]))
# #   
# # }
# # 
# # # calculates a new SMA for each bin, in order to re-bin 
# # # data_binned : output of make_bins
# # # w : proportion of zero counts (over samples) allowed in a bin
# # # D : the matrix to split and from which the first component is calculated. it is probably the raw linear data the first bins were calculated, D > 0
# # # sample_fac : index of which sample each event originated from, length(sample_fac) == nrow(D)
# # get_comp1 <- function(data_binned, w, D, sample_fac) {
# #   
# #   tab <- data_binned$xtab
# #   tab <- do.call(rbind, tab)
# #   empty <- apply(tab, 2, function(x) length(which(x == 0)))
# #   which_empty <- which(empty < (nrow(tab) * w))
# #   
# #   grid <- data_binned$grid
# #   x <- split(data.frame(D), factor(sample_fac))
# #   data_bins <- mapply(function(x, y) split(x, y), x, grid, SIMPLIFY = FALSE)
# #   data_bins <- mapply(function(x) x[which_empty], data_bins, SIMPLIFY = FALSE)
# #   combine_bins <- lapply(1:length(which_empty), function(j) do.call(rbind, lapply(1:length(data_bins), function(i) data_bins[[i]][[j]])))
# #   comp1 <- mapply(function(x) comp_svd(x)$u[, 1], combine_bins, SIMPLIFY = FALSE)
# #   
# #   return(comp1)
# #   
# # }
# # 
# # get_comp1_pca <- function(data_binned, w, D, sample_fac) {
# #   
# #   tab <- data_binned$xtab
# #   tab <- do.call(rbind, tab)
# #   empty <- apply(tab, 2, function(x) length(which(x == 0)))
# #   which_empty <- which(empty < (nrow(tab) * w))
# #   
# #   grid <- data_binned$grid
# #   x <- split(data.frame(D), factor(sample_fac))
# #   data_bins <- mapply(function(x, y) split(x, y), x, grid, SIMPLIFY = FALSE)
# #   data_bins <- mapply(function(x) x[which_empty], data_bins, SIMPLIFY = FALSE)
# #   combine_bins <- lapply(1:length(which_empty), function(j) do.call(rbind, lapply(1:length(data_bins), function(i) data_bins[[i]][[j]])))
# #   comp1 <- mapply(function(x) as.matrix(log(x, 10)) %*% eigen(cov(log(x, 10)))$vectors[, 1], combine_bins, SIMPLIFY = FALSE)
# #   
# #   return(comp1)
# #   
# # }
# 
# 
# 
# ################
# 
# # ################
# # not currently needed
# # ################
# 
# # finds the mean direction, by first normalising all vectors (rows) to have length one
# # X : data matrix, rows are cells
# # angular_mean <- function(X) {
# #   X <- t(apply(X, 1, function(x) x/sqrt(sum(x^2))))
# #   X <- colSums(X) / nrow(X)
# #   return(X)
# # }
# # 
# # # generates data for plotting a layout of the bins, with  size of each bin according to the cell count in a category and colour according to marker expression
# # # BP_output : make_grid output
# # # D : the raw data (linear)
# # # D_comp : the compensated matrix
# # # minimum cell count in bin to include in plot (pooled over categories)
# # # C : sample categories
# # 
# # 
# # bins_MDS <- function(data_binned, D, D_comp, thres, C, sample_fac) {
# #   
# #   # find angular mean of each bin and make distance matrix, MDS
# #   D <- log(D, 10)
# #   gr <- do.call(rbind, data_binned)$combined
# #   bins <- split(data.frame(D), gr)
# #   nr <- mapply(nrow, bins)
# #   bins <- bins[nr > thres] 
# #   
# #   bins_center <- mapply(angular_mean, bins, SIMPLIFY = FALSE)
# #   bins_center <- do.call(rbind, mapply(function(x) x/sqrt(sum(x^2)), bins_center, SIMPLIFY = FALSE))
# #   xprod <- bins_center %*% t(bins_center)
# #   xprod <- acos(xprod)
# #   diag(xprod) <- 0
# #   mds_coords <- cmdscale(xprod)
# #   
# #   # this will be 'col' parameter
# #   bins_comp <- split(data.frame(D_comp), gr)
# #   bins_comp <- bins_comp[nr > thres] 
# #   
# #   bins_center_comp <- mapply(colMeans, bins_comp, SIMPLIFY = FALSE)
# #   bins_center_comp <- do.call(rbind, mapply(function(x) x/sqrt(sum(x^2)), bins_center_comp, SIMPLIFY = FALSE))
# #   
# #   mypal <- colorRampPalette(c('gray4', 'mediumturquoise'))(100)
# #   Q <- data.frame(apply(bins_center_comp, 2, function(x) quantile(x, seq(0, 1, length.out = 100))))
# #   col <- mapply(function(x, Q) sapply(1:length(x), function(j) length(which(Q < x[j]))) + 1, data.frame(bins_center_comp), Q, SIMPLIFY = FALSE)
# #   col <- mapply(function(x) mypal[x], col, SIMPLIFY = FALSE)
# #   names(col) <- markers
# #   
# #   # bins in each donor
# #   gr1 <- data_binned
# #   D1 <- split(data.frame(D), factor(sample_fac))
# #   indiv_bins <- mapply(function(x, y) split(x, y$combined), D1, gr1, SIMPLIFY = FALSE)
# #   
# #   # this will be 'cex' parameter
# #   sample_bincounts <- lapply(unique(C), function(.x) {
# #     B <- indiv_bins[C == .x]
# #     B <- lapply(1:length(which(nr > thres)), function(i) do.call(rbind, mapply(function(x) x[nr > thres][[i]], B, SIMPLIFY = FALSE)))
# #     B <- mapply(nrow, B)/length(which(C == .x))
# #     return(B)
# #   }
# #   )
# #   Q <- quantile(unlist(sample_bincounts), seq(0, 1, length.out = 100))
# #   cex <- mapply(function(x) sapply(1:length(x), function(j) length(which(Q < x[j]))) + 1, sample_bincounts, SIMPLIFY = FALSE)
# #   C <- seq(0.7, 1.3, length.out = 100)^2
# #   cex <- mapply(function(x) C[x], cex, SIMPLIFY = FALSE)
# #   
# #   return(list(mds_coords = mds_coords, col = col, cex = cex))
# # }
# # 
# # ######################