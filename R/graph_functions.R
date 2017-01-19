#' wrapper function, takes a data frame of marker pairs and constructs a graph
#' @param pairs dataframe with two columns, containing the marker pairs of interest
#' @param markers a vector of strings describing all markers. It is not necessary that every element that appears in markers also appears in pairs.
#' @return a graph
#' @export
make_subgraph <- function(pairs, markers) {
  pairs <- graph_from_data_frame(pairs, directed = FALSE, vertices = markers)
  return(pairs)
}

#'  takes a subgraph, returns a list of triangles
#'  @param subgraph output from make_subgraph
#'  @return a list of triplets, describing the triangles in the graph
#'  @export
get_triangles <- function(subgraph) {
  tri <- triangles(subgraph)
  tri <- split(tri, rep(1:(length(tri)/3), rep(3, length(tri)/3)))
  tri <- mapply(sort, tri, SIMPLIFY = FALSE)
  return(tri)
}

#' checks equality of two triangles
#' @param A, B both vectors of length three. Output from get_triangles. Assumes A and B have been ordered, i.e. doesn't check permutations.
#' @return logical
#' @export
tri_equal <- function(A, B) {
  A[1] == B[1] & A[2] == B[2] & A[3] == B[3]
}


#' function to find overlapping triangles between subgraphs
#' @param subgraphs a list of graphs, each is output of make_subgraphs. Length of list is at least two. List must be names.
#' @param triangles A list of triangles, each is output of get_triangles. Order of subgraphs and triangles must be the same.
#' @return a named list of triangles that overlap every pair of graphs.
#' @export
bridge_triangles <- function(subgraphs, triangles) {
  bridge_tri <-
  lapply(1:(length(subgraphs)-1), function(I) lapply((I+1):length(subgraphs), function(J) {

    merge_graphs <- subgraphs[[I]] + subgraphs[[J]]
    all_tri <- get_triangles(merge_graphs)

    z1 <- lapply(1:length(triangles[[I]]), function(j) sapply(1:length(all_tri), function(i) tri_equal(all_tri[[i]], triangles[[I]][[j]])))
    z1 <- Reduce(`+`, z1)
    z2 <- lapply(1:length(triangles[[J]]), function(j) sapply(1:length(all_tri), function(i) tri_equal(all_tri[[i]], triangles[[J]][[j]])))
    z2 <- Reduce(`+`, z2)
    z3 <- as.numeric(z1) + as.numeric(z2)

    bt <- all_tri[which(z3 == 0)]
    bt <- mapply(as.numeric, bt)
    return(bt)
  }
  )
  )

  bridge_tri <- unlist(bridge_tri, recursive = FALSE)
  newnames <- names(subgraphs)
  newnames <- lapply(1:(length(subgraphs)-1), function(I) lapply((I+1):length(subgraphs), function(J) paste(newnames[I], newnames[J], sep = '_')))
  newnames <- unlist(newnames)
  names(bridge_tri) <- newnames
  return(bridge_tri)
}



