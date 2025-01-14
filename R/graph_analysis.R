#' Calculate antibody counts per A-node
#'
#' Computes and returns a data frame of antibody counts per node (vertex) of
#' the A node graph given a component edge list as input. The parameter \code{k}
#' allows to include neighbors (of each node) when computing the counts.
#' \code{k} defines the number of levels when searching neighbors.
#'
#' If \code{k > 0}, \code{\link{edgelist_to_simple_Anode_graph}} is used
#' to calculate an A-node projected graph.
#'
#' @param component_edge_list An object of class \code{tbl_df}
#' @param k Number of neighbors to include
#'
#' @return A matrix with node marker counts
#'
#' @export
#'
node_markers_counts <- function(
  component_edge_list,
  k = 0
) {
  # Check input parameters
  assert_non_empty_object(component_edge_list, "tbl_df")
  assert_single_value(k, "integer")
  assert_x_in_y(c("upia", "upib", "marker"), colnames(component_edge_list))

  component_counts <-
    component_edge_list %>%
    group_by(upia, marker) %>%
    count() %>%
    pivot_wider(id_cols = "upia", names_from = "marker", values_from = "n", values_fill = 0) %>%
    column_to_rownames("upia")

  if (k > 0) {
    # Simplify graph
    simple_graph <-
      component_edge_list %>%
      select(upia, upib) %>%
      distinct() %>%
      edgelist_to_simple_Anode_graph(verbose = FALSE)

    # fetch adjacency matrix
    adj_mat <-
      connect(simple_graph[[1]], order = k) %>%
      as_adjacency_matrix()
    diag(adj_mat) <- 1

    # Filter component_counts
    component_counts <- component_counts[colnames(adj_mat), ]

    return(as.matrix(adj_mat %*% as(as.matrix(component_counts), "dgCMatrix")))
  } else {
    return(as.matrix(component_counts))
  }
}
