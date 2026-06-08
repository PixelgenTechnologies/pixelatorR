#' Coarsened pmds
#'
#' This function implements an algorithm to compute a layout on a coarsened version
#' of the graph. By coarsening the graph, it is possible to process large graphs
#' faster and using less memory.
#'
#' @section Algorithm:
#' 1. Cluster the graph using the Leiden method and create a coarsened version of
#' the graph where each community is represented as a single node. Edge weights represent
#' the number of connections between communities. The resolution parameter controls the
#' granularity of the clustering, with higher values leading to more communities. Fewer
#' communities will lead to faster computation but potentially less accurate layouts.
#' 2. Compute a layout for the coarsened graph using PMDS. The `weight_edges_by` parameter
#' controls how edges are weighted in the coarsened graph. "tp" uses bidirectional transition
#' probabilities as described in `layout_with_weighted_pmds`, while "crossing_edges" uses 1
#' divided by the number of edges crossing between communities.
#' 3. Each node in the original graph is assigned the coordinates of its corresponding
#' community in the coarsened graph. Then, the nodes are iteratively "wiggled" by replacing
#' their coordinates with a weighted average of their neighbors' coordinates + some random
#' jitter (`jitter_sd`).
#'
#' @param g A `tbl_graph` or `igraph` object.
#' @param dim Number of dimensions for the layout (default is 3).
#' @param resolution Resolution parameter for the community detection step.
#' Higher values lead to more communities.
#' @param pivots Number of pivots to use for the PMDS layout step.
#' @param n_iter Number of iterations for the smoothing step.
#' @param jitter_sd Standard deviation of the jitter added during smoothing.
#' @param weight_edges_by How to weight edges in the coarsened graph for the
#' PMDS layout. "tp" uses bidirectional transition probabilities for 5-step
#' random walks using the `layout_with_weighted_pmds` function. "crossing_edges"
#' uses 1 divided by the number of edges crossing between communities (from the
#' coarsening step) as weights.
#' @param seed Random seed for reproducibility.
#' @param verbose Whether to print messages about the coarsening process.
#'
#' @examples
#' library(dplyr)
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
#'   LoadCellGraphs(cells = colnames(.)[4], verbose = FALSE)
#'
#' cg <- CellGraphs(se)[[4]]
#' g <- cg@cellgraph
#'
#' xyz <- layout_with_coarsened_pmds(g) %>%
#'   as_tibble(.name_repair = ~ c("x", "y", "z"))
#'
#' plotly::plot_ly(
#'   xyz,
#'   x = ~x, y = ~y, z = ~z,
#'   mode = "markers",
#'   type = "scatter3d",
#'   marker = list(size = 1)
#' )
#'
#' @return A matrix of layout coordinates for the original graph.
#'
#' @export
#'
layout_with_coarsened_pmds <- function(
  g,
  dim = 3,
  resolution = 1,
  pivots = 200,
  n_iter = 10,
  jitter_sd = 1e-2,
  weight_edges_by = c("tp", "crossing_edges"),
  seed = 123,
  verbose = FALSE
) {
  pixelatorR:::assert_class(g, c("tbl_graph", "igraph"))
  pixelatorR:::assert_single_value(dim, "integer")
  pixelatorR:::assert_within_limits(dim, c(2, 3))
  pixelatorR:::assert_single_value(resolution, "numeric")
  pixelatorR:::assert_within_limits(resolution, c(0.01, 10))
  pixelatorR:::assert_single_value(pivots, "integer")
  pixelatorR:::assert_within_limits(pivots, c(10, min(1000, length(g))))
  pixelatorR:::assert_single_value(n_iter, "integer")
  pixelatorR:::assert_within_limits(n_iter, c(1, 100))
  pixelatorR:::assert_single_value(jitter_sd, "numeric")
  pixelatorR:::assert_within_limits(jitter_sd, c(1e-3, 0.1))
  weight_edges_by <- match.arg(weight_edges_by, c("tp", "crossing_edges"))
  pixelatorR:::assert_single_value(seed, "integer")
  pixelatorR:::assert_single_value(verbose, "bool")
  set.seed(seed)

  # Normalize resolution parameter to fit PNA graphs
  # Here we use the number of edges to scale the resolution parameter
  res <- resolution * igraph::gsize(g) / 1000

  # Run Leiden community detection
  cl <- igraph::cluster_leiden(g, resolution = res, objective_function = "modularity")$membership
  if (verbose) {
    cli::cli_alert_info("Graph size: {length(g)} nodes, {igraph::gsize(g)} edges.")
    cli::cli_alert_info("Detected {length(unique(cl))} communities.")
    cli::cli_alert_info("Average community size: {round(mean(table(cl)))} nodes.")
  }

  # Fetch adjacency matrix
  A_orig <- igraph::as_adjacency_matrix(g)

  # Create a cluster-level adjacency matrix by summing connections between clusters
  mm <- Matrix::sparse.model.matrix(~ 0 + factor(cl)) %>% as("dgCMatrix")
  cl_counts <- Matrix::t(mm) %*% A_orig %*% mm

  # Remove internal cluster connections by setting the diagonal to zero
  Matrix::diag(cl_counts) <- 0
  cl_counts <- Matrix::drop0(cl_counts, tol = 0)

  if (weight_edges_by == "tp") {
    # When using transition probabilities, we need to set the edge weights to 1
    A <- cl_counts
    A@x <- rep(1, length(A@x))
    g_small <- igraph::graph_from_adjacency_matrix(A, mode = "upper") %>%
      as_tbl_graph(directed = FALSE)

    xyz <- layout_with_weighted_pmds(g_small, dim = dim, pivots = min(pivots, length(g_small)))
  }

  if (weight_edges_by == "crossing_edges") {
    # When using crossing edges, we set the edge weights to 1 divided by the number
    # of edges crossing between communities
    g_small <- igraph::graph_from_adjacency_matrix(cl_counts, mode = "upper", weighted = "weight") %>%
      as_tbl_graph(directed = FALSE)
    w <- 1 / (g_small %E>% pull(weight))
    xyz <- fast_pmds(g_small, dim = dim, pivots = min(pivots, length(g_small)), weights = w)
  }

  if (verbose) {
    cli::cli_alert_info("Coarsened graph size: {length(g_small)} nodes, {igraph::gsize(g_small)} edges.")
  }

  # Normalize the layout coordinates such that the median radius is 1
  xyz <- xyz %>%
    normalize_layout_coordinates() %>%
    as.matrix()

  # Create transition probability matrix weighted by within and between clusters
  Matrix::diag(A_orig) <- 1

  # Get the row and column index for every non-zero element
  rows <- A_orig@i + 1
  cols <- rep(seq_len(ncol(A_orig)), diff(A_orig@p))

  # Identify which of these non-zero elements connect within communities
  keep <- cl[rows] == cl[cols]

  # Filter the adjacency matrix to keep only within-community connections
  A_within <- A_orig
  A_within@x <- A_within@x * keep

  # Drop the now-zeroed entries to maintain sparsity
  A_within <- Matrix::drop0(A_within)

  # Calculate the between-community adjacency matrix by subtracting the
  # within-community matrix from the original adjacency matrix
  A_between <- Matrix::drop0(A_orig - A_within)

  # Normalize the within and between adjacency matrices to create transition probabilities
  # Use pmax to avoid division by zero for rows that have no connections
  P_between <- A_between / pmax(Matrix::rowSums(A_between), 1)
  P_within <- (A_within / pmax(Matrix::rowSums(A_within), 1))

  # Combine the within and between transition probabilities and renormalize to ensure rows sum to 1
  P <- P_within + P_between
  P <- P / Matrix::rowSums(P)

  # Extrapolate points to whole cell
  # Repeat xyz such that each cl$membership is matched
  xyz_full <- xyz[cl, ]

  for (i in 1:n_iter) {
    # Add noise
    xyz_full <- xyz_full + matrix(rnorm(prod(dim(xyz_full)), sd = jitter_sd), nrow = nrow(xyz_full))
    # Apply smoothing
    xyz_full <- (P %*% xyz_full)
  }
  xyz_full <- as.matrix(xyz_full)

  # Renormalize coordinates after smoothing
  xyz_full <- xyz_full %>%
    normalize_layout_coordinates() %>%
    as.matrix()
  colnames(xyz_full) <- c("x", "y", "z")[seq_len(dim)]

  return(xyz_full)
}
