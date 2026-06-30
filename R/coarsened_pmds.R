#' Coarsened pMDS
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
#' @param leiden_iterations Number of iterations for the Leiden community detection algorithm.
#' @param leiden_weighted Whether to use weighted edges for the Leiden community detection.
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
  n_iter = 20,
  jitter_sd = 1e-2,
  weight_edges_by = c("crossing_edges", "tp"),
  leiden_iterations = 2,
  leiden_weighted = FALSE,
  seed = 123,
  verbose = FALSE
) {
  assert_class(g, c("tbl_graph", "igraph"))
  assert_single_value(dim, "integer")
  assert_within_limits(dim, c(2, 3))
  assert_single_value(resolution, "numeric")
  assert_within_limits(resolution, c(0.1, 10))
  assert_single_value(pivots, "integer")
  assert_within_limits(pivots, c(10, min(1000, length(g))))
  assert_single_value(n_iter, "integer")
  assert_within_limits(n_iter, c(1, 100))
  assert_single_value(jitter_sd, "numeric")
  assert_within_limits(jitter_sd, c(1e-3, 0.1))
  weight_edges_by <- match.arg(weight_edges_by, c("crossing_edges", "tp"))
  assert_single_value(leiden_iterations, "integer")
  assert_within_limits(leiden_iterations, c(1, 10))
  assert_single_value(leiden_weighted, "bool")
  assert_single_value(seed, "integer")
  assert_single_value(verbose, "bool")
  set.seed(seed)

  # Normalize resolution parameter to fit PNA graphs
  res <- resolution * igraph::gsize(g) / 2000

  # Run Leiden community detection
  if (leiden_weighted) {
    g <- g %>% prob_distance_weights(k = 1, min_weight = 0)
    ew <- -log10(E(g)$bi_prob)
    ew <- ew / mean(ew)
  } else {
    ew <- NULL
  }
  cl <- igraph::cluster_leiden(
    g,
    resolution = res,
    objective_function = "modularity",
    n_iterations = leiden_iterations,
    weights = ew
  )$membership
  if (verbose) {
    cli::cli_alert_info("Graph size: {length(g)} nodes, {igraph::gsize(g)} edges.")
    cli::cli_alert_info("Detected {length(unique(cl))} communities.")
    cli::cli_alert_info("Average community size: {round(mean(table(cl)))} nodes.")
  }

  # Fetch adjacency matrix
  A_orig <- igraph::as_adjacency_matrix(g)

  # Create a cluster-level adjacency matrix cleanly
  mm <- Matrix::sparseMatrix(
    i = seq_along(cl), 
    j = cl, 
    x = 1, 
    dims = c(length(cl), length(unique(cl))),
    giveCsparse = TRUE
  )
  cl_counts <- Matrix::t(mm) %*% A_orig %*% mm

  # Remove internal cluster connections by setting the diagonal to zero
  Matrix::diag(cl_counts) <- 0
  cl_counts <- Matrix::drop0(cl_counts, tol = 0)

  # Drop tidygraph cast: Keep as raw igraph for speed
  g_small <- igraph::graph_from_adjacency_matrix(cl_counts, mode = "upper", weighted = TRUE)
  igraph::E(g_small)$inverse_weight <- 1 / igraph::E(g_small)$weight

  if (weight_edges_by == "tp") {
    xyz <- layout_with_weighted_pmds(g_small, dim = dim, pivots = min(pivots, length(g_small)), seed = seed)
  }

  if (weight_edges_by == "crossing_edges") {
    w <- igraph::E(g_small)$inverse_weight
    xyz <- fast_pmds(g_small, dim = dim, pivots = min(pivots, length(g_small)), weights = w)
  }

  if (verbose) {
    cli::cli_alert_info("Coarsened graph size: {length(g_small)} nodes, {igraph::gsize(g_small)} edges.")
  }

  # Normalize the layout coordinates
  xyz <- normalize_layout_coordinates(xyz, as_tibble = FALSE)

  # Create transition probability matrix weighted by within and between clusters
  Matrix::diag(A_orig) <- 1

  # Compute row sums exactly ONCE up front
  row_sums_orig <- Matrix::rowSums(A_orig)
  # Guard against 0 values if any isolated nodes exist (unlikely with self loops, but safe)
  row_sums_orig <- pmax(row_sums_orig, 1) 

  # Extract triplets efficiently
  A_triplets <- Matrix::summary(A_orig)
  is_within <- cl[A_triplets$i] == cl[A_triplets$j]

  # Construct A_within and A_between directly
  A_within <- Matrix::sparseMatrix(
    i = A_triplets$i[is_within],
    j = A_triplets$j[is_within],
    x = A_triplets$x[is_within],
    dims = dim(A_orig)
  )

  A_between <- Matrix::sparseMatrix(
    i = A_triplets$i[!is_within],
    j = A_triplets$j[!is_within],
    x = A_triplets$x[!is_within],
    dims = dim(A_orig)
  )

  # Highly optimized fast diagonal multiplication for transition matrix P
  # Avoids calculating rowSums 3 times and prevents dense vector allocations
  inv_D <- Matrix::Diagonal(x = 1 / row_sums_orig)
  P <- inv_D %*% (A_within + A_between)

  # Extrapolate points to whole cell
  xyz_full <- xyz[cl, , drop = FALSE]

  n_rows <- nrow(xyz_full)
  n_cols <- ncol(xyz_full)
  total_elements <- n_rows * n_cols

  # Pre-allocation strategy for smoothing
  for (i in seq_len(n_iter)) {
    # Inline vector addition bypasses dense matrix creation over loops
    xyz_full <- xyz_full + rnorm(total_elements, mean = 0, sd = jitter_sd)
    xyz_full <- P %*% xyz_full
  }
  
  # Ensure array typing remains raw matrix 
  if (!is.matrix(xyz_full)) xyz_full <- as.matrix(xyz_full)

  # Renormalize coordinates after smoothing
  xyz_full <- normalize_layout_coordinates(xyz_full, as_tibble = FALSE)
  colnames(xyz_full) <- c("x", "y", "z")[seq_len(dim)]

  return(xyz_full)
}
