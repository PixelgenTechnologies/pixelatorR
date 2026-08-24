#' Spectral Layout via Graph Laplacian
#'
#' Computes a graph layout using the spectral geometry of the graph Laplacian.
#' Nodes that are well-connected in the graph are placed close together, while
#' structurally distant nodes are separated. Layout coordinates are always
#' passed through \code{\link{normalize_layout_coordinates}} before being returned.
#'
#' @section eigen:
#' The layout coordinates are the eigenvectors corresponding to the
#' (\code{dim}) smallest non-trivial eigenvalues of the Laplacian or the
#' normalized Laplacian, also known as Laplacian Eigenmaps or the spectral
#' embedding of the graph.
#'
#' 1. Compute the sparse adjacency matrix \eqn{A}. When
#' `normalize_laplacian = FALSE`, form the unnormalized graph Laplacian
#' \eqn{L = D - A}, where \eqn{D} is the diagonal degree matrix. When
#' `normalize_laplacian = TRUE`, form the normalized graph Laplacian
#' \eqn{L = I - D^{-1/2}AD^{-1/2}}, where \eqn{D^{-1/2}} is the inverse-square-root
#' diagonal degree matrix. The Laplacian is positive semi-definite, so all
#' eigenvalues are non-negative.
#' 2. Use an iterative sparse eigensolver (\code{RSpectra::eigs_sym}) to compute
#' the \code{dim + 1} smallest-magnitude eigenvalues and their eigenvectors. For
#' a connected graph, the smallest eigenvalue is exactly 0 and its eigenvector
#' is the trivial constant vector.
#' 3. Discard the trivial eigenvector (eigenvalue \eqn{\approx 0}) and set the
#' remaining \code{dim} eigenvectors as layout coordinates. If
#' `normalize_laplacian = TRUE`, the coordinates are corrected using
#' \eqn{D^{-1/2}} (mapping back to the random-walk Laplacian embedding).
#'
#' @section svd:
#' The layout coordinates are the singular vectors corresponding to the
#' (\code{dim}) largest non-trivial singular values of the normalized
#' biadjacency matrix. This solver requires a bipartite graph
#' (`attr(g, "type") == "bipartite"`) with a `node_type` vertex attribute and
#' named vertices, and always uses the normalized Laplacian
#' (`normalize_laplacian = TRUE`).
#'
#' Spectral graph layouts require finding the eigenvectors corresponding to the
#' smallest non-zero eigenvalues of the Laplacian matrix. For bipartite graphs,
#' finding the bottom eigenvalues of the normalized Laplacian is equivalent to
#' finding the largest singular values of the normalized biadjacency matrix.
#' This allows replacing the eigenvalue solver with an SVD solver tailored for
#' rectangular data, which often converges faster.
#'
#' 1. Compute the sparse biadjacency matrix \eqn{B}, the diagonal degree
#' matrices \eqn{D_1} (row sums of \eqn{B}) and \eqn{D_2} (column sums of
#' \eqn{B}), then form the normalized biadjacency matrix
#' \eqn{B_{\mathrm{norm}} = D_1^{-1/2} B D_2^{-1/2}}.
#' 2. Use singular value decomposition (\code{irlba::irlba}) to compute the
#' \code{dim + 1} largest-magnitude singular values and their singular vectors
#' from \eqn{B_{\mathrm{norm}}}. For a connected graph, the largest singular
#' value is exactly 1.
#' 3. Discard the trivial singular vector (singular value \eqn{\approx 1}) and
#' return the remaining \code{dim} singular vectors as layout coordinates. Left
#' and right singular vectors are combined so that every node receives
#' coordinates. The coordinates are corrected using \eqn{D_1^{-1/2}} and
#' \eqn{D_2^{-1/2}} (mapping back to the random-walk Laplacian embedding).
#'
#' @param g A \code{tbl_graph} or \code{igraph} object.
#' @param dim Number of dimensions for the layout (default is 3). Must be 2 or 3.
#' @param normalize_laplacian Whether to use the normalized graph Laplacian.
#' Only supported for \code{solver = "eigen"}; the SVD solver always uses the
#' normalized Laplacian.
#' @param solver One of \code{"svd"} (default) or \code{"eigen"}. See the
#' sections above. \code{"svd"} requires a named bipartite graph.
#' @param jitter_sd Adds a small amount of jitter to the node coordinates.
#' Degree-1 nodes can have nearly identical solutions and therefore overlay each
#' other when plotting. Setting \code{jitter_sd = 1e-2} is typically sufficient
#' to displace these while having a negligible effect on higher-degree node
#' coordinates. By default, no jitter is applied.
#' @param seed Random seed for reproducibility.
#' @param verbose Whether to print messages about the graph and
#' eigenvalues/singular values.
#' @param ... Additional arguments passed to \code{RSpectra::eigs_sym} or
#' \code{irlba::irlba}.
#'
#' @examples
#' library(dplyr)
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
#'   LoadCellGraphs(cells = colnames(.)[4], verbose = FALSE)
#'
#' cg <- CellGraphs(se)[[4]]
#' g <- cg@cellgraph
#'
#' xyz <- layout_with_spectral(g) %>%
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
#' @return A numeric matrix of layout coordinates with \code{dim} columns named
#' \code{"x"}, \code{"y"}, and (when \code{dim = 3}) \code{"z"}, with one row
#' per node in \code{g}.
#'
#' @seealso [layout_with_weighted_pmds()], [layout_with_coarsened_pmds()],
#' [normalize_layout_coordinates()], [ComputeLayout()]
#'
#' @export
#'
layout_with_spectral <- function(
  g,
  dim = 3,
  normalize_laplacian = TRUE,
  solver = c("svd", "eigen"),
  jitter_sd = 0,
  seed = 123,
  verbose = FALSE,
  ...
) {
  assert_class(g, c("tbl_graph", "igraph"))
  assert_single_value(dim, "integer")
  if (!dim %in% c(2, 3)) {
    cli::cli_abort("{.var dim} must be either 2 or 3. Got {dim}")
  }
  assert_single_value(normalize_laplacian, "bool")
  solver <- match.arg(solver, c("svd", "eigen"))
  if (!normalize_laplacian && solver == "svd") {
    cli::cli_abort(
      "Using {.var normalize_laplacian = FALSE} is not supported for {.var solver = 'svd'}"
    )
  }
  if (solver == "svd") {
    graph_type <- attr(g, "type")
    if (is.null(graph_type)) {
      cli::cli_abort("The graph is missing a {.str type} attribute")
    }
    if (!(graph_type == "bipartite")) {
      cli::cli_abort("The graph must be bipartite when {.var solver = 'svd'}")
    }
    node_type_attr <- igraph::vertex_attr(g, "node_type")
    if (is.null(node_type_attr)) {
      cli::cli_abort("The graph is missing node attribute {.str node_type}")
    }
    node_types <- node_type_attr == "umi1"
  }
  assert_single_value(jitter_sd, "numeric")
  assert_within_limits(jitter_sd, c(0, 0.1))
  assert_single_value(seed, "integer")
  assert_single_value(verbose, "bool")
  if (solver == "svd" && is.null(igraph::vertex_attr(g, "name"))) {
    cli::cli_abort(
      c(
        "i" = "Graph nodes must be named when `solver = 'svd'`.",
        "x" = "Graph attribute {.var name} is missing."
      )
    )
  }
  set.seed(seed)

  if (verbose) {
    cli::cli_alert_info("Graph size: {length(g)} nodes, {igraph::gsize(g)} edges.")
  }

  # Sparse adjacency matrix
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE)

  coords <- switch(solver,
    "eigen" = .spectral_decomposition_eigen(A, normalize_laplacian, dim, verbose, ...),
    "svd" = .spectral_decomposition_svd(A, node_types, dim, verbose, ...)
  )

  coords <- normalize_layout_coordinates(coords, as_df = FALSE)

  if (jitter_sd > 0) {
    total_elements <- nrow(coords) * ncol(coords)
    coords <- coords + rnorm(total_elements, mean = 0, sd = jitter_sd)
  }

  return(coords)
}


#' Internal function to run spectral decomposition with RSpectra
#'
#' @param A An adjacency matrix (dgCMatrix)
#' @param normalize_laplacian Either TRUE or FALSE
#' @param dim The desired number of eigenvectors
#' @param verbose Print messages
#' @param ... Additional arguments passed to \code{RSpectra::eigs_sym}
#'
#' @noRd
.spectral_decomposition_eigen <- function(
  A,
  normalize_laplacian,
  dim,
  verbose,
  ...
) {
  expect_RSpectra()

  # Compute node degrees
  d <- Matrix::rowSums(A)

  if (!normalize_laplacian) {
    # Sparse diagonal degree matrix
    D <- Matrix::Diagonal(x = d)

    # Unnormalized graph Laplacian: L = D - A
    L <- D - A
  } else {
    # Create D^(-1/2)
    d_inv_sqrt <- ifelse(d > 0, 1 / sqrt(d), 0)

    # Make it a sparse diagonal matrix
    D_inv_sqrt <- Matrix::Diagonal(x = d_inv_sqrt)

    # Create the Identity matrix (I)
    I <- Matrix::Diagonal(n = nrow(A))

    # Calculate L = I - D^(-1/2) * A * D^(-1/2)
    L <- I - (D_inv_sqrt %*% A %*% D_inv_sqrt)
  }

  # Compute (dim + 1) smallest-magnitude eigenvalues and eigenvectors.
  # The smallest eigenvalue (≈ 0) yields the trivial constant eigenvector,
  # so we need one extra to obtain `dim` useful coordinates.
  k <- dim + 1
  eigen_res <- RSpectra::eigs_sym(L, k = k, which = "SA", ...)

  # Sort strictly from smallest to largest eigenvalue
  order_idx <- order(eigen_res$values)

  if (verbose) {
    cli::cli_alert_info(
      "Eigenvalues (sorted): {paste(round(eigen_res$values[order_idx], 6), collapse = ', ')}"
    )
  }

  # Drop the first column (trivial eigenvector for eigenvalue ≈ 0)
  coords <- eigen_res$vectors[, order_idx[-1], drop = FALSE]

  if (normalize_laplacian) {
    # Apply the random-walk degree correction to remove banding effects
    coords <- as.matrix(D_inv_sqrt %*% coords)
  }

  return(coords)
}

#' Internal function to run spectral decomposition with SVD
#'
#' @param A An adjacency matrix (dgCMatrix)
#' @param node_types A logical vector specifying node types for bipartite sets
#' @param dim The desired number of singular vectors
#' @param verbose Print messages
#' @param ... Additional arguments passed to \code{irlba::irlba}
#'
#' @noRd
.spectral_decomposition_svd <- function(
  A,
  node_types,
  dim,
  verbose,
  ...
) {
  expect_irlba()

  # Extract the biadjacency matrix B
  B <- A[node_types, !node_types]

  # Compute node degrees for both sets
  d1 <- Matrix::rowSums(B)
  d2 <- Matrix::colSums(B)

  # Compute the inverse square root of degrees
  d1_inv_sqrt <- ifelse(d1 > 0, 1 / sqrt(d1), 0)
  d2_inv_sqrt <- ifelse(d2 > 0, 1 / sqrt(d2), 0)

  # Create sparse diagonal scaling matrices
  D1_inv_sqrt_mat <- Matrix::Diagonal(x = d1_inv_sqrt)
  D2_inv_sqrt_mat <- Matrix::Diagonal(x = d2_inv_sqrt)

  # Construct the normalized biadjacency matrix
  B_norm <- D1_inv_sqrt_mat %*% B %*% D2_inv_sqrt_mat

  # Run Implicitly Restarted Lanczos Bidiagonalization (irlba).
  # We look for dim + 1 largest singular values. The 1st represents the
  # trivial graph component, so we need the remaining vectors for the layout.
  k <- dim + 1
  svd_res <- irlba::irlba(B_norm, nv = k, ...)

  # Extract the relevant singular vectors (drop the trivial first vector)
  U_coords <- svd_res$u[, 2:k, drop = FALSE]
  V_coords <- svd_res$v[, 2:k, drop = FALSE]

  # Apply the random-walk degree correction
  layout_V1 <- as.matrix(D1_inv_sqrt_mat %*% U_coords)
  rownames(layout_V1) <- rownames(B)

  layout_V2 <- as.matrix(D2_inv_sqrt_mat %*% V_coords)
  rownames(layout_V2) <- colnames(B)

  # Combine coordinate sets and restore the original node order of A
  coords <- rbind(layout_V1, layout_V2)
  coords <- coords[rownames(A), , drop = FALSE]

  # Ensure columns follow decreasing singular-value order
  order_idx <- order(svd_res$d[-1], decreasing = TRUE)

  if (verbose) {
    cli::cli_alert_info(
      "Singular values (sorted): {paste(round(svd_res$d[-1][order_idx], 6), collapse = ', ')}"
    )
  }

  coords <- coords[, order_idx, drop = FALSE]

  return(coords)
}
