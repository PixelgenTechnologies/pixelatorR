#' Compute local proximity scores
#'
#' This function computes a local proximity score for each node based on the
#' specified markers, using a simple neighborhood count based statistic.
#' The score is computed as the log2 ratio of the observed local proximity
#' statistic to the expected proximity statistic, where the expected statistic
#' is computed based on the analytical formula or from a null model of random
#' sampling. See details about the method in the section below.
#'
#' @section Local Neighborhood Enrichment (LNE):
#'
#' The local proximity statistic is defined by the UMI counts for the set of
#' `markers` within a neighborhood of size `k` around each node. The statistic
#' comes in two flavors depending on the `mode`:
#'
#' - all: As in "all `markers` must be present".
#' - any: As in "any marker in `markers` can be present".
#'
#' The neighborhood around a node reaches out to its `k` nearest neighbors,
#' meaning that the maximum possible distance between any two nodes in the
#' neighborhood is `2*k` edges.
#'
#' The observed statistic is calculated as either:
#'
#' all: \deqn{O_{i,m} = min(x_{i,m})}
#' any: \deqn{O_{i,m} = sum(x_{i,m})}
#'
#' where \eqn{x_{i,m}} is the vector of UMI counts for the selected `markers`
#' \eqn{m} in the neighborhood of node \eqn{i}.
#'
#' The first mode (all) is the default and computes the minimum UMI count
#' for `markers` in the neighborhood. This is used when all markers need
#' be present in the neighborhood and is useful for identifying colocation
#' of an entire set (module) of markers.
#'
#' The second mode (any) computes the sum of UMI counts for `markers`
#' in the neighborhood. This is useful when the presence of any marker
#' is sufficient to define clustering. For instance, if we want to detect
#' a region where at least one of the `markers` is present.
#'
#' Analytical formula for the expected value (\eqn{E_{i,m}}) of the local
#' proximity statistic for a single node \eqn{i} and markers \eqn{m} is given by:
#'
#' all: \deqn{E_{i,m} = min(f_{umi1} \in M_{umi1}) \times N_{umi1}(i) + min(f_{umi2} \in M_{umi2}) \times N_{umi2}(i)}
#' any: \deqn{E_{i,m} = sum(f_{umi1} \in M_{umi1}) \times N_{umi1}(i) + sum(f_{umi2} \in M_{umi2}) \times N_{umi2}(i)}
#'
#' where \eqn{M_{umi1}}, \eqn{M_{umi2}} are the sets of marker frequencies for the selected
#' `markers` for node types umi1 and umi2. The two node types (umi1/umi2) in a PNA graph
#' can have different marker frequencies and are therefore treated separately in the calculation.
#' \eqn{N_{umi1}(i)}, \eqn{N_{umi2}(i)} are the neighborhood sizes for the two node types
#' for node \eqn{i}.
#'
#' The local proximity score for node \eqn{i} is then computed as:
#' \deqn{LNE_{i,m} = log2(\frac{O_{i,m}+k-1}{E_{i,m}+k-1})}
#'
#' where \eqn{k-1} is a pseudocount used to flatten weak signals.
#'
#' A higher `k` will result in a more global proximity score, while a
#' lower `k` will result in a more local proximity score. Hence, `k` determines
#' the "spatial resolution".
#'
#' `markers` can include any number of protein IDs:
#'
#' - A single marker: local clustering
#' - Two markers: local proximity
#' - Three or more markers: local proximity for multiple markers
#'
#' Note that the score measures colocation of the specified markers. Higher values
#' indicate that the markers appear together more frequently than expected, but it
#' does not imply that they are spatially clustered together. For example, protein
#' A might be polarized while protein B is uniformly distributed, resulting in
#' high local proximity scores for A/B in the polarized region.
#'
#' @param object A `CellGraph` object containing the cell graph and counts.
#' @param markers A character vector specifying the markers to use for the
#' local proximity score.
#' @param method A character string specifying the method to use for
#' computing the local proximity score. Options are "analytical" or "permutation".
#' "analytical" is generally preferred as it is much faster. The "permutation"
#' should converge to the analytical solution as the number of iterations increases.
#' @param mode A character string specifying the mode of computation. "all" will
#' compute a joint proximity score for all selected markers. "self-clustering"
#' will compute a proximity score for each marker independently. The latter option
#' is only available for the "analytical" method which is significantly faster than
#' "permutation". With option "any", the local statistics are aggregated for
#' `markers`. In other words, the score will measure enrichment of any of the
#' `markers`.
#' @param iterations An integer specifying the number of iterations to run.
#' @param k An integer specifying the neighborhood size to consider. The minimum
#' allowed value is 2. A neighborhood size of 3 or 4 is generally recommended.
#' Note that increasing `k` lowers the spatial resolution of the score, but can
#' help detecting weaker proximity patterns in larger regions of the graph.
#' @param A_k An optional pre-computed expanded adjacency matrix (reachability matrix)
#' for neighborhood size `k`.
#' @param seed An integer for random seed setting.
#' @param ... Additional parameters
#'
#' @return A vector or matrix of local proximity scores for each node in the graph.
#'
#' @examples
#' library(ggplot2)
#' library(dplyr)
#' library(tibble)
#'
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
#'   LoadCellGraphs(cells = colnames(.)[2], add_layouts = TRUE)
#' cg <- CellGraphs(se)[[2]]
#'
#' # Compute local scores using analytical method
#' log2_ratio_analytical <- local_proximity(
#'   object = cg,
#'   markers = "B2M"
#' )
#'
#' # Compute local scores using permutation method
#' log2_ratio_permuted_20 <- local_proximity(
#'   object = cg,
#'   markers = "B2M",
#'   method = "permutation",
#'   iterations = 20
#' )
#' log2_ratio_permuted_200 <- local_proximity(
#'   object = cg,
#'   markers = "B2M",
#'   method = "permutation",
#'   iterations = 2e2
#' )
#'
#' gg <- tibble(
#'   log2_ratio_analytical = log2_ratio_analytical,
#'   log2_ratio_permuted_20 = log2_ratio_permuted_20,
#'   log2_ratio_permuted_200 = log2_ratio_permuted_200
#' )
#'
#' p1 <- ggplot(gg, aes(log2_ratio_analytical, log2_ratio_permuted_20)) +
#'   geom_point() +
#'   ggtitle("Analytical vs Permutation (20 iterations)") +
#'   theme_bw() +
#'   geom_abline(linetype = "dashed", color = "red")
#' p2 <- ggplot(gg, aes(log2_ratio_analytical, log2_ratio_permuted_200)) +
#'   geom_point() +
#'   ggtitle("Analytical vs Permutation (200 iterations)") +
#'   theme_bw() +
#'   geom_abline(linetype = "dashed", color = "red")
#'
#' # Permutation method converges to analytical with many iterations;
#' # however, the permutation method relies on fewer assumptions and
#' # is generally safer
#' p1
#' p2
#'
#' @export
#'
local_proximity <- function(
  object,
  markers,
  method = c("analytical", "permutation"),
  mode = c("all", "any", "self-clustering"),
  iterations = 50L,
  k = 3L,
  A_k = NULL,
  seed = 123,
  ...
) {
  assert_class(object, "CellGraph")
  assert_single_value(iterations, "integer")
  assert_vector(markers, "character", n = 1L)
  method <- match.arg(method, choices = c("analytical", "permutation"))
  mode <- match.arg(mode, choices = c("all", "any", "self-clustering"))
  assert_within_limits(iterations, c(1L, 10000L))
  assert_single_value(k, "integer")
  assert_within_limits(k, c(2L, 10L))
  assert_single_value(seed, "integer")
  assert_class(A_k, classes = "dgCMatrix", allow_null = TRUE)
  if (method == "permutation") {
    set.seed(seed)
  }

  if (mode == "self-clustering" && method == "permutation") {
    cli::cli_abort(
      c("x" = "self-clustering mode is not supported for permutation method.")
    )
  }

  if (!"node_type" %in% (object@cellgraph %N>% as_tibble() %>% names())) {
    cli::cli_abort(
      c("x" = "Node column {.val node_type} is missing from the graph")
    )
  }

  # Sort nodes by node type for more efficient downstream processing
  g <- object@cellgraph %N>% arrange(node_type)

  # Validate node types
  node_types <- table(g %>% pull(node_type))
  if (!all(sort(names(node_types)) == c("umi1", "umi2"))) {
    cli::cli_abort(
      c("x" = "Node types must be named {.val umi1} and {.val umi2}.")
    )
  }

  # Count nodes of each type
  nA <- node_types["umi1"] %>% as.integer()
  nB <- node_types["umi2"] %>% as.integer()

  # Sort count matrix according to the new node sorting
  counts <- object@counts[g %>% pull(name), , drop = FALSE]

  # Validate markers
  assert_x_in_y(markers, colnames(counts))

  # Get adjacency matrix and expand it to desired neighborhood size
  if (is.null(A_k)) {
    A <- igraph::as_adjacency_matrix(g)
    A_k <- expand_adjacency_matrix(A, k = k)
  } else {
    A_k <- A_k[g %>% pull(name), g %>% pull(name)]
  }

  # Include center node in the neighborhood
  Matrix::diag(A_k) <- 1L

  # Compute observed statistic
  x <- counts[, markers, drop = FALSE]
  if (mode == "self-clustering") {
    # Count the number of nodes in the neighborhood for each marker
    local_stat_obs <- as.matrix(A_k %*% x)
  }
  if (mode == "all") {
    # Take the minimum number of marker node counts for each neighborhood
    local_stat_obs <- do.call(pmin, as.data.frame(as.matrix(A_k %*% x)))
  }
  if (mode == "any") {
    # Take the sum of marker node counts for each neighborhood
    local_stat_obs <- Matrix::rowSums(A_k %*% x)
  }

  # Compute expected statistic
  local_stat_exp <- switch(method,
    "analytical" = .analytical_lp_exp(A_k, x, nA, nB, mode, ...),
    "permutation" = .permuted_lp_exp(A_k, x, nA, nB, iterations, mode)
  )

  # Compute log2 ratio of observed to expected statistic
  log2_ratio <- log2((local_stat_obs + k - 1) / (local_stat_exp + k - 1))

  # Sort results to match input graph
  log2_ratio <- switch(mode,
    "self-clustering" = log2_ratio[object@cellgraph %N>% pull(name), ],
    "all" = log2_ratio[object@cellgraph %N>% pull(name)],
    "any" = log2_ratio[object@cellgraph %N>% pull(name)]
  )

  return(log2_ratio)
}

#' @param A_k Expanded adjacency matrix
#' @param counts Count matrix of nodes
#' @param nA Number of nodes of type umi1
#' @param nB Number of nodes of type umi2
#' @param mode Mode of local proximity calculation
#' @param ... Additional parameters
#'
#' @noRd
.analytical_lp_exp <- function(A_k, counts, nA, nB, mode, ...) {
  dots <- list(...)
  if (all(c("freqs1", "freqs2") %in% names(dots))) {
    freq1 <- dots$freqs1
    freq2 <- dots$freqs2
  }
  # For the analytical method, we compute the expected node counts
  # using the marker frequencies multiplied by the neighborhood size.
  # UMI1/UMI2 are handled separately since they have different frequencies.
  degree_exp_A <- A_k[, 1:nA] %>% Matrix::rowSums()
  degree_exp_B <- A_k[, (nA + 1):(nA + nB)] %>% Matrix::rowSums()

  # Compute marker frequencies for each node type
  if (!exists("freq1", inherits = FALSE)) {
    freq1 <- (counts[1:nA, , drop = FALSE] %>% Matrix::colSums()) / nA
  }
  if (mode != "self-clustering") {
    if (mode == "all") {
      freq1 <- freq1 %>% min()
    }
    if (mode == "any") {
      freq1 <- freq1 %>% sum()
    }
    local_stat_exp_A <- freq1 * degree_exp_A
  } else {
    local_stat_exp_A <- outer(degree_exp_A, freq1)
  }
  if (!exists("freq2", inherits = FALSE)) {
    freq2 <- (counts[(nA + 1):(nA + nB), , drop = FALSE] %>% Matrix::colSums()) / nB
  }
  if (mode != "self-clustering") {
    if (mode == "all") {
      freq2 <- freq2 %>% min()
    }
    if (mode == "any") {
      freq2 <- freq2 %>% sum()
    }
    local_stat_exp_B <- freq2 * degree_exp_B
  }
  if (mode == "self-clustering") {
    local_stat_exp_B <- outer(degree_exp_B, freq2)
  }

  # Combine the expected statistics for both node types
  local_stat_exp <- local_stat_exp_A + local_stat_exp_B

  return(local_stat_exp)
}

#' @param A_k Expanded adjacency matrix
#' @param x Count matrix of nodes for selected markers
#' @param nA Number of nodes of type umi1
#' @param nB Number of nodes of type umi2
#' @param iterations Number of iterations for permutation step
#' @param mode Mode of local proximity calculation
#' @noRd
.permuted_lp_exp <- function(A_k, x, nA, nB, iterations, mode) {
  # For the permutation method, we compute the expected node counts
  # by permuting the node marker labels, computing the neighborhood marker
  # node counts and averaging across all iterations.  UMI1/UMI2 are handled
  # separately since they have different frequencies.
  # Accumulate a running sum to avoid storing an n_nodes x iterations matrix.
  local_stat_sum <- numeric(nrow(x))
  names(local_stat_sum) <- rownames(x)

  # Permute the counts matrix and compute statistic
  for (i in 1L:iterations) {
    # Shuffle nodes independently
    x_shuffled <- x[c(sample.int(nA), sample.int(nB) + nA), ]
    if (mode == "all") {
      jcs <- as.data.frame(as.matrix(A_k %*% x_shuffled))
      local_stat_sum <- local_stat_sum + do.call(pmin, jcs)
    }
    if (mode == "any") {
      local_stat_sum <- local_stat_sum + Matrix::rowSums(A_k %*% x_shuffled)
    }
  }

  # Compute expected statistic
  local_stat_exp <- local_stat_sum / iterations
  return(local_stat_exp)
}
