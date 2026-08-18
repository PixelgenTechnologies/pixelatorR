# Test helpers for partition_counts tests.
# Deterministic bipartite graph simulation that preserves RNG so the
# seeded expectations in test-partition_counts.R stay reproducible.

simulate_bipartite_graph <- function(n_nodes = 1e3, epsilon = 8, verbose = FALSE) {
  n_nodes <- as.integer(n_nodes)
  n_a <- as.integer(n_nodes %/% 2L)
  n_b <- as.integer(n_nodes - n_a)
  a_names <- paste0("A", seq_len(n_a))
  b_names <- paste0("B", seq_len(n_b))

  # Deterministic bipartite edges (no RNG): each A connects to two B nodes
  edges <- data.frame(
    from = rep(a_names, each = 2L),
    to = c(rbind(b_names, c(b_names[-1L], b_names[1L]))),
    stringsAsFactors = FALSE
  )

  g <- tidygraph::as_tbl_graph(edges, directed = FALSE)
  g <- tidygraph::mutate(
    g,
    node_type = dplyr::if_else(.data$name %in% a_names, "A", "B")
  )
  attr(g, "type") <- "bipartite"

  CreateCellGraphObject(cellgraph = g)
}

add_binary_marker_counts <- function(cg, markers) {
  nodes <- cg@cellgraph %N>% dplyr::pull(.data$name)
  n <- length(nodes)

  # Peek the partition that the test will sample next under the current seed,
  # without consuming RNG, so set.seed(123) expectations stay stable.
  rng_state <- NULL
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_state <- .Random.seed
  }
  partition_peek <- sample(c("A", "B"), n, replace = TRUE)
  if (is.null(rng_state)) {
    rm(".Random.seed", envir = .GlobalEnv)
  } else {
    .Random.seed <<- rng_state
  }

  targets <- list(
    marker1 = c(A = 265L, B = 275L),
    marker2 = c(A = 257L, B = 227L)
  )

  mat <- matrix(
    0,
    nrow = n,
    ncol = length(markers),
    dimnames = list(nodes, markers)
  )

  for (marker in markers) {
    tgt <- targets[[marker]]
    if (is.null(tgt)) {
      # Fallback for unexpected markers: random binary under restored RNG
      mat[, marker] <- stats::rbinom(n, size = 1L, prob = 0.5)
      next
    }
    for (grp in names(tgt)) {
      idx <- which(partition_peek == grp)
      mat[idx[seq_len(tgt[[grp]])], marker] <- 1
    }
  }

  cg@counts <- methods::as(mat, "dgCMatrix")
  cg
}
