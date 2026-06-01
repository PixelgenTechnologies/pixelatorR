#' @rdname ComputeProximityScores
#' @method ComputeProximityScores CellGraph
#'
#' @examples
#' library(dplyr)
#' cg <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
#'   LoadCellGraphs(cells = colnames(.)[1], verbose = FALSE) %>%
#'   CellGraphs() %>%
#'   .[[1]]
#'
#' ComputeProximityScores(cg) %>%
#'   filter(join_count_z > 3)
#'
#' @export
#'
ComputeProximityScores.CellGraph <- function(
  object,
  mode = c("analytical", "permutation"),
  k = NULL,
  iterations = 100L,
  calc_z_score = TRUE,
  calc_log2_ratio = TRUE,
  min_marker_count = 10,
  seed = 123,
  ...
) {
  assert_class(object, "CellGraph")
  mode <- match.arg(mode, c("analytical", "permutation"))
  if (mode == "permutation") {
    expect_matrixStats()
  }
  assert_single_value(k, "integer", allow_null = TRUE)
  if (is.null(k)) {
    k <- 1L
  } else {
    lifecycle::signal_stage("experimental", "ComputeProximityScores(k = )")
    assert_within_limits(k, c(1L, 6L))
  }
  assert_single_value(iterations, "integer")
  assert_within_limits(iterations, c(1L, 10000L))
  assert_single_value(calc_z_score, "bool")
  assert_single_value(calc_log2_ratio, "bool")
  assert_single_value(min_marker_count, "integer")
  assert_single_value(seed, "integer")

  if (!"node_type" %in% (object@cellgraph %N>% as_tibble() %>% names())) {
    cli::cli_abort(
      c("x" = "Node column {.val node_type} is missing from the graph")
    )
  }

  g <- object@cellgraph %N>% arrange(node_type)
  node_types_vector <- g %>% pull(node_type)
  node_types <- table(node_types_vector)
  if (!all(sort(names(node_types)) == c("umi1", "umi2"))) {
    cli::cli_abort(
      c("x" = "Node types must be named {.val umi1} and {.val umi2}.")
    )
  }
  nA <- node_types["umi1"] %>% as.integer()
  nB <- node_types["umi2"] %>% as.integer()

  counts <- object@counts[g %>% pull(name), ]

  # Get adjacency matrix. We only need the upper triangle
  if (k > 1L) {
    A <- igraph::as_adjacency_matrix(g)
    A <- expand_adjacency_matrix(A, k = k)
    A <- Matrix::triu(A, k = 1)
  } else {
    A <- igraph::as_adjacency_matrix(g, type = "upper")
  }

  # Ignore proteins with less than min_marker_count counts
  counts <- counts[, Matrix::colSums(counts) >= min_marker_count, drop = FALSE]
  if (ncol(counts) == 0) {
    cli::cli_abort(
      c("x" = "No proteins with at least {.val {min_marker_count}} counts found.")
    )
  }

  markers <- colnames(counts)

  # Compute observed join counts
  jc_obs <- Matrix::t(counts) %*% A %*% counts

  if (mode == "analytical") {
    N_edges <- sum(A)
    fA <- Matrix::colSums(counts[node_types_vector == "umi1", , drop = FALSE]) / nA
    fB <- Matrix::colSums(counts[node_types_vector == "umi2", , drop = FALSE]) / nB
    join_count_expected_mean <- outer(fA, fB) * N_edges
    join_count_expected_mean <- join_count_expected_mean + t(join_count_expected_mean)
    diag(join_count_expected_mean) <- diag(join_count_expected_mean) / 2
    join_count_expected_mean <- join_count_expected_mean[upper.tri(join_count_expected_mean, diag = TRUE)]

    S1 <- (sum((A + Matrix::t(A))^2)) / 2
    S2 <- sum((Matrix::colSums(A) + Matrix::rowSums(A))^2)

    join_count_expected_var <- (
      (2 * S1 * outer(fA, fB)) +
        (S2 - (2 * S1)) * (outer(fA, fB) * outer(fA, fB, "+")) +
        (4 * (S1 - S2) * outer(fA^2, fB^2))
    ) / 4

    join_count_expected_var <- (join_count_expected_var + t(join_count_expected_var))
    diag(join_count_expected_var) <- diag(join_count_expected_var) / 2
    join_count_expected_sd <- sqrt(join_count_expected_var)
    join_count_expected_sd <- join_count_expected_sd[upper.tri(join_count_expected_sd, diag = TRUE)]
  }
  if (mode == "permutation") {
    set.seed(seed)
    N <- ncol(counts)

    # Initialize the array to store the permuted join counts
    ar <- array(data = 0, dim = c(N, N, iterations))
    # Permute the counts matrix and compute join counts
    for (i in 1L:iterations) {
      # Shuffle nodes independently
      counts <- counts[c(sample.int(nA), sample.int(nB) + nA), ]
      jc_perm <- Matrix::t(counts) %*% A %*% counts
      ar[, , i] <- jc_perm %>% as.matrix()
    }

    # Add the transpose to ignore direction, e.g. A/B = B/A
    ar <- ar + aperm(ar, c(2, 1, 3))

    # Cast to 2D matrix and fast compute the mean and sd
    ar_2d <- matrix(ar, N^2, iterations)
    avgs <- ar_2d %>%
      Matrix::rowMeans() %>%
      matrix(ncol = N)
    diag(avgs) <- diag(avgs) / 2
    sds <- ar_2d %>%
      matrixStats::rowSds() %>%
      matrix(ncol = N)
    diag(sds) <- diag(sds) / 2

    join_count_expected_mean <- avgs[upper.tri(avgs, diag = TRUE)]
    join_count_expected_sd <- sds[upper.tri(sds, diag = TRUE)]
  }

  # Collect stats
  jc_obs <- jc_obs + Matrix::t(jc_obs)
  Matrix::diag(jc_obs) <- Matrix::diag(jc_obs) / 2
  arr_inds <- which(upper.tri(jc_obs, diag = TRUE), arr.ind = TRUE)
  jc_obs <- jc_obs[upper.tri(jc_obs, diag = TRUE)]

  old_locale <- Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  on.exit(Sys.setlocale("LC_COLLATE", old_locale), add = TRUE)

  # Create final tibble
  proximity_scores <- arr_inds %>%
    as_tibble() %>%
    mutate(
      marker_1_unordered = markers[row],
      marker_2_unordered = markers[col],
      join_count = jc_obs,
      join_count_expected_mean = join_count_expected_mean,
      join_count_expected_sd = join_count_expected_sd
    ) %>%
    select(-row, -col) %>%
    mutate(
      marker_1 = pmin(marker_1_unordered, marker_2_unordered),
      marker_2 = pmax(marker_1_unordered, marker_2_unordered)
    ) %>%
    select(-marker_1_unordered, -marker_2_unordered)

  if (calc_z_score) {
    proximity_scores <- proximity_scores %>%
      mutate(join_count_z = (join_count - join_count_expected_mean) / pmax(join_count_expected_sd, 1))
  }
  if (calc_log2_ratio) {
    proximity_scores <- proximity_scores %>%
      mutate(log2_ratio = log2(pmax(join_count, 1) / pmax(join_count_expected_mean, 1)))
  }

  return(proximity_scores)
}

#' @param cl Number of threads to use for parallel processing.
#' Only used on unix systems. If `NULL`, sequential processing is used.
#'
#' @rdname ComputeProximityScores
#' @method ComputeProximityScores list
#'
#' @export
#'
ComputeProximityScores.list <- function(
  object,
  mode = c("analytical", "permutation"),
  k = 1L,
  iterations = 100L,
  calc_z_score = TRUE,
  calc_log2_ratio = TRUE,
  min_marker_count = 10,
  seed = 123,
  cl = NULL,
  ...
) {
  assert_single_value(cl, type = "integer", allow_null = TRUE)
  if (!all(sapply(object, class) == "CellGraph")) {
    cli::cli_abort(
      c("x" = "All elements in {.var object} must be of class {.cls CellGraph}")
    )
  }
  if (is.null(names(object))) {
    cli::cli_abort(
      c("x" = "The list {.var object} must be named")
    )
  } else {
    if (length(object) != length(unique(names(object)))) {
      cli::cli_abort(
        c("x" = "The list {.var object} must have unique names")
      )
    }
  }

  proximity_scores <- pbapply::pblapply(names(object), function(nm) {
    ComputeProximityScores(object[[nm]], mode, k, iterations, calc_z_score, calc_log2_ratio, min_marker_count, seed) %>%
      mutate(component = nm)
  }, cl = cl) %>%
    bind_rows()

  return(proximity_scores)
}


#' @param cells A vector of cell names to compute proximity scores for.
#'
#' @rdname ComputeProximityScores
#' @method ComputeProximityScores PNAAssay
#'
#' @export
#'
ComputeProximityScores.PNAAssay <- function(
  object,
  mode = c("analytical", "permutation"),
  k = 1L,
  cells = NULL,
  iterations = 100L,
  calc_z_score = TRUE,
  calc_log2_ratio = TRUE,
  min_marker_count = 10,
  seed = 123,
  cl = NULL,
  ...
) {
  assert_vector(cells, type = "character", allow_null = TRUE, n = 1)
  assert_x_in_y(cells, colnames(object), allow_null = TRUE)

  cg_list <- CellGraphs(object)
  if (!is.null(cells)) {
    cg_list <- cg_list[cells]
  }

  if (!all(sapply(cg_list, class) == "CellGraph")) {
    cli::cli_abort(
      c(
        "i" = "Make sure to run {.fn LoadCellGraphs} first.",
        "x" = "All elements in {.var cg_list} must be of class {.cls CellGraph}"
      )
    )
  }
  if (is.null(names(cg_list))) {
    cli::cli_abort(
      c("x" = "The list {.var cg_list} must be named")
    )
  } else {
    if (length(cg_list) != length(unique(names(cg_list)))) {
      cli::cli_abort(
        c("x" = "The list {.var cg_list} must have unique names")
      )
    }
  }

  proximity_scores <- cg_list %>%
    ComputeProximityScores(
      mode,
      k,
      iterations,
      calc_z_score,
      calc_log2_ratio,
      min_marker_count,
      seed,
      cl
    )

  return(proximity_scores)
}


#' @rdname ComputeProximityScores
#' @method ComputeProximityScores PNAAssay5
#'
#' @export
#'
ComputeProximityScores.PNAAssay5 <- ComputeProximityScores.PNAAssay


#' @param assay Name of the \code{PNAAssay} or \code{PNAAssay5} to use for
#' computing proximity scores.
#'
#' @rdname ComputeProximityScores
#' @method ComputeProximityScores Seurat
#'
#' @examples
#' library(ggplot2)
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
#'   LoadCellGraphs(cells = colnames(.)[1:2], verbose = FALSE)
#'
#' # Compute proximity scores for selected cells
#' proximity_scores <- ComputeProximityScores(se, cells = colnames(se)[1:2])
#'
#' # Compare with available proximity scores
#' proximity_scores %>%
#'   left_join(ProximityScores(se %>% subset(cells = colnames(se)[1:2])),
#'     by = c("marker_1", "marker_2", "component"),
#'     suffix = c("_post", "_pre")
#'   ) %>%
#'   na.omit() %>%
#'   ggplot(aes(log2_ratio_pre, log2_ratio_post)) +
#'   geom_abline() +
#'   geom_point() +
#'   theme_bw() +
#'   labs(x = "log2_ratio pre-computed", y = "log2_ratio computed") +
#'   facet_grid(~component)
#'
#' @export
#'
ComputeProximityScores.Seurat <- function(
  object,
  mode = c("analytical", "permutation"),
  k = 1L,
  assay = NULL,
  cells = NULL,
  iterations = 100L,
  calc_z_score = TRUE,
  calc_log2_ratio = TRUE,
  min_marker_count = 10,
  seed = 123,
  cl = NULL,
  ...
) {
  if (!is.null(assay)) {
    assert_single_value(assay, type = "string")
  } else {
    # Use default assay if assay = NULL
    assay <- DefaultAssay(object)
  }

  pixel_assay <- object[[assay]]
  assert_pixel_assay(pixel_assay)

  proximity_scores <-
    ComputeProximityScores(
      pixel_assay,
      mode,
      k,
      cells,
      iterations,
      calc_z_score,
      calc_log2_ratio,
      min_marker_count,
      seed,
      cl
    )
  return(proximity_scores)
}
