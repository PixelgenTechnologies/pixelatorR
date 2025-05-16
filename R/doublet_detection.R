#' @include generics.R
NULL

#' Approximate nearest neighbors using Annoy
#'
#' @description
#' This function takes a matrix, and returns approximate Euclidean nearest
#' neighbors and distances of row items given the number of trees (n_trees)
#' and number of nearest neighbors (n_nn).
#'
#'
#' @param dat A matrix with data to find nearest neighbors. Rows are cells, and columns are features.
#' @param cells A character vector with cell names to find nearest neighbors for. If NULL, all cells are used.
#' @param n_trees An integer with the number of trees to build in the Annoy index.
#' @param n_nn An integer with the number of nearest neighbors to find.
#' @param search_k An integer with the number of nodes to search in the Annoy index. Default is `n_trees * n_nn`.
#' @param annoy_alg An character specifying which distance algorithm to use. Default is
#'                  \code{\link[RcppAnnoy]{AnnoyEuclidean}} (`"euclidean"`).
#'                  Available options are `"euclidean"`, `"angular"`, `"manhattan"`, and `"hamming"`.
#'
#'
#' @examples
#' dat <- matrix(rnorm(1000), ncol = 10)
#' FindAnnoyNeighbors(dat, n_trees = 50, n_nn = 10)
#'
#' @return A tibble with the following columns:
#' \itemize{
#'    \item{id: The cell name}
#'    \item{index: The index of the cell in the Annoy index}
#'    \item{item: The index of the nearest neighbor}
#'    \item{distance: The distance (Euclidean by default) to the nearest neighbor}
#'    \item{nn: The rank of the nearest neighbor}
#'    \item{neighbor: The name of the nearest neighbor}
#' }
#'
#' @export
#'
FindAnnoyNeighbors <- function(
  dat,
  cells = NULL,
  n_trees = 50,
  n_nn = 10,
  search_k = NULL,
  annoy_alg = c("euclidean", "angular", "manhattan", "hamming")
) {
  stopifnot(
    "'dat' must be a matrix" =
      inherits(dat, what = "matrix"),
    "'n_trees' must be an integer" =
      inherits(n_trees, what = "numeric") && (n_trees > 0),
    "'n_nn' must be an integer" =
      inherits(n_nn, what = "numeric") && (n_nn > 0),
    "'search_k' must be an integer" =
      is.null(search_k) |
        inherits(search_k, what = "numeric"),
    "'dat' must have row names if 'cells' are specified" =
      is.null(cells) |
        (inherits(cells, what = "character") &&
          all(cells %in% rownames(dat)))
  )

  expect_RcppAnnoy()

  annoy_alg <-
    switch(match.arg(
      annoy_alg,
      c(
        "euclidean",
        "angular",
        "manhattan",
        "hamming"
      )
    ),
    euclidean = RcppAnnoy::AnnoyEuclidean,
    angular = RcppAnnoy::AnnoyAngular,
    manhattan = RcppAnnoy::AnnoyManhattan,
    hamming = RcppAnnoy::AnnoyHamming
    )

  nrows <- nrow(dat)
  ncols <- ncol(dat)
  annoy_index <- methods::new(annoy_alg, ncols)

  if (is.null(rownames(dat))) {
    rownames(dat) <- seq_len(nrows)
  }

  if (is.null(search_k)) {
    search_k <- -1 # -1 invokes default for RcppAnnoy
  }

  for (i in 1:nrows) {
    annoy_index$addItem(i - 1, dat[i, ])
  }

  annoy_index$build(n_trees)

  if (is.null(cells)) {
    cells <- rownames(dat)
  }

  tibble(id = rownames(dat)) %>%
    mutate(index = row_number() - 1) %>%
    filter(id %in% cells) %>%
    group_by_all() %>%
    reframe({

      annoy_index$getNNsByItemList(index,
                                   n = n_nn,
                                   search_k = search_k,
                                   include_distances = TRUE
      ) %>%
        as_tibble() %>%
        mutate(nn = row_number())

    }) %>%
    ungroup() %>%
    mutate(neighbor = rownames(dat)[item + 1])

}

#' Simulate doublets
#'
#' Simulates doublets by summing pairs of randomly selected cells from a matrix.
#'
#' @param count_data A matrix with data to simulate doublets. Rows are features, and columns are cells.
#' @param ref_cells1,ref_cells2 A character vector with cell names or indices to use as the first and second reference
#'                              populations. If NULL, all cells are used.
#' @param n_sim An integer with the number of doublets to simulate.
#' @param seed An integer with the seed for the random number generator.
#' @param method A character with the method to use to simulate doublets. Default is "average".
#'
#'
#' @examples
#' dat <- matrix(rnorm(1000), ncol = 10)
#' SimulateDoublets(dat, n_sim = 1000)
#'
#' @return A matrix with simulated doublets. Rows are features, and columns are simulated doublets.
#'
#' @export
#'
SimulateDoublets <- function(
  count_data,
  ref_cells1 = NULL,
  ref_cells2 = NULL,
  n_sim = 1000,
  seed = 37,
  method = c("average", "sum")
) {
  assert_class(count_data, c("matrix", "Matrix"))
  assert_class(ref_cells1, "character", allow_null = TRUE)
  assert_class(ref_cells2, "character", allow_null = TRUE)
  assert_class(n_sim, "numeric")
  assert_within_limits(n_sim, c(1, Inf))
  assert_class(seed, "numeric")
  method <- match.arg(method, c("average", "sum"))

  # Create doublet template data
  if (is.null(ref_cells1) & is.null(ref_cells2)) {
    ref_pop1 <- count_data
    ref_pop2 <- count_data
  } else if (!is.null(ref_cells1) & !is.null(ref_cells2)) {
    ref_pop1 <- count_data[, ref_cells1]
    ref_pop2 <- count_data[, ref_cells2]
  } else {
    cli_abort("Please provide either `ref_cells1` and `ref_cells2` or both.")
  }

  n_cells1 <-
    ncol(ref_pop1)
  n_cells2 <-
    ncol(ref_pop2)



  cli_alert_info(glue("Simulating {n_sim} doublets"))
  set.seed(seed)
  random_doublets <-
    tibble(
      i1 = sample(seq_len(n_cells1), n_sim, replace = TRUE),
      i2 = sample(seq_len(n_cells2), n_sim, replace = TRUE)
    )


  simulated_doublets <-
    (ref_pop1[, random_doublets$i1] +
       ref_pop2[, random_doublets$i2])

  if (method == "average") {
    simulated_doublets <-
      simulated_doublets / 2
  }

  colnames(simulated_doublets) <-
    paste0("sim_", seq_len(n_sim))

  return(simulated_doublets)
}


#' @rdname PredictDoublets
#' @method PredictDoublets Matrix
#'
#' @export
#'
PredictDoublets.Matrix <- function(
  object,
  ref_cells1 = NULL,
  ref_cells2 = NULL,
  simulation_rate = 3,
  n_neighbor = 100,
  npcs = 10,
  p_adjust_method = "BH",
  p_threshold = 0.01,
  seed = 37,
  ...
) {
  assert_class(object, "Matrix")
  assert_class(ref_cells1, "character", allow_null = TRUE)
  assert_class(ref_cells2, "character", allow_null = TRUE)
  assert_class(simulation_rate, "numeric")
  assert_within_limits(simulation_rate, c(0, Inf))
  assert_class(n_neighbor, "numeric")
  assert_within_limits(n_neighbor, c(1, Inf))
  assert_class(npcs, "numeric")
  assert_within_limits(npcs, c(1, Inf))
  assert_class(p_adjust_method, "character")
  assert_x_in_y(p_adjust_method, stats::p.adjust.methods)
  assert_class(p_threshold, "numeric")
  assert_within_limits(p_threshold, c(0, 1))
  assert_class(seed, "numeric")

  expect_pcaMethods()

  n_cells <-
    ncol(object)

  cell_names <-
    tibble(
      original_names = colnames(object),
      temp_names = paste0("cell_", seq_len(n_cells))
    )

  colnames(object) <-
    cell_names$temp_names

  simulated_doublets <-
    SimulateDoublets(
      count_data = object,
      ref_cells1 = ref_cells1,
      ref_cells2 = ref_cells2,

      n_sim = round(simulation_rate * n_cells),
      seed = seed,
      method = "sum"
    )

  pca_prep <-
    . %>%
    as.matrix() %>%
    log1p() %>%
    sweep(2, apply(., 2, mean), `-`) %>%
    t() %>%
    scale()


  cli_alert_info(glue("Performing PCA on combined data"))
  pca_scores <-
    cbind(
      object,
      simulated_doublets
    ) %>%
    pca_prep() %>%
    pcaMethods::pca(
      method = "ppca",
      nPcs = npcs
    ) %>%
    slot("scores")

  cli_alert_info(glue("Predicting doublets using {n_neighbor} nearest neighbors"))
  cell_nn <-
    pca_scores %>%
    FindAnnoyNeighbors(
      cells = cell_names$temp_names,
      n_nn = n_neighbor + 1
    )

  expected_doublet_rate <- simulation_rate / (1 + simulation_rate)
  expected_doublet_nns <- n_neighbor * expected_doublet_rate

  doublet_prediction <-
    cell_nn %>%
    filter(id != neighbor) %>%
    mutate(simulated = !neighbor %in% cell_names$temp_names) %>%
    group_by(id) %>%
    summarise(
      doublet_nns = sum(simulated),
      doublet_nn_rate = doublet_nns / n_neighbor,
      .groups = "drop"
    ) %>%
    mutate(
      doublet_p = stats::pbinom(doublet_nns - 1, n_neighbor, expected_doublet_rate, lower.tail = FALSE),
      doublet_p_adj = p.adjust(doublet_p, p_adjust_method),
      doublet_prediction = ifelse(doublet_p_adj < p_threshold, "doublet", "singlet")
    ) %>%
    left_join(cell_names,
      by = c("id" = "temp_names")
    ) %>%
    select(-id) %>%
    column_to_rownames("original_names")

  return(doublet_prediction)
}


#' @rdname PredictDoublets
#' @method PredictDoublets Seurat
#'
#' @export
#'
PredictDoublets.Seurat <- function(
  object,
  ref_cells1 = NULL,
  ref_cells2 = NULL,
  simulation_rate = 3,
  n_neighbor = 100,
  npcs = 10,
  p_adjust_method = "BH",
  p_threshold = 0.01,
  seed = 37,
  assay = NULL,
  layer = "counts",
  ...
) {
  assert_class(object, "Seurat")

  LayerData(object, assay = assay, layer = layer) %>%
    PredictDoublets(
      simulation_rate = simulation_rate,
      n_neighbor = n_neighbor,
      npcs = npcs,
      p_adjust_method = p_adjust_method,
      p_threshold = p_threshold,
      seed = seed
    )
}


