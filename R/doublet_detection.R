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
#' @param x A numeric matrix with data to find nearest neighbors. Rows are cells, and columns are features.
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
#' x <- matrix(rnorm(1000), ncol = 10)
#' FindAnnoyNeighbors(x, n_trees = 50, n_nn = 10)
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
  x,
  cells = NULL,
  n_trees = 50L,
  n_nn = 10L,
  search_k = NULL,
  annoy_alg = c("euclidean", "angular", "manhattan", "hamming")
) {
  assert_class(x, c("matrix", "Matrix"))
  assert_class(cells, c("character", "integer"), allow_null = TRUE)
  assert_x_in_y(cells, rownames(x), allow_null = TRUE)
  assert_class(n_trees, c("numeric", "integer"))
  assert_within_limits(n_trees, c(1, Inf))
  assert_class(n_nn, c("numeric", "integer"))
  assert_within_limits(n_nn, c(1, Inf))
  assert_class(search_k, c("numeric", "integer"), allow_null = TRUE)
  assert_class(annoy_alg, "character")
  assert_x_in_y(annoy_alg, c("euclidean", "angular", "manhattan", "hamming"))

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

  nrows <- nrow(x)
  ncols <- ncol(x)
  annoy_index <- methods::new(annoy_alg, ncols)

  if (is.null(rownames(x))) {
    rownames(x) <- seq_len(nrows)
  }

  if (is.null(search_k)) {
    search_k <- -1 # -1 invokes default for RcppAnnoy
  }

  for (i in 1:nrows) {
    annoy_index$addItem(i - 1, x[i, ])
  }

  annoy_index$build(n_trees)

  if (is.null(cells)) {
    cells <- rownames(x)
  }

  tibble(id = rownames(x)) %>%
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
    mutate(neighbor = rownames(x)[item + 1])
}

#' Simulate doublets
#'
#' Simulates doublets by summing pairs of randomly selected cells from a matrix.
#'
#' @param count_data A count matrix. Rows are features, and columns are cells.
#' @param ref_cells1,ref_cells2 A character vector with cell names or indices to use as the first and second reference
#'                              populations. If NULL, all cells are used.
#' @param n_sim An integer with the number of doublets to simulate.
#' @param seed An integer with the seed for the random number generator.
#' @param method A character with the method to use to simulate doublets. Options are "average" (average of the two
#'               cells) or "sum" (sum of the two cells).
#' @param simulated_cell_prefix A character with the prefix to use for the simulated doublet cell names.
#' @param return_id If TRUE, returns a list with the simulated doublet counts and a tibble withe their original cell
#'                  IDs. If FALSE, returns only the simulated doublet counts.
#' @param verbose Print messages.
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
  method = c("average", "sum"),
  simulated_cell_prefix = "simcell_",
  return_id = FALSE,
  verbose = TRUE
) {
  assert_class(count_data, c("matrix", "Matrix"))
  assert_class(ref_cells1, c("character", "integer"), allow_null = TRUE)
  assert_class(ref_cells2, c("character", "integer"), allow_null = TRUE)
  assert_class(n_sim, "numeric")
  assert_within_limits(n_sim, c(1, Inf))
  assert_class(seed, "numeric")
  method <- match.arg(method, c("average", "sum"))

  # Create doublet template data
  if (is.null(ref_cells1) && is.null(ref_cells2)) {
    ref_pop1 <- count_data
    ref_pop2 <- count_data
  } else if (!is.null(ref_cells1) && !is.null(ref_cells2)) {
    ref_pop1 <- count_data[, ref_cells1]
    ref_pop2 <- count_data[, ref_cells2]
  } else {
    cli::cli_abort("Please provide either `ref_cells1` and `ref_cells2` or both.")
  }

  n_cells1 <-
    ncol(ref_pop1)
  n_cells2 <-
    ncol(ref_pop2)


  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info(glue("Simulating {n_sim} doublets"))
  }
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
    paste0(simulated_cell_prefix, seq_len(n_sim))

  if (return_id) {
    return(
      list(
        counts = simulated_doublets,
        ids = random_doublets %>%
          mutate(
            id1 = colnames(count_data)[i1],
            id2 = colnames(count_data)[i2]
          )
      )
    )
  } else {
    return(simulated_doublets)
  }
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
  simulation_rate = 1,
  n_neighbor = 100,
  npcs = 10,
  p_adjust_method = "BH",
  p_threshold = 0.05,
  seed = 37,
  iter = 1,
  return_trials = FALSE,
  verbose = TRUE,
  ...
) {
  assert_class(object, "Matrix")
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

  if (verbose && check_global_verbosity()) {
    cli::cli_alert_info(glue("Predicting doublets using {n_neighbor} nearest neighbors"))
  }

  if (verbose && check_global_verbosity() && iter > 1) {
    pb <- cli_progress_bar("Predicting doublets", total = iter, .auto_close = FALSE)
  }

  pca_prep <-
    . %>%
    as.matrix() %>%
    log1p() %>%
    sweep(2, apply(., 2, mean), `-`) %>%
    t() %>%
    scale()

  cell_nn_res <-
    lapply(seq_len(iter),
           function(i) {

             if (verbose && check_global_verbosity() && iter > 1) {
               cli_progress_update(id = pb)
             }

             simulated_doublets <-
               SimulateDoublets(
                 count_data = object,
                 ref_cells1 = ref_cells1,
                 ref_cells2 = ref_cells2,
                 n_sim = round(simulation_rate * n_cells),
                 seed = seed + i - 1,
                 method = "sum",
                 verbose = FALSE
               )

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

             cell_nn <-
               pca_scores %>%
               FindAnnoyNeighbors(
                 cells = colnames(object),
                 n_nn = n_neighbor + 1
               ) %>%
               filter(id != neighbor) %>%
               mutate(simulated = !neighbor %in% colnames(object)) %>%
               group_by(id) %>%
               summarise(
                 doublet_nns = sum(simulated),
                 doublet_nn_rate = doublet_nns / n_neighbor,
                 .groups = "drop"
               )

           })
  if (verbose && check_global_verbosity() && iter > 1) {
    cli_progress_done()
  }

  expected_doublet_rate <- simulation_rate / (1 + simulation_rate)
  expected_doublet_nns <- n_neighbor * expected_doublet_rate

  calc_stats <-
    . %>%
    mutate(
      doublet_p = stats::pbinom(doublet_nns - 1, n_neighbor, expected_doublet_rate, lower.tail = FALSE),
      doublet_p_adj = p.adjust(doublet_p, p_adjust_method),
      logratio = log2(doublet_nns / expected_doublet_nns),
      doublet_prediction = ifelse(doublet_p_adj < p_threshold, "doublet", "singlet"),
    )

  doublet_prediction_trials <-
    cell_nn_res %>%
    bind_rows(.id = "trial") %>%
    mutate(trial = as.numeric(trial)) %>%
    calc_stats() %>%
    arrange(factor(id, levels = colnames(object)),
            trial)

  if(return_trials) {
    return(doublet_prediction_trials)
  }

  doublet_prediction <-
    doublet_prediction_trials %>%
    group_by(id) %>%
    summarise(doublet_nns = round(median(doublet_nns), 0),
              doublet_vote = sum(doublet_prediction == "doublet") / n()) %>%
    calc_stats() %>%
    arrange(factor(id, levels = colnames(object)))

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
  simulation_rate = 1,
  n_neighbor = 100,
  npcs = 10,
  p_adjust_method = "BH",
  p_threshold = 0.05,
  seed = 37,
  iter = 1,
  assay = NULL,
  layer = "counts",
  verbose = TRUE,
  ...
) {
  assert_class(object, "Seurat")
  assert_class(assay, "character", allow_null = TRUE)
  assert_class(layer, "character", allow_null = TRUE)

  LayerData(object, assay = assay, layer = layer) %>%
    PredictDoublets(
      ref_cells1 = ref_cells1,
      ref_cells2 = ref_cells2,
      simulation_rate = simulation_rate,
      n_neighbor = n_neighbor,
      npcs = npcs,
      p_adjust_method = p_adjust_method,
      p_threshold = p_threshold,
      seed = seed,
      iter = iter,
      verbose = verbose
    )
}
