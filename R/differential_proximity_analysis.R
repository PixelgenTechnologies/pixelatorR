#' @include generics.R
NULL

.datatable.aware <- TRUE


#' @rdname DifferentialProximityAnalysis
#' @method DifferentialProximityAnalysis data.frame
#'
#' @examples
#' library(dplyr)
#' example_data <- tidyr::expand_grid(
#'   marker_1 = c("HLA-ABC", "B2M", "CD4", "CD8", "CD20", "CD19", "CD45", "CD43") |>
#'     rep(each = 50),
#'   marker_2 = c("HLA-ABC", "B2M", "CD4", "CD8", "CD20", "CD19", "CD45", "CD43") |>
#'     rep(each = 50)
#' ) |>
#'   mutate(
#'     join_count_z = rnorm(n(), sd = 10)
#'   )
#'
#' example_data <- example_data |>
#'   mutate(sampleID = "ctrl") |>
#'   bind_rows(
#'     example_data |> mutate(join_count_z = join_count_z + 1) |>
#'       mutate(sampleID = "treatment")
#'   )
#'
#' # Compute statistics
#' dp_results <- DifferentialProximityAnalysis(
#'   example_data,
#'   contrast_column = "sampleID",
#'   reference = "ctrl",
#'   proximity_metric = "join_count_z",
#'   metric_type = "self"
#' )
#'
#' @export
#'
DifferentialProximityAnalysis.data.frame <- function(
  object,
  contrast_column,
  reference,
  targets = NULL,
  group_vars = NULL,
  proximity_metric = "log2_ratio",
  metric_type = c("all", "self", "co"),
  backend = c("dplyr", "data.table"),
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  min_cells_per_group = 10,
  verbose = TRUE,
  ...
) {
  .validate_da_input(
    object, contrast_column, reference, targets, group_vars,
    proximity_metric, min_cells_per_group, FALSE,
    cl = NULL, data_type = "proximity"
  )
  if (is.na(min_cells_per_group) || min_cells_per_group <= 0) {
    cli::cli_abort(
      c("x" = "{.arg min_cells_per_group} must be a positive number (> 0).")
    )
  }
  backend <- match.arg(backend, choices = c("dplyr", "data.table"))
  if (backend == "data.table") {
    expect_dtplyr()
  }
  metric_type <- match.arg(metric_type, choices = c("all", "self", "co"))
  object <- switch(metric_type,
    all = object,
    self = object |> filter(marker_1 == marker_2),
    co = object |> filter(marker_1 != marker_2)
  )
  if (nrow(object) == 0) {
    cli::cli_abort("No data found for the specified metric type.")
  }

  # Define targets as all groups except the reference if not specified
  targets <- targets %||% setdiff(unique(object[, contrast_column, drop = TRUE]), reference)

  # Check multiple choice args
  p_adjust_method <- match.arg(p_adjust_method,
    choices = c(
      "bonferroni", "holm", "hochberg",
      "hommel", "BH", "BY", "fdr"
    )
  )

  # Keep relevant columns
  object <- object |>
    select(any_of(c("marker_1", "marker_2", proximity_metric, "component", contrast_column, group_vars)))

  # Group data by marker pair, contrast column and optional group variables
  test_groups <- object |>
    group_by(pick(all_of(c("marker_1", "marker_2", contrast_column, group_vars))))

  if (min_cells_per_group > 0) {
    if (backend == "data.table") {
      test_groups <- dtplyr::lazy_dt(test_groups)
    }

    # Compress the massive table into a tiny count table
    valid_groups <- test_groups |>
      count(name = "cell_count") |>
      filter(cell_count >= min_cells_per_group) |>
      ungroup(!!sym(contrast_column)) |>
      mutate(n_groups = n(), ref_present = sum(!!sym(contrast_column) == reference)) |>
      filter(n_groups > 1, ref_present == 1) |>
      select(-n_groups, -ref_present)

    # Filter the massive original table using the valid combinations found above
    test_groups <- test_groups |>
      ungroup() |>
      semi_join(valid_groups, by = c("marker_1", "marker_2", contrast_column, group_vars)) |>
      as_tibble() |>
      group_by(pick(all_of(c("marker_1", "marker_2", contrast_column, group_vars))))
    if (nrow(test_groups) == 0) {
      cli::cli_abort(glue(
        "Found no groups with at least {min_cells_per_group} observations."
      ))
    }
  }

  # Compute group keys (drop marker columns; this method groups by marker pair
  # internally but iterates over contrast/group-var combinations only)
  test_groups_keys <- test_groups |>
    group_keys() |>
    select(-marker_1, -marker_2) |>
    filter(!!sym(contrast_column) %in% targets) |>
    distinct() |>
    mutate_all(as.character)

  if (nrow(test_groups_keys) == 0) {
    cli::cli_abort("Found no valid target data.")
  }

  # Print the test set up and start a progress bar
  .print_da_setup(
    keys = test_groups_keys,
    contrast_column = contrast_column,
    reference = reference,
    targets = targets,
    group_vars = group_vars,
    verbose = verbose,
    test_groups = test_groups,
    min_cells_per_group = min_cells_per_group
  )

  # Create a list to store the results in
  keys <- apply(test_groups_keys, 1, paste, collapse = "_")
  diff_prox_res <- rep(list(NULL), length(keys)) |>
    set_names(keys)

  # iterate over group pairs
  for (i in seq_len(nrow(test_groups_keys))) {
    # Fetch the current group keys
    # If group_vars=NULL, the key is simply one of targets.
    # Otherwise, the key is target + additional group_vars.
    cur_keys <- test_groups_keys[i, ]
    key <- as.character(cur_keys) |> paste(collapse = "_")
    groups_keep <- c(reference, cur_keys[, contrast_column, drop = TRUE])

    if (verbose && check_global_verbosity()) {
      cli_progress_update(
        set = i,
        status = .format_da_comparison_status(
          cur_keys = cur_keys,
          contrast_column = contrast_column,
          reference = reference,
          group_vars = group_vars
        )
      )
    }

    # Filter data to the current comparison (contrast + group_vars)
    cur_test_groups <- .filter_cur_group(
      test_groups = test_groups,
      cur_keys = cur_keys,
      contrast_column = contrast_column,
      reference = reference,
      group_vars = group_vars
    )

    # Group by marker pair and compute ranks
    cur_test_groups_rank <- cur_test_groups |>
      (\(.) {
        if (backend == "data.table") {
          dtplyr::lazy_dt(.) |>
            group_by(marker_1, marker_2) |>
            mutate(r = data.table::frank(-!!sym(proximity_metric), ties.method = "average"))
        } else {
          group_by(., marker_1, marker_2) |>
            mutate(r = rank(-!!sym(proximity_metric), ties.method = "average"))
        }
      })()

    # Compute the number of ties per marker pair
    cur_test_groups_ties <- cur_test_groups_rank |>
      group_by(marker_1, marker_2, r) |>
      summarize(nties = n(), .groups = "drop") |>
      group_by(marker_1, marker_2) |>
      summarize(nties_const = sum(nties^3 - nties), .groups = "drop") |>
      collect()

    # Compute the U statistic, Z score and p-value
    cur_test_groups_u <- cur_test_groups_rank |>
      group_by(!!sym(contrast_column), .add = TRUE) |>
      summarize(rank_sum = sum(r), n = n(), med = median(!!sym(proximity_metric)), .groups = "drop") |>
      pivot_wider(names_from = !!sym(contrast_column), values_from = c("rank_sum", "n", "med")) |>
      select(all_of(c(
        "marker_1", "marker_2",
        paste0("rank_sum_", groups_keep),
        paste0("n_", groups_keep),
        paste0("med_", groups_keep)
      ))) |>
      rename(rs_ref = 3L, rs_tgt = 4L, n_ref = 5L, n_tgt = 6L, med_ref = 7L, med_tgt = 8L) |>
      mutate(
        med_diff = med_tgt - med_ref,
        u = rs_tgt - n_tgt * (n_tgt + 1) / 2,
        auc = u / (n_ref * n_tgt)
      ) |>
      # Include the number of ties in the computation which
      # were computed in the previous step
      (\(.) {
        if (nrow(cur_test_groups_ties) > 0) {
          left_join(., cur_test_groups_ties, by = c("marker_1", "marker_2"))
        } else {
          mutate(., nties_const = 0)
        }
      })() |>
      mutate(z = u - (n_ref * n_tgt) / 2) |>
      mutate(
        sigma = sqrt(
          (n_ref * n_tgt / 12) *
            ((n_ref + n_tgt + 1) - nties_const /
              ((n_ref + n_tgt) * (n_ref + n_tgt - 1))
            )
        )
      ) |>
      # Formula for two.sided test which is currently the only option
      mutate(z = (z - sign(z) * 0.5) / sigma) |>
      mutate(p_val = 2 * pmin(pnorm(z), pnorm(z, lower.tail = FALSE))) |>
      collect() |>
      na.omit()

    # Tidy up the results
    diff_prox_res[[key]] <- tibble(
      target = cur_keys[, contrast_column, drop = TRUE],
      reference = reference,
      n_tgt = cur_test_groups_u$n_tgt,
      n_ref = cur_test_groups_u$n_ref,
      med_tgt = cur_test_groups_u$med_tgt,
      med_ref = cur_test_groups_u$med_ref,
      diff_median = cur_test_groups_u$med_diff,
      statistic = cur_test_groups_u$u,
      auc = cur_test_groups_u$auc,
      p = cur_test_groups_u$p_val,
      marker_1 = cur_test_groups_u$marker_1,
      marker_2 = cur_test_groups_u$marker_2
    ) |>
      filter(!is.na(n_ref), !is.na(n_tgt)) |>
      .append_group_vars(cur_keys = cur_keys, group_vars = group_vars)
  }

  if (verbose && check_global_verbosity()) {
    cli_progress_done()
  }

  if (!is.null(group_vars) && verbose) {
    cli::cli_alert_info("Adjusting p-values per group defined by group column{?s} {.str {group_vars}}.")
  }

  .finalize_da_results(diff_prox_res, p_adjust_method, group_vars)
}


#' @param group_data A data.frame with a column for the contrast and optional group variables.
#' The rownames of this data.frame should correspond to the columns names of the matrix `object`.
#' @param diff_threshold Minimum difference in proximity metric to consider a pair of groups for
#' testing. Default is 0.1. This parameter is only used when the `method` argument is set to "seurat".
#' @param only_pos If \code{TRUE}, test only marker pairs whose median is higher in the
#' target group than in the reference. Pairs that do not pass are not tested. Default is \code{FALSE}.
#' @param min_pct Minimum percentage of cells in either group that must express a marker pair for it
#' to be considered for testing. Default is 0. This parameter is only used when the `method` argument
#' is set to "seurat".
#' @param min_diff_pct Minimum difference in percentage of cells expressing a marker pair between the
#' two groups for it to be considered for testing. Default is -Inf. This parameter is only used when
#' the `method` argument is set to "seurat".
#'
#' @param group_data A tibble with a column for the contrast and optional group variables.
#' The rownames of this tibble should correspond to the columns names of the matrix `object`.
#'
#' @rdname DifferentialProximityAnalysis
#' @method DifferentialProximityAnalysis Matrix
#'
#' @export
#'
DifferentialProximityAnalysis.Matrix <- function(
  object,
  group_data,
  contrast_column,
  reference,
  targets = NULL,
  group_vars = NULL,
  proximity_metric = "log2_ratio",
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  diff_threshold = 0.01,
  min_pct = 0,
  min_diff_pct = -Inf,
  min_cells_per_group = 10,
  only_pos = FALSE,
  verbose = TRUE,
  ...
) {
  assert_single_value(only_pos, type = "bool")
  assert_single_value(min_cells_per_group, type = "numeric")
  if (is.na(min_cells_per_group) || min_cells_per_group <= 0) {
    cli::cli_abort(
      c("x" = "{.arg min_cells_per_group} must be a positive number (> 0).")
    )
  }
  # Validate group_data layout against the matrix `object`
  .validate_matrix_group_data(
    object = object,
    group_data = group_data,
    contrast_column = contrast_column,
    reference = reference,
    targets = targets,
    group_vars = group_vars
  )

  # Define targets as all groups except the reference if not specified
  targets <- targets %||% setdiff(unique(group_data[, contrast_column, drop = TRUE]), reference)

  # Check multiple choice args
  p_adjust_method <- match.arg(p_adjust_method,
    choices = c(
      "bonferroni", "holm", "hochberg",
      "hommel", "BH", "BY", "fdr"
    )
  )

  # Lift the component identifier into a column for downstream filtering
  group_data <- group_data |>
    as_tibble(rownames = "component")

  # Group data by contrast column and optional group variables, and compute keys
  grouped <- .compute_test_groups_keys(
    data = group_data,
    group_cols = c(contrast_column, group_vars),
    contrast_column = contrast_column,
    targets = targets
  )
  test_groups <- grouped$test_groups
  test_groups_keys <- grouped$keys

  # Print the test set up and start a progress bar
  .print_da_setup(
    keys = test_groups_keys,
    contrast_column = contrast_column,
    reference = reference,
    targets = targets,
    group_vars = group_vars,
    verbose = verbose,
    test_groups = test_groups,
    min_cells_per_group = min_cells_per_group
  )

  # Create a list to store the results in
  keys <- apply(test_groups_keys, 1, paste, collapse = "_")
  diff_prox_res <- rep(list(NULL), length(keys)) |>
    set_names(keys)
  n_passed_cells <- 0L

  # iterate over group pairs
  for (i in seq_len(nrow(test_groups_keys))) {
    # Fetch the current group keys
    cur_keys <- test_groups_keys[i, ]
    key <- as.character(cur_keys) |> paste(collapse = "_")

    # Filter data to the current comparison (contrast + group_vars)
    cur_test_groups <- .filter_cur_group(
      test_groups = test_groups,
      cur_keys = cur_keys,
      contrast_column = contrast_column,
      reference = reference,
      group_vars = group_vars
    )

    components_reference <- cur_test_groups |>
      filter(!!sym(contrast_column) == reference) |>
      pull(component)
    components_target <- cur_test_groups |>
      filter(!!sym(contrast_column) != reference) |>
      pull(component)

    n_comps_reference <- length(components_reference)
    n_comps_target <- length(components_target)
    if (verbose && check_global_verbosity()) {
      cli_progress_update(
        set = i,
        status = .format_da_comparison_status(
          cur_keys = cur_keys,
          contrast_column = contrast_column,
          reference = reference,
          group_vars = group_vars,
          n_tgt = n_comps_target,
          n_ref = n_comps_reference
        )
      )
    }
    if (n_comps_reference < min_cells_per_group || n_comps_target < min_cells_per_group) {
      target <- cur_keys[, contrast_column, drop = TRUE]
      if (ncol(cur_keys) > 1) {
        cur_reference <- paste0(reference, " (", paste(cur_keys[, group_vars, drop = TRUE], collapse = ", "), ")")
        cur_target <- paste0(target, " (", paste(cur_keys[, group_vars, drop = TRUE], collapse = ", "), ")")
      } else {
        cur_reference <- reference
        cur_target <- target
      }
      if (n_comps_reference < min_cells_per_group && n_comps_target < min_cells_per_group) {
        cli::cli_warn(
          "Skipping {cur_target} vs {cur_reference} because both
           groups have fewer than {.val {min_cells_per_group}} cells."
        )
      } else if (n_comps_reference < min_cells_per_group && n_comps_target >= min_cells_per_group) {
        cli::cli_warn(
          "Skipping {cur_target} vs {cur_reference} because the reference
           group has fewer than {.val {min_cells_per_group}} cells."
        )
      } else if (n_comps_target < min_cells_per_group && n_comps_reference >= min_cells_per_group) {
        cli::cli_warn(
          "Skipping {cur_target} vs {cur_reference} because the target group
           has fewer than {.val {min_cells_per_group}} cells."
        )
      }
      next
    }

    n_passed_cells <- n_passed_cells + 1L
    de_results <- .wilcox_de_test(
      object,
      cells_1 = components_target,
      cells_2 = components_reference,
      min_diff = diff_threshold,
      min_pct = min_pct,
      min_diff_pct = min_diff_pct,
      only_pos = only_pos
    )

    # Filters in .wilcox_de_test() can leave no pairs. Skip the comparison
    # instead of building a zero-row tibble with recycled columns.
    if (nrow(de_results) == 0) {
      next
    }

    # Tidy up the results
    diff_prox_res[[key]] <- tibble(
      target = cur_keys[, contrast_column, drop = TRUE],
      reference = reference,
      pct_tgt = de_results$pct_1,
      pct_ref = de_results$pct_2,
      med_tgt = de_results$med_1,
      med_ref = de_results$med_2,
      diff_median = de_results$difference,
      p = de_results$p_val,
      pair = rownames(de_results)
    ) |>
      .append_group_vars(cur_keys = cur_keys, group_vars = group_vars)
  }

  if (verbose && check_global_verbosity()) {
    cli_progress_done()
  }

  valid_results <- sapply(diff_prox_res, function(x) !is.null(x))
  if (!any(valid_results)) {
    .abort_no_valid_da_results(
      n_total = nrow(test_groups_keys),
      n_passed_cells = n_passed_cells,
      min_cells_per_group = min_cells_per_group,
      contrast_column = contrast_column,
      group_vars = group_vars,
      diff_threshold = diff_threshold,
      min_pct = min_pct,
      min_diff_pct = min_diff_pct,
      only_pos = only_pos
    )
  }

  .finalize_da_results(diff_prox_res, p_adjust_method, group_vars) |>
    tidyr::separate(pair, into = c("marker_1", "marker_2"), sep = ":")
}


#' @param lazy If TRUE, the proximity scores will be loaded lazily and filtered using the
#' `duckdb` backend. Scores are fetched with \code{\link{ProximityScores}(object, lazy = TRUE)}.
#' Set this when the \code{Seurat} object was created with \code{load_proximity_scores = FALSE}.
#' @param min_exp_join_count Minimum number of join counts required for a marker pair to be
#' included in the analysis. Dropped protein pairs (those with fewer than `min_exp_join_count` counts)
#' will be treated as missing entries. With `method = "seurat"`, these are treated as having a
#' proximity score of 0.
#' @param method One of "seurat" or "legacy". The former uses the Seurat framework for
#' differential testing, while the latter uses a custom implementation. The main difference
#' between the two methods is that missing observations are handled differently. With
#' the "seurat" method, all missing values are set to 0, while the "legacy" method ignores
#' missing values. For the former, this means that the number of observations per group is
#' always equal to the number of cells in that group, while for the latter ("legacy"), the
#' number of observations per group can be less than the number of cells in that group. Ignoring
#' missing values can lead to confusing estimates of group statistics and misses comparisons
#' where one of the two test groups have no observations.
#' The "legacy" method is provided for backward compatibility and may be removed in future versions.
#'
#' @param assay Name of assay to use
#'
#' @rdname DifferentialProximityAnalysis
#' @method DifferentialProximityAnalysis Seurat
#'
#' @export
#'
DifferentialProximityAnalysis.Seurat <- function(
  object,
  contrast_column,
  reference,
  targets = NULL,
  assay = NULL,
  group_vars = NULL,
  lazy = FALSE,
  min_exp_join_count = 0,
  min_cells_per_group = 10,
  diff_threshold = 0.01,
  min_pct = 0,
  min_diff_pct = -Inf,
  proximity_metric = "log2_ratio",
  metric_type = c("all", "self", "co"),
  backend = c("dplyr", "data.table"),
  method = c("seurat", "legacy"),
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
) {
  # Validate input parameters
  assert_col_in_data(contrast_column, object[[]])
  assert_single_value(lazy, type = "bool")

  assert_single_value(min_cells_per_group, type = "numeric")
  if (is.na(min_cells_per_group) || min_cells_per_group <= 0) {
    cli::cli_abort(
      c("x" = "{.arg min_cells_per_group} must be a positive number (> 0).")
    )
  }

  # Use default assay if assay = NULL
  assay <- .validate_or_set_assay(object, assay)

  method <- match.arg(method, choices = c("seurat", "legacy"))

  # Add group data
  group_data <- .select_group_data(object[[]], contrast_column, group_vars)

  # Fast single-pass DuckDB path when lazy = TRUE and method = "seurat"
  if (isTRUE(lazy) && identical(method, "seurat")) {
    fs_map <- tryCatch(FSMap(object, assay = assay), error = function(e) NULL)
    if (!inherits(fs_map, "pxg_fs_map")) {
      return(.run_lazy_duckdb_dpa_seurat(
        object = object,
        group_data_tbl = group_data,
        contrast_column = contrast_column,
        reference = reference,
        targets = targets,
        assay = assay,
        group_vars = group_vars,
        min_exp_join_count = min_exp_join_count,
        min_cells_per_group = min_cells_per_group,
        diff_threshold = diff_threshold,
        min_pct = min_pct,
        min_diff_pct = min_diff_pct,
        proximity_metric = proximity_metric,
        metric_type = metric_type,
        p_adjust_method = p_adjust_method,
        verbose = verbose,
        ...
      ))
    }
  }

  active_targets <- targets %||% setdiff(unique(group_data[, contrast_column, drop = TRUE]), reference)
  cells_keep <- group_data |>
    filter(.data[[contrast_column]] %in% c(reference, active_targets)) |>
    pull(component)

  # Fetch proximity scores
  proximity_data <- .load_proximity_for_da(
    object = object,
    assay = assay,
    lazy = lazy,
    min_exp_join_count = min_exp_join_count,
    proximity_metric = proximity_metric,
    metric_type = metric_type,
    method = method,
    cells_keep = if (identical(method, "legacy")) cells_keep else NULL
  )

  if (method == "seurat") {
    group_data <- group_data |>
      column_to_rownames("component")
  }

  # Add group data to colocalization table
  if (method == "legacy") {
    proximity_data <- proximity_data |>
      left_join(y = group_data, copy = TRUE, by = "component") |>
      collect()
  }

  # Run differential proximity analysis
  proximity_test_results <- switch(method,
    legacy = DifferentialProximityAnalysis(
      object = proximity_data,
      contrast_column = contrast_column,
      reference = reference,
      targets = targets,
      group_vars = group_vars,
      proximity_metric = proximity_metric,
      backend = backend,
      p_adjust_method = p_adjust_method,
      min_cells_per_group = min_cells_per_group,
      verbose = verbose
    ),
    seurat = DifferentialProximityAnalysis(
      object = proximity_data,
      group_data = group_data,
      contrast_column = contrast_column,
      reference = reference,
      targets = targets,
      group_vars = group_vars,
      proximity_metric = proximity_metric,
      diff_threshold = diff_threshold,
      p_adjust_method = p_adjust_method,
      min_cells_per_group = min_cells_per_group,
      min_pct = min_pct,
      min_diff_pct = min_diff_pct,
      verbose = verbose,
      ...
    )
  )

  return(proximity_test_results)
}

#' Compute the median difference and the fraction of non-zero scores in two groups.
#'
#' @param object A matrix where rows are features and columns are cells.
#' @param cells_1 Cell names for the target group.
#' @param cells_2 Cell names for the reference group.
#'
#' @return A data.frame with columns \code{difference}, \code{med_1},
#'   \code{med_2}, \code{pct_1}, and \code{pct_2}.
#'
#' @noRd
.median_difference <- function(object, cells_1, cells_2) {
  pct_1 <- round(
    Matrix::rowSums(object[, cells_1, drop = FALSE] != 0) /
      length(x = cells_1),
    digits = 3
  )
  pct_2 <- round(
    Matrix::rowSums(object[, cells_2, drop = FALSE] != 0) /
      length(x = cells_2),
    digits = 3
  )
  med_1 <- sparseMatrixStats::rowMedians(object[, cells_1, drop = FALSE])
  med_2 <- sparseMatrixStats::rowMedians(object[, cells_2, drop = FALSE])
  med_diff <- med_1 - med_2
  fc_results <- as.data.frame(x = cbind(med_diff, med_1, med_2, pct_1, pct_2))
  colnames(fc_results) <- c("difference", "med_1", "med_2", "pct_1", "pct_2")
  return(fc_results)
}


#' Compute row-wise tie constants for a sparse matrix.
#'
#' Counts the large zero tie block from sparsity and only inspects stored
#' non-zero values to count additional ties.
#'
#' @param object A sparse matrix with features in rows.
#'
#' @return A numeric vector of tie constants, one per row.
#'
#' @noRd
.sparse_tie_constant <- function(object) {
  object <- methods::as(object, "RsparseMatrix")
  n_rows <- nrow(object)
  n_cols <- ncol(object)
  out <- numeric(n_rows)
  p <- object@p
  x <- object@x

  for (i in seq_len(n_rows)) {
    start <- p[i] + 1L
    end <- p[i + 1L]
    if (start <= end) {
      vals <- x[start:end]
      n_nonzero <- length(vals)
    } else {
      vals <- numeric()
      n_nonzero <- 0L
    }

    n_zero <- as.numeric(n_cols - n_nonzero)
    tie_zero <- if (n_zero > 1) n_zero^3 - n_zero else 0
    tie_nonzero <- 0
    if (n_nonzero > 1L) {
      runs <- as.numeric(rle(sort(vals))$lengths)
      runs <- runs[runs > 1]
      if (length(runs) > 0L) {
        tie_nonzero <- sum(runs^3 - runs)
      }
    }
    out[i] <- tie_zero + tie_nonzero
  }

  out
}


#' Compute row-wise Wilcoxon p-values with a sparse-exact tie correction.
#'
#' @param data_use A matrix where rows are features and columns are ordered as
#'   the first group followed by the second group.
#' @param n_cells_1 Number of cells in the first group.
#'
#' @return A numeric vector of p-values, one per row of
#'   \\code{data_use}.
#'
#' @noRd
.wilcox_p_sparse_exact <- function(data_use, n_cells_1) {
  expect_matrixStats()
  data_dense <- as.matrix(data_use)
  n_cells_1 <- as.numeric(n_cells_1)
  n_cells_2 <- as.numeric(ncol(data_use)) - n_cells_1

  ranks <- matrixStats::rowRanks(
    data_dense,
    ties.method = "average",
    preserveShape = TRUE
  )
  rank_sum_1 <- rowSums(ranks[, seq_len(n_cells_1), drop = FALSE])
  u <- rank_sum_1 - n_cells_1 * (n_cells_1 + 1) / 2

  tie_const <- .sparse_tie_constant(data_use)
  sigma <- sqrt(
    (n_cells_1 * n_cells_2 / 12) *
      ((n_cells_1 + n_cells_2 + 1) - tie_const /
        ((n_cells_1 + n_cells_2) * (n_cells_1 + n_cells_2 - 1))
      )
  )
  z <- u - (n_cells_1 * n_cells_2) / 2

  p_val <- rep(1, length(z))
  ok <- is.finite(sigma) & sigma > 0
  z[ok] <- (z[ok] - sign(z[ok]) * 0.5) / sigma[ok]
  p_val[ok] <- 2 * pmin(
    pnorm(z[ok]),
    pnorm(z[ok], lower.tail = FALSE)
  )
  p_val
}


#' Compute row-wise Wilcoxon p-values with the legacy fallback path.
#'
#' Uses the existing \\pkg{limma} loop when available, otherwise falls back to
#' \\code{stats::wilcox.test()}.
#'
#' @param data_use A matrix where rows are features and columns are ordered as
#'   the first group followed by the second group.
#' @param cells_1 Cell names for the first group.
#' @param cells_2 Cell names for the second group.
#'
#' @return A numeric vector of p-values, one per row of \\code{data_use}.
#'
#' @noRd
.wilcox_p_fallback <- function(data_use, cells_1, cells_2) {
  if (requireNamespace("limma", quietly = TRUE)) {
    j <- seq_along(cells_1)
    return(sapply(seq_len(nrow(data_use)), function(i) {
      min(2 * min(limma::rankSumTestWithCorrelation(
        index = j,
        statistics = data_use[i, ]
      )), 1)
    }))
  }

  rlang::inform(
    c(
      "i" = "For a faster implementation of the Wilcoxon Rank Sum Test,
      please install the limma package: pak::pak('limma')",
      "v" = "This message will only appear once per R session."
    ),
    .frequency = "once",
    .frequency_id = "limma_not_installed"
  )
  sapply(seq_len(nrow(data_use)), function(i) {
    x1 <- as.numeric(data_use[i, cells_1])
    x2 <- as.numeric(data_use[i, cells_2])
    suppressWarnings(stats::wilcox.test(x1, x2, exact = FALSE)$p.value)
  })
}

#' Utility function to perform Wilcoxon rank-sum test on two groups of cells.
#'
#' This function is adapted from Seurat::FindMarkers, mainly to avoid overhead
#' costs associated with pbsapply.
#'
#' @param data_use A matrix where rows are features (e.g. marker pairs) and columns are cells.
#' @param cells_1 A vector of cell names for group 1
#' @param cells_2 A vector of cell names for group 2
#' @param min_diff Minimum difference in median expression between the two groups to consider a feature
#' @param min_pct Minimum percentage of cells expressing a feature in either group to consider it
#' @param min_diff_pct Minimum difference in percentage of cells expressing a feature between the two
#' @param only_pos If \code{TRUE}, test only features with a positive median difference.
#'
#' @return A data.frame of Wilcoxon results, or an empty data.frame when no feature passes.
#'
#' @noRd
.wilcox_de_test <- function(
  data_use,
  cells_1,
  cells_2,
  min_diff = 0.01,
  min_pct = 0.01,
  min_diff_pct = -Inf,
  only_pos = FALSE
) {
  expect_sparseMatrixStats()

  fc_results <- .median_difference(
    object = data_use,
    cells_1 = cells_1,
    cells_2 = cells_2
  )

  alpha_min <- pmax(fc_results$pct_1, fc_results$pct_2)
  names(alpha_min) <- rownames(fc_results)
  features <- names(which(alpha_min >= min_pct))

  if (length(features) == 0) {
    cli::cli_warn(
      "No features pass min_pct threshold; returning empty data.frame"
    )
    return(fc_results[features, ])
  }

  alpha_diff <- alpha_min - pmin(fc_results$pct_1, fc_results$pct_2)
  features <- names(which(alpha_min >= min_pct & alpha_diff >= min_diff_pct))
  if (length(features) == 0) {
    cli::cli_warn(
      "No features pass min_diff_pct threshold; returning empty data.frame"
    )
    return(fc_results[features, ])
  }

  total_diff <- fc_results[, 1]
  names(total_diff) <- rownames(fc_results)
  if (only_pos) {
    # cells_1 is the target, so a positive difference is higher in the target.
    features_diff <- names(which(total_diff > 0 & total_diff >= min_diff))
  } else {
    features_diff <- names(which(abs(total_diff) >= min_diff))
  }

  features <- intersect(features, features_diff)
  if (length(features) == 0) {
    if (only_pos) {
      cli::cli_warn(
        "No features have a higher median in the target group; returning empty data.frame"
      )
    } else {
      cli::cli_warn("No features pass min_diff threshold; returning empty data.frame")
    }
    return(fc_results[features, ])
  }

  data_use <- data_use[features, c(cells_1, cells_2), drop = FALSE]

  p_val <- tryCatch(
    .wilcox_p_sparse_exact(data_use, n_cells_1 = length(cells_1)),
    error = function(e) NULL
  )
  if (is.null(p_val)) {
    p_val <- .wilcox_p_fallback(data_use, cells_1 = cells_1, cells_2 = cells_2)
  }

  de_results <- cbind(
    data.frame(p_val, row.names = rownames(data_use)),
    fc_results[rownames(data_use), ]
  )
  return(de_results)
}


#' @param assay Name of the assay to use.
#' @param lazy If \code{TRUE}, proximity scores are fetched with
#' \code{\link{ProximityScores}(object, lazy = TRUE)} and filtered with the
#' \code{duckdb} backend before testing. Set this when the \code{Seurat}
#' object was created with \code{load_proximity_scores = FALSE} and the
#' proximity scores were not stored in the object.
#' @param min_exp_join_count Minimum expected join count for a marker pair to
#' be included. Pairs below the threshold are treated as missing and set to 0.
#'
#' @rdname FindAllProximityMarkers
#' @method FindAllProximityMarkers Seurat
#'
#' @export
#'
FindAllProximityMarkers.Seurat <- function(
  object,
  group_by,
  idents = NULL,
  assay = NULL,
  lazy = FALSE,
  min_exp_join_count = 0,
  min_cells_per_group = 10,
  diff_threshold = 0.01,
  min_pct = 0,
  min_diff_pct = -Inf,
  only_pos = FALSE,
  proximity_metric = "log2_ratio",
  metric_type = c("all", "self", "co"),
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
) {
  call <- caller_env()
  .reject_find_all_proximity_args(..., call = call)
  assert_single_value(group_by, type = "string", call = call)
  assert_col_in_data(group_by, object[[]], call = call)
  assert_col_class(
    group_by,
    object[[]],
    classes = c("character", "factor"),
    call = call
  )
  .validate_find_all_common_args(
    min_cells_per_group = min_cells_per_group,
    diff_threshold = diff_threshold,
    min_pct = min_pct,
    min_diff_pct = min_diff_pct,
    only_pos = only_pos,
    proximity_metric = proximity_metric,
    verbose = verbose,
    call = call
  )

  assay <- .validate_or_set_assay(object, assay, call = call)
  assert_single_value(lazy, type = "bool", call = call)
  metric_type <- match.arg(metric_type, choices = c("all", "self", "co"))
  p_adjust_method <- match.arg(
    p_adjust_method,
    choices = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr")
  )

  meta <- object[[]]
  label_chr <- as.character(meta[[group_by]])
  names(label_chr) <- rownames(meta)
  idents <- .resolve_one_vs_rest_idents(
    labels = meta[[group_by]],
    idents = idents,
    group_by = group_by,
    call = call
  )
  rest_label <- .make_rest_label(label_chr)
  .warn_missing_group_labels(label_chr, group_by)

  if (isTRUE(lazy)) {
    fs_map <- tryCatch(FSMap(object, assay = assay), error = function(e) NULL)
    if (!inherits(fs_map, "pxg_fs_map")) {
      return(.run_lazy_duckdb_one_vs_rest_seurat(
        object = object,
        assay = assay,
        label_chr = label_chr,
        idents = idents,
        rest_label = rest_label,
        group_by = group_by,
        min_exp_join_count = min_exp_join_count,
        min_cells_per_group = min_cells_per_group,
        diff_threshold = diff_threshold,
        min_pct = min_pct,
        min_diff_pct = min_diff_pct,
        only_pos = only_pos,
        proximity_metric = proximity_metric,
        metric_type = metric_type,
        p_adjust_method = p_adjust_method,
        verbose = verbose,
        call = call
      ))
    }
  }

  proximity_matrix <- .load_proximity_for_da(
    object = object,
    assay = assay,
    lazy = lazy,
    min_exp_join_count = min_exp_join_count,
    proximity_metric = proximity_metric,
    metric_type = metric_type,
    method = "seurat"
  )

  .run_matrix_one_vs_rest(
    proximity_matrix = proximity_matrix,
    label_chr = label_chr,
    idents = idents,
    rest_label = rest_label,
    group_by = group_by,
    min_cells_per_group = min_cells_per_group,
    diff_threshold = diff_threshold,
    min_pct = min_pct,
    min_diff_pct = min_diff_pct,
    only_pos = only_pos,
    proximity_metric = proximity_metric,
    p_adjust_method = p_adjust_method,
    verbose = verbose,
    call = call
  )
}


#' @rdname FindAllProximityMarkers
#' @method FindAllProximityMarkers data.frame
#'
#' @examples
#' library(dplyr)
#' set.seed(1)
#' example_data <- tidyr::expand_grid(
#'   marker_1 = c("CD19", "CD3E"),
#'   marker_2 = c("CD19", "CD3E"),
#'   component = paste0("cell", seq_len(12))
#' ) |>
#'   mutate(
#'     join_count_z = rnorm(n()),
#'     seurat_clusters = dplyr::case_when(
#'       component %in% paste0("cell", 1:4) ~ "0",
#'       component %in% paste0("cell", 5:8) ~ "1",
#'       TRUE ~ "2"
#'     )
#'   )
#'
#' # Compare each cluster with the cells in every other cluster
#' FindAllProximityMarkers(
#'   example_data,
#'   group_by = "seurat_clusters",
#'   proximity_metric = "join_count_z",
#'   min_cells_per_group = 4,
#'   verbose = FALSE
#' )
#'
#' @export
#'
FindAllProximityMarkers.data.frame <- function(
  object,
  group_by,
  idents = NULL,
  min_cells_per_group = 10,
  diff_threshold = 0.01,
  min_pct = 0,
  min_diff_pct = -Inf,
  only_pos = FALSE,
  proximity_metric = "log2_ratio",
  metric_type = c("all", "self", "co"),
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
) {
  call <- caller_env()
  .reject_find_all_proximity_args(..., call = call)
  assert_single_value(group_by, type = "string", call = call)
  .validate_find_all_common_args(
    min_cells_per_group = min_cells_per_group,
    diff_threshold = diff_threshold,
    min_pct = min_pct,
    min_diff_pct = min_diff_pct,
    only_pos = only_pos,
    proximity_metric = proximity_metric,
    verbose = verbose,
    call = call
  )
  metric_type <- match.arg(metric_type, choices = c("all", "self", "co"))
  p_adjust_method <- match.arg(
    p_adjust_method,
    choices = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr")
  )

  if (group_by %in% c("marker_1", "marker_2", "component", proximity_metric)) {
    cli::cli_abort(
      "{.arg group_by} must be a grouping column, not {.val {group_by}}.",
      call = call
    )
  }

  # Copy so the caller's table is left unchanged, including data.table inputs.
  object <- as_tibble(object) |> ungroup()
  assert_col_in_data("component", object, call = call)
  assert_col_in_data("marker_1", object, call = call)
  assert_col_in_data("marker_2", object, call = call)
  assert_col_in_data(group_by, object, call = call)
  assert_col_class(group_by, object, classes = c("character", "factor"), call = call)
  assert_col_in_data(proximity_metric, object, call = call)
  assert_col_class(proximity_metric, object, classes = "numeric", call = call)

  n_missing_component <- sum(is.na(object$component))
  if (n_missing_component > 0) {
    cli::cli_warn(
      "Dropping {n_missing_component} row{?s} with missing component values."
    )
    object <- object |> filter(!is.na(component))
  }

  labels_by_cell <- object |>
    distinct(component, !!sym(group_by))
  if (anyDuplicated(labels_by_cell$component)) {
    cli::cli_abort(
      "Each component must have a single value in {.val {group_by}}.",
      call = call
    )
  }

  label_chr <- as.character(labels_by_cell[[group_by]])
  names(label_chr) <- as.character(labels_by_cell$component)
  idents <- .resolve_one_vs_rest_idents(
    labels = labels_by_cell[[group_by]],
    idents = idents,
    group_by = group_by,
    call = call
  )
  rest_label <- .make_rest_label(label_chr)
  .warn_missing_group_labels(label_chr, group_by)

  scored <- switch(metric_type,
    all = object,
    self = object |> filter(marker_1 == marker_2),
    co = object |> filter(marker_1 != marker_2)
  )
  if (nrow(scored) == 0) {
    cli::cli_abort(
      "No data found for the specified metric type.",
      call = call
    )
  }
  proximity_matrix <- scored |>
    ProximityScoresToAssay(values_from = proximity_metric) |>
    .pad_proximity_matrix(cells = names(label_chr))

  .run_matrix_one_vs_rest(
    proximity_matrix = proximity_matrix,
    label_chr = label_chr,
    idents = idents,
    rest_label = rest_label,
    group_by = group_by,
    min_cells_per_group = min_cells_per_group,
    diff_threshold = diff_threshold,
    min_pct = min_pct,
    min_diff_pct = min_diff_pct,
    only_pos = only_pos,
    proximity_metric = proximity_metric,
    p_adjust_method = p_adjust_method,
    verbose = verbose,
    call = call
  )
}


#' Load and filter proximity scores once for differential testing.
#'
#' @param object A Seurat object.
#' @param assay Name of the assay to read.
#' @param lazy If \code{TRUE}, fetch scores from the PXL file.
#' @param min_exp_join_count Minimum expected join count for a marker pair.
#' @param proximity_metric Numeric column used when building the matrix.
#' @param metric_type One of \code{"all"}, \code{"self"}, or \code{"co"}.
#' @param method \code{"seurat"} returns a matrix and \code{"legacy"} returns a table.
#' @param cells_keep Optional character vector of cell IDs to retain before
#'   materializing the lazy table.
#'
#' @return A sparse matrix of marker pairs by cells, or a proximity table for the legacy method.
#'
#' @noRd
.load_proximity_for_da <- function(
  object,
  assay,
  lazy,
  min_exp_join_count,
  proximity_metric,
  metric_type,
  method,
  cells_keep = NULL
) {
  method <- match.arg(method, choices = c("seurat", "legacy"))
  assert_single_value(lazy, type = "bool")
  if (!lazy) {
    proximity_slot <- slot(object[[assay]], name = "proximity")
    if (is.null(proximity_slot) || ncol(proximity_slot) == 0) {
      cli::cli_abort(
        c(
          "x" = "Proximity scores are missing from the {.cls Seurat} object.",
          "i" = paste(
            "Create the object with {.code load_proximity_scores = TRUE},",
            "or set {.code lazy = TRUE} to fetch them with {.fn ProximityScores}."
          )
        ),
        call = caller_env()
      )
    }
  }
  proximity_data <- ProximityScores(object, assay = assay, lazy = lazy)
  metric_type <- match.arg(metric_type, choices = c("all", "self", "co"))
  proximity_data <- switch(metric_type,
    all = proximity_data,
    self = proximity_data |> filter(marker_1 == marker_2),
    co = proximity_data |> filter(marker_1 != marker_2)
  ) |>
    filter(join_count_expected_mean >= min_exp_join_count)
  if (!is.null(cells_keep)) {
    proximity_data <- proximity_data |> filter(component %in% !!cells_keep)
  }
  proximity_data <- proximity_data |> compute()
  if (method == "seurat") {
    proximity_data <- proximity_data |>
      ProximityScoresToAssay(values_from = proximity_metric) |>
      .pad_proximity_matrix(cells = colnames(object))
  }
  proximity_data
}


#' Compute exact Wilcoxon rank-sum p-values from sorted non-zero target and reference lists.
#'
#' Computes tie-corrected two-sided Wilcoxon rank-sum p-values (with continuity
#' correction) identical to \code{.wilcox_p_sparse_exact()} by shifting positive
#' ranks by the implicit zero block \code{n_zero_tot} and accounting for both
#' the zero tie block (\code{n_zero_tot^3 - n_zero_tot}) and non-zero ties.
#'
#' @param neg_tgt List of sorted negative numeric vectors in the target group.
#' @param pos_tgt List of sorted positive numeric vectors in the target group.
#' @param neg_ref List of sorted negative numeric vectors in the reference group.
#' @param pos_ref List of sorted positive numeric vectors in the reference group.
#' @param n_tgt Integer vector of total cells in the target group per row.
#' @param n_ref Integer vector of total cells in the reference group per row.
#'
#' @return A numeric vector of two-sided Wilcoxon rank-sum p-values.
#'
#' @noRd
.wilcox_from_sorted_lists <- function(
  neg_tgt,
  pos_tgt,
  neg_ref,
  pos_ref,
  n_tgt,
  n_ref
) {
  n_rows <- length(n_tgt)
  p_val <- rep(1.0, n_rows)
  for (r in seq_len(n_rows)) {
    nt <- as.numeric(n_tgt[[r]])
    nr <- as.numeric(n_ref[[r]])
    n_tot <- nt + nr
    v_neg_t <- neg_tgt[[r]] %||% numeric()
    v_pos_t <- pos_tgt[[r]] %||% numeric()
    v_neg_r <- neg_ref[[r]] %||% numeric()
    v_pos_r <- pos_ref[[r]] %||% numeric()

    k_t <- length(v_neg_t) + length(v_pos_t)
    k_r <- length(v_neg_r) + length(v_pos_r)
    n_zero_t <- nt - k_t
    n_zero_r <- nr - k_r
    n_zero_tot <- n_zero_t + n_zero_r
    n_neg_tot <- length(v_neg_t) + length(v_neg_r)

    vals_nz <- c(v_neg_t, v_pos_t, v_neg_r, v_pos_r)
    n_nz <- length(vals_nz)
    if (n_nz == 0L) next

    rk_nz <- rank(vals_nz, ties.method = "average")
    pos_mask <- vals_nz > 0
    rk_nz[pos_mask] <- rk_nz[pos_mask] + n_zero_tot

    sum_rk_nz_t <- if (k_t > 0L) sum(rk_nz[seq_len(k_t)]) else 0.0
    avg_rk_zero <- n_neg_tot + (n_zero_tot + 1.0) / 2.0
    rank_sum_1 <- sum_rk_nz_t + n_zero_t * avg_rk_zero
    u <- rank_sum_1 - nt * (nt + 1.0) / 2.0

    tie_zero <- if (n_zero_tot > 1.0) n_zero_tot^3 - n_zero_tot else 0.0
    tie_nonzero <- 0.0
    if (n_nz > 1L) {
      runs <- as.numeric(rle(sort(vals_nz))$lengths)
      runs <- runs[runs > 1.0]
      if (length(runs) > 0L) {
        tie_nonzero <- sum(runs^3 - runs)
      }
    }
    tie_const <- tie_zero + tie_nonzero
    sigma <- sqrt(
      (nt * nr / 12.0) * ((n_tot + 1.0) - tie_const / (n_tot * (n_tot - 1.0)))
    )
    if (is.finite(sigma) && sigma > 0.0) {
      z <- u - (nt * nr) / 2.0
      z <- (z - sign(z) * 0.5) / sigma
      p_val[[r]] <- 2.0 * min(stats::pnorm(z), stats::pnorm(z, lower.tail = FALSE))
    }
  }
  p_val
}


#' Execute single-pass zero-padded Wilcoxon differential proximity comparisons in DuckDB.
#'
#' @param lazy_prox A \code{tbl_lazy} proximity table.
#' @param cells_tbl A \code{data.frame} mapping \code{component} to \code{grp_key} and \code{side} (\code{"tgt"} or \code{"ref"}).
#' @param group_sizes_tbl A \code{data.frame} with columns \code{grp_key}, \code{n_tgt}, and \code{n_ref}.
#' @param proximity_metric Name of the numeric proximity column to test.
#' @param diff_threshold Minimum absolute median difference required to test a pair.
#' @param min_pct Minimum fraction of cells with a non-zero score in either group.
#' @param min_diff_pct Minimum difference in that fraction between the two groups.
#' @param only_pos If \code{TRUE}, test only pairs with a higher median in the target.
#' @param call Calling environment used in error messages.
#'
#' @return A \code{data.frame} of filtered marker-pair summary statistics and Wilcoxon p-values across all comparisons.
#'
#' @noRd
.run_duckdb_wilcox_comparisons <- function(
  lazy_prox,
  cells_tbl,
  group_sizes_tbl,
  proximity_metric,
  diff_threshold,
  min_pct,
  min_diff_pct,
  only_pos,
  call = caller_env()
) {
  assert_single_value(proximity_metric, type = "string", call = call)
  assert_col_in_data(proximity_metric, lazy_prox, call = call)
  assert_col_class(proximity_metric, lazy_prox, classes = "numeric", call = call)

  con <- lazy_prox$src$con
  base_sql <- as.character(dbplyr::sql_render(lazy_prox))
  val_col <- as.character(DBI::dbQuoteIdentifier(con, proximity_metric))

  cells_tbl_name <- paste0("__pxl_da_cells_", .generate_random_string(8))
  groups_tbl_name <- paste0("__pxl_da_groups_", .generate_random_string(8))
  DBI::dbWriteTable(con, cells_tbl_name, cells_tbl, overwrite = TRUE, temporary = TRUE)
  DBI::dbWriteTable(con, groups_tbl_name, group_sizes_tbl, overwrite = TRUE, temporary = TRUE)
  on.exit(
    {
      if (DBI::dbIsValid(con)) {
        try(DBI::dbRemoveTable(con, cells_tbl_name), silent = TRUE)
        try(DBI::dbRemoveTable(con, groups_tbl_name), silent = TRUE)
      }
    },
    add = TRUE
  )

  sql_query <- glue::glue(
    "
    WITH prox_base AS (
      SELECT row_number() OVER () AS __rn, *
      FROM ({base_sql}) __q
    ),
    pair_order AS (
      SELECT marker_1, marker_2, MIN(__rn) AS first_rn
      FROM prox_base
      GROUP BY marker_1, marker_2
      HAVING COUNT(*) FILTER (WHERE {val_col} IS NOT NULL AND {val_col} != 0) > 0
    ),
    nz AS (
      SELECT
        c.grp_key,
        c.side,
        p.marker_1,
        p.marker_2,
        CAST(p.{val_col} AS DOUBLE) AS val
      FROM prox_base p
      INNER JOIN {cells_tbl_name} c ON p.component = c.component
      WHERE p.{val_col} IS NOT NULL AND p.{val_col} != 0
    ),
    pair_counts AS (
      SELECT
        g.grp_key,
        po.marker_1,
        po.marker_2,
        po.first_rn,
        g.n_tgt,
        g.n_ref,
        COALESCE(agg.k_tgt, 0) AS k_tgt,
        COALESCE(agg.k_ref, 0) AS k_ref,
        COALESCE(agg.nn_tgt, 0) AS nn_tgt,
        COALESCE(agg.np_tgt, 0) AS np_tgt,
        COALESCE(agg.nn_ref, 0) AS nn_ref,
        COALESCE(agg.np_ref, 0) AS np_ref,
        agg.neg_tgt,
        agg.pos_tgt,
        agg.neg_ref,
        agg.pos_ref
      FROM {groups_tbl_name} g
      CROSS JOIN pair_order po
      LEFT JOIN (
        SELECT
          grp_key,
          marker_1,
          marker_2,
          CAST(COUNT(*) FILTER (WHERE side = 'tgt') AS INTEGER) AS k_tgt,
          CAST(COUNT(*) FILTER (WHERE side = 'ref') AS INTEGER) AS k_ref,
          CAST(COUNT(*) FILTER (WHERE side = 'tgt' AND val < 0) AS INTEGER) AS nn_tgt,
          CAST(COUNT(*) FILTER (WHERE side = 'tgt' AND val > 0) AS INTEGER) AS np_tgt,
          CAST(COUNT(*) FILTER (WHERE side = 'ref' AND val < 0) AS INTEGER) AS nn_ref,
          CAST(COUNT(*) FILTER (WHERE side = 'ref' AND val > 0) AS INTEGER) AS np_ref,
          list_sort(LIST(val) FILTER (WHERE side = 'tgt' AND val < 0)) AS neg_tgt,
          list_sort(LIST(val) FILTER (WHERE side = 'tgt' AND val > 0)) AS pos_tgt,
          list_sort(LIST(val) FILTER (WHERE side = 'ref' AND val < 0)) AS neg_ref,
          list_sort(LIST(val) FILTER (WHERE side = 'ref' AND val > 0)) AS pos_ref
        FROM nz
        GROUP BY grp_key, marker_1, marker_2
      ) agg
      ON g.grp_key = agg.grp_key AND po.marker_1 = agg.marker_1 AND po.marker_2 = agg.marker_2
    ),
    pair_medians AS (
      SELECT
        grp_key, marker_1, marker_2, first_rn, n_tgt, n_ref, k_tgt, k_ref,
        neg_tgt, pos_tgt, neg_ref, pos_ref,
        0.5 * (
          CASE WHEN ((n_tgt + 1) // 2) <= nn_tgt THEN list_extract(neg_tgt, (n_tgt + 1) // 2)
               WHEN ((n_tgt + 1) // 2) <= (n_tgt - np_tgt) THEN 0.0
               ELSE list_extract(pos_tgt, ((n_tgt + 1) // 2) - (n_tgt - np_tgt)) END +
          CASE WHEN ((n_tgt + 2) // 2) <= nn_tgt THEN list_extract(neg_tgt, (n_tgt + 2) // 2)
               WHEN ((n_tgt + 2) // 2) <= (n_tgt - np_tgt) THEN 0.0
               ELSE list_extract(pos_tgt, ((n_tgt + 2) // 2) - (n_tgt - np_tgt)) END
        ) AS med_tgt,
        0.5 * (
          CASE WHEN ((n_ref + 1) // 2) <= nn_ref THEN list_extract(neg_ref, (n_ref + 1) // 2)
               WHEN ((n_ref + 1) // 2) <= (n_ref - np_ref) THEN 0.0
               ELSE list_extract(pos_ref, ((n_ref + 1) // 2) - (n_ref - np_ref)) END +
          CASE WHEN ((n_ref + 2) // 2) <= nn_ref THEN list_extract(neg_ref, (n_ref + 2) // 2)
               WHEN ((n_ref + 2) // 2) <= (n_ref - np_ref) THEN 0.0
               ELSE list_extract(pos_ref, ((n_ref + 2) // 2) - (n_ref - np_ref)) END
        ) AS med_ref
      FROM pair_counts
    )
    SELECT * FROM pair_medians ORDER BY grp_key, first_rn
    "
  )

  agg_df <- DBI::dbGetQuery(con, sql_query)
  all_grp_keys <- as.character(group_sizes_tbl$grp_key)
  empty_reasons <- set_names(rep(NA_character_, length(all_grp_keys)), all_grp_keys)

  if (nrow(agg_df) == 0L) {
    empty_reasons[] <- "min_pct"
  } else {
    agg_df$pct_tgt <- round(agg_df$k_tgt / agg_df$n_tgt, digits = 3)
    agg_df$pct_ref <- round(agg_df$k_ref / agg_df$n_ref, digits = 3)
    agg_df$diff_median <- agg_df$med_tgt - agg_df$med_ref

    pass_pct <- pmax(agg_df$pct_tgt, agg_df$pct_ref) >= min_pct
    pass_diff_pct <- pass_pct &
      ((pmax(agg_df$pct_tgt, agg_df$pct_ref) - pmin(agg_df$pct_tgt, agg_df$pct_ref)) >= min_diff_pct)
    pass_diff <- if (isTRUE(only_pos)) {
      pass_diff_pct & (agg_df$diff_median > 0) & (agg_df$diff_median >= diff_threshold)
    } else {
      pass_diff_pct & (abs(agg_df$diff_median) >= diff_threshold)
    }

    for (gk in all_grp_keys) {
      idx <- which(agg_df$grp_key == gk)
      if (!any(pass_pct[idx])) {
        empty_reasons[[gk]] <- "min_pct"
      } else if (!any(pass_diff_pct[idx])) {
        empty_reasons[[gk]] <- "min_diff_pct"
      } else if (!any(pass_diff[idx])) {
        empty_reasons[[gk]] <- if (isTRUE(only_pos)) "only_pos" else "min_diff"
      }
    }
    agg_df <- agg_df[pass_diff, , drop = FALSE]
  }

  if (nrow(agg_df) == 0L) {
    out_empty <- tibble(
      grp_key = character(),
      marker_1 = character(),
      marker_2 = character(),
      pct_tgt = numeric(),
      pct_ref = numeric(),
      med_tgt = numeric(),
      med_ref = numeric(),
      diff_median = numeric(),
      p = numeric()
    )
    attr(out_empty, "empty_reasons") <- empty_reasons
    return(out_empty)
  }

  agg_df$p <- .wilcox_from_sorted_lists(
    neg_tgt = agg_df$neg_tgt,
    pos_tgt = agg_df$pos_tgt,
    neg_ref = agg_df$neg_ref,
    pos_ref = agg_df$pos_ref,
    n_tgt = agg_df$n_tgt,
    n_ref = agg_df$n_ref
  )
  attr(agg_df, "empty_reasons") <- empty_reasons
  agg_df
}


#' Emit the same empty-comparison warning as `.wilcox_de_test()`.
#'
#' @param reason Character string indicating which filter dropped all features
#'   (\code{"min_pct"}, \code{"min_diff_pct"}, \code{"only_pos"}, or \code{"min_diff"}).
#' @param only_pos Logical indicating whether only positive median differences were tested.
#'
#' @return Nothing. Called for its warning side effect.
#'
#' @noRd
.warn_empty_wilcox_comparison <- function(reason, only_pos = FALSE) {
  if (identical(reason, "min_pct")) {
    cli::cli_warn("No features pass min_pct threshold; returning empty data.frame")
  } else if (identical(reason, "min_diff_pct")) {
    cli::cli_warn("No features pass min_diff_pct threshold; returning empty data.frame")
  } else if (identical(reason, "only_pos") || isTRUE(only_pos)) {
    cli::cli_warn(
      "No features have a higher median in the target group; returning empty data.frame"
    )
  } else {
    cli::cli_warn("No features pass min_diff threshold; returning empty data.frame")
  }
}


#' Fast single-pass DuckDB path for \code{DifferentialProximityAnalysis.Seurat(..., lazy = TRUE, method = "seurat")}.
#'
#' @param object A \code{Seurat} object.
#' @param group_data_tbl A tibble of cell metadata with a \code{component} column.
#' @param contrast_column The contrast column name.
#' @param reference The reference level.
#' @param targets Target levels of the contrast column.
#' @param assay Name of the assay to read.
#' @param group_vars Optional grouping variables.
#' @param min_exp_join_count Minimum expected join count for a marker pair.
#' @param min_cells_per_group Minimum cells required on each side of a comparison.
#' @param diff_threshold Minimum absolute median difference required to test a pair.
#' @param min_pct Minimum fraction of cells with a non-zero score in either group.
#' @param min_diff_pct Minimum difference in that fraction between the two groups.
#' @param proximity_metric Name of the numeric proximity column to test.
#' @param metric_type One of \code{"all"}, \code{"self"}, or \code{"co"}.
#' @param p_adjust_method Method for multiple-testing adjustment.
#' @param verbose Logical indicating whether to print setup and progress messages.
#' @param only_pos If \code{TRUE}, test only pairs with a higher median in the target.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{tbl_df} of differential proximity results.
#'
#' @noRd
.run_lazy_duckdb_dpa_seurat <- function(
  object,
  group_data_tbl,
  contrast_column,
  reference,
  targets,
  assay,
  group_vars,
  min_exp_join_count,
  min_cells_per_group,
  diff_threshold,
  min_pct,
  min_diff_pct,
  proximity_metric,
  metric_type,
  p_adjust_method,
  verbose,
  only_pos = FALSE,
  ...
) {
  call <- caller_env()
  assert_single_value(only_pos, type = "bool", call = call)
  metric_type <- match.arg(metric_type, choices = c("all", "self", "co"))
  p_adjust_method <- match.arg(
    p_adjust_method,
    choices = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr")
  )

  group_data_df <- column_to_rownames(group_data_tbl, "component")
  dummy_mat <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    dims = c(1L, ncol(object)),
    dimnames = list("pair1", colnames(object))
  )
  .validate_matrix_group_data(
    object = dummy_mat,
    group_data = group_data_df,
    contrast_column = contrast_column,
    reference = reference,
    targets = targets,
    group_vars = group_vars,
    call = call
  )
  assert_single_value(min_cells_per_group, type = "numeric", call = call)
  if (is.na(min_cells_per_group) || min_cells_per_group <= 0) {
    cli::cli_abort(
      c("x" = "{.arg min_cells_per_group} must be a positive number (> 0)."),
      call = call
    )
  }

  targets <- targets %||% setdiff(unique(group_data_df[, contrast_column, drop = TRUE]), reference)
  grouped <- .compute_test_groups_keys(
    data = group_data_tbl,
    group_cols = c(contrast_column, group_vars),
    contrast_column = contrast_column,
    targets = targets
  )
  test_groups <- grouped$test_groups
  test_groups_keys <- grouped$keys

  proximity_data <- ProximityScores(object, assay = assay, lazy = TRUE)
  con <- proximity_data$src$con
  on.exit(
    {
      if (!is.null(con) && DBI::dbIsValid(con)) {
        DBI::dbDisconnect(con, shutdown = TRUE)
      }
    },
    add = TRUE
  )
  assert_single_value(proximity_metric, type = "string", call = call)
  assert_col_in_data(proximity_metric, proximity_data, call = call)
  assert_col_class(proximity_metric, proximity_data, classes = "numeric", call = call)

  proximity_data <- switch(metric_type,
    all = proximity_data,
    self = proximity_data |> filter(marker_1 == marker_2),
    co = proximity_data |> filter(marker_1 != marker_2)
  ) |>
    filter(join_count_expected_mean >= min_exp_join_count)

  .print_da_setup(
    keys = test_groups_keys,
    contrast_column = contrast_column,
    reference = reference,
    targets = targets,
    group_vars = group_vars,
    verbose = verbose,
    test_groups = test_groups,
    min_cells_per_group = min_cells_per_group
  )

  keys <- apply(test_groups_keys, 1, paste, collapse = "_")
  diff_prox_res <- set_names(rep(list(NULL), length(keys)), keys)

  cell_rows <- list()
  group_size_rows <- list()
  valid_key_indices <- integer()

  for (i in seq_len(nrow(test_groups_keys))) {
    cur_keys <- test_groups_keys[i, ]
    grp_id <- paste0("g_", sprintf("%06d", i))
    cur_test_groups <- .filter_cur_group(
      test_groups = test_groups,
      cur_keys = cur_keys,
      contrast_column = contrast_column,
      reference = reference,
      group_vars = group_vars
    )
    components_reference <- cur_test_groups |>
      filter(!!sym(contrast_column) == reference) |>
      pull(component)
    components_target <- cur_test_groups |>
      filter(!!sym(contrast_column) != reference) |>
      pull(component)

    n_comps_reference <- length(components_reference)
    n_comps_target <- length(components_target)
    if (n_comps_reference < min_cells_per_group || n_comps_target < min_cells_per_group) {
      target <- cur_keys[, contrast_column, drop = TRUE]
      if (ncol(cur_keys) > 1) {
        cur_reference <- paste0(reference, " (", paste(cur_keys[, group_vars, drop = TRUE], collapse = ", "), ")")
        cur_target <- paste0(target, " (", paste(cur_keys[, group_vars, drop = TRUE], collapse = ", "), ")")
      } else {
        cur_reference <- reference
        cur_target <- target
      }
      if (n_comps_reference < min_cells_per_group && n_comps_target < min_cells_per_group) {
        cli::cli_warn(
          "Skipping {cur_target} vs {cur_reference} because both
           groups have fewer than {.val {min_cells_per_group}} cells."
        )
      } else if (n_comps_reference < min_cells_per_group && n_comps_target >= min_cells_per_group) {
        cli::cli_warn(
          "Skipping {cur_target} vs {cur_reference} because the reference
           group has fewer than {.val {min_cells_per_group}} cells."
        )
      } else if (n_comps_target < min_cells_per_group && n_comps_reference >= min_cells_per_group) {
        cli::cli_warn(
          "Skipping {cur_target} vs {cur_reference} because the target group
           has fewer than {.val {min_cells_per_group}} cells."
        )
      }
      next
    }

    valid_key_indices <- c(valid_key_indices, i)
    cell_rows[[length(cell_rows) + 1L]] <- tibble(
      component = c(components_target, components_reference),
      grp_key = rep(grp_id, n_comps_target + n_comps_reference),
      side = c(rep("tgt", n_comps_target), rep("ref", n_comps_reference))
    )
    group_size_rows[[length(group_size_rows) + 1L]] <- tibble(
      grp_key = grp_id,
      n_tgt = as.integer(n_comps_target),
      n_ref = as.integer(n_comps_reference)
    )
  }

  if (length(valid_key_indices) > 0L) {
    agg_df <- .run_duckdb_wilcox_comparisons(
      lazy_prox = proximity_data,
      cells_tbl = bind_rows(cell_rows),
      group_sizes_tbl = bind_rows(group_size_rows),
      proximity_metric = proximity_metric,
      diff_threshold = diff_threshold,
      min_pct = min_pct,
      min_diff_pct = min_diff_pct,
      only_pos = only_pos,
      call = call
    )

    empty_reasons <- attr(agg_df, "empty_reasons")
    for (i in valid_key_indices) {
      cur_keys <- test_groups_keys[i, ]
      key <- paste(as.character(cur_keys), collapse = "_")
      grp_id <- paste0("g_", sprintf("%06d", i))
      if (verbose && check_global_verbosity()) {
        sz_row <- group_size_rows[[match(i, valid_key_indices)]]
        cli_progress_update(
          set = i,
          status = .format_da_comparison_status(
            cur_keys = cur_keys,
            contrast_column = contrast_column,
            reference = reference,
            group_vars = group_vars,
            n_tgt = sz_row$n_tgt,
            n_ref = sz_row$n_ref
          )
        )
      }
      sub_df <- agg_df[agg_df$grp_key == grp_id, , drop = FALSE]
      if (nrow(sub_df) == 0L) {
        .warn_empty_wilcox_comparison(
          reason = empty_reasons[[grp_id]],
          only_pos = only_pos
        )
        next
      }
      diff_prox_res[[key]] <- tibble(
        target = cur_keys[, contrast_column, drop = TRUE],
        reference = reference,
        pct_tgt = sub_df$pct_tgt,
        pct_ref = sub_df$pct_ref,
        med_tgt = sub_df$med_tgt,
        med_ref = sub_df$med_ref,
        diff_median = sub_df$diff_median,
        p = sub_df$p,
        marker_1 = sub_df$marker_1,
        marker_2 = sub_df$marker_2
      ) |>
        .append_group_vars(cur_keys = cur_keys, group_vars = group_vars)
    }
  }

  if (verbose && check_global_verbosity()) {
    cli_progress_done()
  }

  valid_results <- sapply(diff_prox_res, function(x) !is.null(x))
  if (!any(valid_results)) {
    .abort_no_valid_da_results(
      n_total = nrow(test_groups_keys),
      n_passed_cells = length(valid_key_indices),
      min_cells_per_group = min_cells_per_group,
      contrast_column = contrast_column,
      group_vars = group_vars,
      diff_threshold = diff_threshold,
      min_pct = min_pct,
      min_diff_pct = min_diff_pct,
      only_pos = only_pos,
      call = call
    )
  }

  .finalize_da_results(diff_prox_res, p_adjust_method, group_vars)
}


#' Fast single-pass DuckDB path for \code{FindAllProximityMarkers.Seurat(..., lazy = TRUE)}.
#'
#' @param object A \code{Seurat} object.
#' @param assay Name of the assay to read.
#' @param label_chr Named character vector of group labels keyed by cell ID.
#' @param idents Character vector of group levels to test.
#' @param rest_label Label used for the pooled remainder group.
#' @param group_by Metadata column name defining the groups.
#' @param min_exp_join_count Minimum expected join count for a marker pair.
#' @param min_cells_per_group Minimum cells required on each side of a comparison.
#' @param diff_threshold Minimum absolute median difference required to test a pair.
#' @param min_pct Minimum fraction of cells with a non-zero score in either group.
#' @param min_diff_pct Minimum difference in that fraction between the two groups.
#' @param only_pos If \code{TRUE}, test only pairs with a higher median in the target.
#' @param proximity_metric Name of the numeric proximity column to test.
#' @param metric_type One of \code{"all"}, \code{"self"}, or \code{"co"}.
#' @param p_adjust_method Method for multiple-testing adjustment.
#' @param verbose Logical indicating whether to show a progress bar.
#' @param call Calling environment used in error messages.
#'
#' @return A \code{tbl_df} of one-versus-rest differential proximity results.
#'
#' @noRd
.run_lazy_duckdb_one_vs_rest_seurat <- function(
  object,
  assay,
  label_chr,
  idents,
  rest_label,
  group_by,
  min_exp_join_count,
  min_cells_per_group,
  diff_threshold,
  min_pct,
  min_diff_pct,
  only_pos,
  proximity_metric,
  metric_type,
  p_adjust_method,
  verbose,
  call = caller_env()
) {
  proximity_data <- ProximityScores(object, assay = assay, lazy = TRUE)
  con <- proximity_data$src$con
  on.exit(
    {
      if (!is.null(con) && DBI::dbIsValid(con)) {
        DBI::dbDisconnect(con, shutdown = TRUE)
      }
    },
    add = TRUE
  )
  assert_single_value(proximity_metric, type = "string", call = call)
  assert_col_in_data(proximity_metric, proximity_data, call = call)
  assert_col_class(proximity_metric, proximity_data, classes = "numeric", call = call)

  proximity_data <- switch(metric_type,
    all = proximity_data,
    self = proximity_data |> filter(marker_1 == marker_2),
    co = proximity_data |> filter(marker_1 != marker_2)
  ) |>
    filter(join_count_expected_mean >= min_exp_join_count)

  cells <- colnames(object)
  cell_labels <- .labels_for_cells(label_chr, cells, group_by)

  if (verbose && check_global_verbosity()) {
    cli::cli_progress_bar(
      format = paste0(
        "Running one-versus-rest proximity tests ",
        "[{cli::pb_bar}] {cli::pb_percent} | ",
        "{cli::pb_current}/{cli::pb_total} | {cli::pb_status}"
      ),
      total = length(idents),
      clear = FALSE,
      status = paste0("starting (rest = ", rest_label, ")")
    )
  }

  cell_rows <- list()
  group_size_rows <- list()
  valid_indices <- integer()

  for (i in seq_along(idents)) {
    ident <- idents[[i]]
    contrast <- .assign_one_vs_rest(cell_labels, ident, rest_label)
    cells_target <- cells[!is.na(contrast) & contrast == ident]
    cells_reference <- cells[!is.na(contrast) & contrast == rest_label]
    n_tgt <- length(cells_target)
    n_ref <- length(cells_reference)
    if (n_tgt < min_cells_per_group || n_ref < min_cells_per_group) {
      .warn_skip_one_vs_rest(
        ident = ident,
        rest_label = rest_label,
        n_tgt = n_tgt,
        n_ref = n_ref,
        min_cells_per_group = min_cells_per_group
      )
      next
    }
    valid_indices <- c(valid_indices, i)
    grp_id <- paste0("g_", sprintf("%06d", i))
    cell_rows[[length(cell_rows) + 1L]] <- tibble(
      component = c(cells_target, cells_reference),
      grp_key = rep(grp_id, n_tgt + n_ref),
      side = c(rep("tgt", n_tgt), rep("ref", n_ref))
    )
    group_size_rows[[length(group_size_rows) + 1L]] <- tibble(
      grp_key = grp_id,
      n_tgt = as.integer(n_tgt),
      n_ref = as.integer(n_ref)
    )
  }

  results <- vector("list", length(idents))
  if (length(valid_indices) > 0L) {
    agg_df <- .run_duckdb_wilcox_comparisons(
      lazy_prox = proximity_data,
      cells_tbl = bind_rows(cell_rows),
      group_sizes_tbl = bind_rows(group_size_rows),
      proximity_metric = proximity_metric,
      diff_threshold = diff_threshold,
      min_pct = min_pct,
      min_diff_pct = min_diff_pct,
      only_pos = only_pos,
      call = call
    )

    empty_reasons <- attr(agg_df, "empty_reasons")
    for (i in valid_indices) {
      ident <- idents[[i]]
      grp_id <- paste0("g_", sprintf("%06d", i))
      if (verbose && check_global_verbosity()) {
        cli::cli_progress_update(set = i)
      }
      sub_df <- agg_df[agg_df$grp_key == grp_id, , drop = FALSE]
      if (nrow(sub_df) == 0L) {
        .warn_empty_wilcox_comparison(
          reason = empty_reasons[[grp_id]],
          only_pos = only_pos
        )
        next
      }
      results[[i]] <- tibble(
        target = ident,
        pct_tgt = sub_df$pct_tgt,
        pct_ref = sub_df$pct_ref,
        med_tgt = sub_df$med_tgt,
        med_rest = sub_df$med_ref,
        diff_median = sub_df$diff_median,
        p = sub_df$p,
        p_adj = p.adjust(sub_df$p, method = p_adjust_method),
        marker_1 = sub_df$marker_1,
        marker_2 = sub_df$marker_2
      )
    }
  }

  out <- bind_rows(results)
  if (nrow(out) == 0L) {
    cli::cli_abort(
      c(
        "x" = "No one-versus-rest proximity tests produced results.",
        "i" = "Check {.arg group_by} and {.arg min_cells_per_group}."
      ),
      call = call
    )
  }
  if (verbose && check_global_verbosity()) {
    cli::cli_progress_done()
  }
  out
}


#' Add zero columns for cells that have no proximity scores.
#'
#' @param proximity_matrix Sparse matrix of marker pairs by cells.
#' @param cells Cell ids, in the column order of the returned matrix.
#'
#' @return A sparse matrix with one column per cell in \code{cells}.
#'
#' @noRd
.pad_proximity_matrix <- function(proximity_matrix, cells) {
  cells <- as.character(cells)
  missing_components <- setdiff(cells, colnames(proximity_matrix))
  if (length(missing_components) > 0) {
    m_missing <- Matrix::rsparsematrix(
      nrow = nrow(proximity_matrix),
      ncol = length(missing_components),
      density = 0
    )
    rownames(m_missing) <- rownames(proximity_matrix)
    colnames(m_missing) <- missing_components
    proximity_matrix <- cbind(proximity_matrix, m_missing)
  }
  proximity_matrix[, cells, drop = FALSE]
}


#' Stop when the caller passes unsupported arguments to `FindAllProximityMarkers()`.
#'
#' @param ... Arguments supplied by the caller.
#' @param call Calling environment used in the error message.
#'
#' @return Nothing. Signals an error when any extra argument is present.
#'
#' @noRd
.reject_find_all_proximity_args <- function(..., call = caller_env()) {
  supplied_names <- names(list(...))
  blocked <- c(
    "contrast_column", "reference", "targets", "group_vars", "group_data"
  )
  supplied <- intersect(supplied_names, blocked)
  if (length(supplied) > 0) {
    cli::cli_abort(
      c(
        "x" = "Cannot pass {.arg {supplied}} to {.fn FindAllProximityMarkers}.",
        "i" = paste0(
          "{.arg group_by} selects the groups. Each level is compared to ",
          "all other levels."
        )
      ),
      call = call
    )
  }
  legacy <- intersect(supplied_names, c("method", "backend"))
  if (length(legacy) > 0) {
    cli::cli_abort(
      c(
        "x" = "Cannot pass {.arg {legacy}} to {.fn FindAllProximityMarkers}.",
        "i" = paste0(
          "Each comparison uses the matrix Wilcoxon test on one proximity ",
          "matrix."
        )
      ),
      call = call
    )
  }
  rlang::check_dots_empty(call = call)
}


#' Check the arguments shared by the \code{FindAllProximityMarkers()} methods.
#'
#' @param min_cells_per_group Minimum cells required on each side of a comparison.
#' @param diff_threshold Minimum absolute median difference required to test a pair.
#' @param min_pct Minimum fraction of cells with a non-zero score in either group.
#' @param min_diff_pct Minimum difference in that fraction between the two groups.
#' @param only_pos If \code{TRUE}, test only pairs with a higher median in the target.
#' @param proximity_metric Name of the numeric proximity column to test.
#' @param verbose Print progress messages.
#' @param call Calling environment used in error messages.
#'
#' @return Nothing. Signals an error when an argument is invalid.
#'
#' @noRd
.validate_find_all_common_args <- function(
  min_cells_per_group,
  diff_threshold,
  min_pct,
  min_diff_pct,
  only_pos,
  proximity_metric,
  verbose,
  call = caller_env()
) {
  assert_single_value(min_cells_per_group, type = "numeric", call = call)
  if (is.na(min_cells_per_group) || min_cells_per_group <= 0) {
    cli::cli_abort(
      c("x" = "{.arg min_cells_per_group} must be a positive number (> 0)."),
      call = call
    )
  }
  assert_single_value(diff_threshold, type = "numeric", call = call)
  if (is.na(diff_threshold) || diff_threshold < 0) {
    cli::cli_abort(
      c("x" = "{.arg diff_threshold} must be a non-negative number."),
      call = call
    )
  }
  assert_single_value(min_pct, type = "numeric", call = call)
  if (is.na(min_pct) || min_pct < 0 || min_pct > 1) {
    cli::cli_abort(
      c("x" = "{.arg min_pct} must be a number between 0 and 1."),
      call = call
    )
  }
  assert_single_value(min_diff_pct, type = "numeric", call = call)
  if (is.na(min_diff_pct)) {
    cli::cli_abort(
      c("x" = "{.arg min_diff_pct} must be a number."),
      call = call
    )
  }
  assert_single_value(only_pos, type = "bool", call = call)
  assert_single_value(proximity_metric, type = "string", call = call)
  assert_single_value(verbose, type = "bool", call = call)
}


#' Coerce requested group levels to character.
#'
#' @param idents Levels to test. Factors and numbers are coerced.
#'
#' @return A character vector, or \code{idents} unchanged when it is already character.
#'
#' @noRd
.coerce_idents <- function(idents) {
  if (is.factor(idents) || is.numeric(idents)) {
    return(as.character(idents))
  }
  idents
}


#' Choose which group levels to test against the pooled remainder.
#'
#' @param labels Group label for each cell.
#' @param idents Levels to test, or \code{NULL} to test every level.
#' @param group_by Name of the grouping column, used in error messages.
#' @param call Calling environment used in error messages.
#'
#' @return A character vector of levels to test.
#'
#' @noRd
.resolve_one_vs_rest_idents <- function(
  labels,
  idents,
  group_by,
  call = caller_env()
) {
  if (is.factor(labels)) {
    all_idents <- levels(droplevels(labels))
  } else {
    all_idents <- unique(as.character(labels))
  }
  all_idents <- all_idents[!is.na(all_idents)]
  if (length(all_idents) < 2) {
    n_groups <- length(all_idents)
    cli::cli_abort(
      c(
        "i" = "Group variable {.val {group_by}} must have at least 2 groups.",
        "x" = "Group variable {.val {group_by}} has {n_groups} unique group{?s}."
      ),
      call = call
    )
  }

  if (is.null(idents)) {
    return(all_idents)
  }

  idents <- .coerce_idents(idents)
  assert_vector(idents, type = "character", n = 1, arg = "idents", call = call)
  if (anyDuplicated(idents)) {
    cli::cli_abort("{.arg idents} must be unique.", call = call)
  }
  missing <- setdiff(idents, all_idents)
  if (length(missing) > 0) {
    cli::cli_abort(
      "Not all {.arg idents} were found in {.val {group_by}}: {.val {missing}}.",
      call = call
    )
  }
  idents
}


#' Choose a label for the pooled remainder that is not already a group level.
#'
#' @param labels Group labels already in use.
#'
#' @return \code{"rest"}, with extra underscores appended until the label is unused.
#'
#' @noRd
.make_rest_label <- function(labels) {
  used <- unique(as.character(labels))
  used <- used[!is.na(used)]
  candidate <- "rest"
  while (candidate %in% used) {
    candidate <- paste0(candidate, "_")
  }
  candidate
}


#' Label one level as the target and every other known level as the remainder.
#'
#' @param labels Group label for each cell.
#' @param ident Level to test.
#' @param rest_label Label for the pooled remainder.
#'
#' @return A character vector. Missing labels stay missing.
#'
#' @noRd
.assign_one_vs_rest <- function(labels, ident, rest_label) {
  labels <- as.character(labels)
  contrast <- rep(NA_character_, length(labels))
  known <- !is.na(labels)
  contrast[known & labels == ident] <- ident
  contrast[known & labels != ident] <- rest_label
  contrast
}


#' Warn when some cells have no group label.
#'
#' @param label_chr Group label for each cell.
#' @param group_by Name of the grouping column.
#'
#' @return Nothing. Called for its warning.
#'
#' @noRd
.warn_missing_group_labels <- function(label_chr, group_by) {
  n_missing <- sum(is.na(label_chr))
  if (n_missing > 0) {
    cli::cli_warn(
      "Excluding {n_missing} cell{?s} with missing {.val {group_by}} labels."
    )
  }
}


#' Align named cell labels to a set of cells.
#'
#' @param label_chr Group labels named by cell id.
#' @param cells Cell ids to keep, in the desired order.
#' @param group_by Name of the grouping column, used in the warning.
#'
#' @return An unnamed character vector, with \code{NA} where a cell has no label.
#'
#' @noRd
.labels_for_cells <- function(label_chr, cells, group_by) {
  cells <- as.character(cells)
  missing <- setdiff(cells, names(label_chr))
  if (length(missing) > 0) {
    n_missing <- length(missing)
    cli::cli_warn(
      "Excluding {n_missing} component{?s} with no {.val {group_by}} label."
    )
  }
  unname(label_chr[cells])
}


#' Run one-versus-rest Wilcoxon tests on one proximity matrix.
#'
#' @param proximity_matrix Sparse matrix of marker pairs by cells.
#' @param label_chr Group labels named by cell id.
#' @param idents Levels to test.
#' @param rest_label Label for the pooled remainder.
#' @param group_by Name of the grouping column, used in messages.
#' @param min_cells_per_group Minimum cells required on each side of a comparison.
#' @param diff_threshold Minimum absolute median difference required to test a pair.
#' @param min_pct Minimum fraction of cells with a non-zero score in either group.
#' @param min_diff_pct Minimum difference in that fraction between the two groups.
#' @param only_pos If \code{TRUE}, test only pairs with a higher median in the target.
#' @param proximity_metric Name of the proximity metric being tested.
#' @param p_adjust_method Method passed to \code{p.adjust()} within each comparison.
#' @param verbose Print progress messages.
#' @param call Calling environment used in the error when no comparison succeeds.
#'
#' @return A tibble of test results stacked across comparisons.
#'
#' @noRd
.run_matrix_one_vs_rest <- function(
  proximity_matrix,
  label_chr,
  idents,
  rest_label,
  group_by,
  min_cells_per_group,
  diff_threshold,
  min_pct,
  min_diff_pct,
  only_pos,
  proximity_metric,
  p_adjust_method,
  verbose,
  call = caller_env()
) {
  cells <- colnames(proximity_matrix)
  cell_labels <- .labels_for_cells(label_chr, cells, group_by)

  if (verbose && check_global_verbosity()) {
    cli::cli_progress_bar(
      format = paste0(
        "Running one-versus-rest proximity tests ",
        "[{cli::pb_bar}] {cli::pb_percent} | ",
        "{cli::pb_current}/{cli::pb_total} | {cli::pb_status}"
      ),
      total = length(idents),
      clear = FALSE,
      status = paste0("starting (rest = ", rest_label, ")")
    )
  }

  results <- vector("list", length(idents))
  for (i in seq_along(idents)) {
    ident <- idents[[i]]
    contrast <- .assign_one_vs_rest(cell_labels, ident, rest_label)
    is_target <- !is.na(contrast) & contrast == ident
    is_reference <- !is.na(contrast) & contrast == rest_label
    cells_target <- cells[is_target]
    cells_reference <- cells[is_reference]
    n_tgt <- length(cells_target)
    n_ref <- length(cells_reference)
    if (n_tgt < min_cells_per_group || n_ref < min_cells_per_group) {
      .warn_skip_one_vs_rest(
        ident = ident,
        rest_label = rest_label,
        n_tgt = n_tgt,
        n_ref = n_ref,
        min_cells_per_group = min_cells_per_group
      )
      next
    }

    if (verbose && check_global_verbosity()) {
      cli::cli_progress_update(
        set = i,
        status = paste0(ident, " vs ", rest_label, " (", n_tgt, " vs ", n_ref, " cells)")
      )
    }

    de_results <- .wilcox_de_test(
      proximity_matrix,
      cells_1 = cells_target,
      cells_2 = cells_reference,
      min_diff = diff_threshold,
      min_pct = min_pct,
      min_diff_pct = min_diff_pct,
      only_pos = only_pos
    )

    if (nrow(de_results) == 0) {
      next
    }

    results[[i]] <- tibble(
      target = ident,
      pct_tgt = de_results$pct_1,
      pct_ref = de_results$pct_2,
      med_tgt = de_results$med_1,
      med_rest = de_results$med_2,
      diff_median = de_results$difference,
      p = de_results$p_val,
      p_adj = p.adjust(de_results$p_val, method = p_adjust_method),
      marker_1 = sub(":.*$", "", rownames(de_results)),
      marker_2 = sub("^[^:]*:", "", rownames(de_results))
    )
  }

  out <- bind_rows(results)
  if (nrow(out) == 0) {
    cli::cli_abort(
      c(
        "x" = "No one-versus-rest proximity tests produced results.",
        "i" = "Check {.arg group_by} and {.arg min_cells_per_group}."
      ),
      call = call
    )
  }
  if (verbose && check_global_verbosity()) {
    cli::cli_progress_done()
  }
  out
}


#' Warn that a one-versus-rest comparison has too few cells.
#'
#' @param ident Level being tested.
#' @param rest_label Label for the pooled remainder.
#' @param n_tgt Number of cells in \code{ident}.
#' @param n_ref Number of cells in the pooled remainder.
#' @param min_cells_per_group Minimum cells required on each side.
#'
#' @return Nothing. Called for its warning.
#'
#' @noRd
.warn_skip_one_vs_rest <- function(
  ident,
  rest_label,
  n_tgt,
  n_ref,
  min_cells_per_group
) {
  if (n_tgt < min_cells_per_group && n_ref < min_cells_per_group) {
    reason <- "both groups have"
  } else if (n_tgt < min_cells_per_group) {
    reason <- "the target group has"
  } else {
    reason <- "the reference group has"
  }
  cli::cli_warn(
    paste0(
      "Skipping {.val {ident}} vs {.val {rest_label}} because {reason} ",
      "fewer than {.val {min_cells_per_group}} cells."
    )
  )
}


#' Add a proximity transition column to differential proximity results
#'
#' Classifies the reference/rest and target group medians into proximity regimes
#' around zero (\code{"anti"} for \code{< -tau}, \code{"random"} for
#' \code{[-tau, tau]}, and \code{"coloc"} for \code{> tau}) and combines them
#' into a \code{transition} column formatted as \code{"<ref_state> -> <tgt_state>"}.
#'
#' Works on output tables from both \code{\link{DifferentialProximityAnalysis}}
#' (using \code{med_ref} and \code{med_tgt}) and
#' \code{\link{FindAllProximityMarkers}} (using \code{med_rest} and
#' \code{med_tgt}).
#'
#' @param data A \code{data.frame} or \code{tbl_df} returned by
#'   \code{\link{DifferentialProximityAnalysis}} or
#'   \code{\link{FindAllProximityMarkers}}.
#' @param tau Non-negative numeric threshold defining the neutral
#'   \code{"random"} zone \code{[-tau, tau]} around 0. Default is \code{0.15}.
#' @param tgt_col Name of the target median column. Default is \code{"med_tgt"}.
#' @param ref_col Name of the reference median column. When \code{NULL}
#'   (default), automatically uses \code{"med_ref"} or \code{"med_rest"} if
#'   present in \code{data}.
#' @param transition_col Name of the output transition column. Default is
#'   \code{"transition"}.
#' @param separator Separator placed between the reference and target states.
#'   Default is \code{" -> "}.
#'
#' @return \code{data} with a character column named by
#'   \code{transition_col}. An existing column of that name is replaced.
#'
#' @examples
#' library(dplyr)
#' df <- tibble(
#'   target = c("B", "CD4 T", "NK"),
#'   med_tgt = c(0.8, 0.9, -0.4),
#'   med_rest = c(0.0, 0.5, 0.0),
#'   diff_median = med_tgt - med_rest
#' )
#' add_proximity_transition(df, tau = 0.15)
#'
#' @export
#'
add_proximity_transition <- function(
  data,
  tau = 0.15,
  tgt_col = "med_tgt",
  ref_col = NULL,
  transition_col = "transition",
  separator = " -> "
) {
  assert_class(data, "data.frame")
  assert_single_value(tau, type = "numeric")
  if (is.na(tau) || tau < 0) {
    cli::cli_abort(
      c("x" = "{.arg tau} must be a non-negative number.")
    )
  }
  assert_single_value(tgt_col, type = "string")
  assert_col_in_data(tgt_col, data)
  assert_col_class(tgt_col, data, classes = "numeric")
  assert_single_value(transition_col, type = "string")
  assert_single_value(separator, type = "string")

  if (is.null(ref_col)) {
    candidates <- intersect(c("med_ref", "med_rest"), colnames(data))
    if (length(candidates) == 0) {
      cli::cli_abort(
        c(
          "x" = "Could not find a reference median column in {.arg data}.",
          "i" = "Expected {.val med_ref} or {.val med_rest}, or supply {.arg ref_col}."
        )
      )
    }
    ref_col <- candidates[[1]]
  } else {
    assert_single_value(ref_col, type = "string")
    assert_col_in_data(ref_col, data)
  }
  assert_col_class(ref_col, data, classes = "numeric")

  tgt_vals <- data[[tgt_col]]
  ref_vals <- data[[ref_col]]

  state_tgt <- case_when(
    is.na(tgt_vals) ~ NA_character_,
    tgt_vals > tau ~ "coloc",
    tgt_vals < -tau ~ "anti",
    TRUE ~ "random"
  )
  state_ref <- case_when(
    is.na(ref_vals) ~ NA_character_,
    ref_vals > tau ~ "coloc",
    ref_vals < -tau ~ "anti",
    TRUE ~ "random"
  )

  data[[transition_col]] <- ifelse(
    is.na(state_ref) | is.na(state_tgt),
    NA_character_,
    paste0(state_ref, separator, state_tgt)
  )
  data
}
