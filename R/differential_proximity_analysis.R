#' @include generics.R
NULL

.datatable.aware <- TRUE


#' @rdname DifferentialProximityAnalysis
#' @method DifferentialProximityAnalysis data.frame
#'
#' @examples
#' library(dplyr)
#' example_data <- tidyr::expand_grid(
#'   marker_1 = c("HLA-ABC", "B2M", "CD4", "CD8", "CD20", "CD19", "CD45", "CD43") %>%
#'     rep(each = 50),
#'   marker_2 = c("HLA-ABC", "B2M", "CD4", "CD8", "CD20", "CD19", "CD45", "CD43") %>%
#'     rep(each = 50)
#' ) %>%
#'   mutate(
#'     join_count_z = rnorm(n(), sd = 10)
#'   )
#'
#' example_data <- example_data %>%
#'   mutate(sampleID = "ctrl") %>%
#'   bind_rows(
#'     example_data %>% mutate(join_count_z = join_count_z + 1) %>%
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
  backend <- match.arg(backend, choices = c("dplyr", "data.table"))
  if (backend == "data.table") {
    expect_dtplyr()
  }
  metric_type <- match.arg(metric_type, choices = c("all", "self", "co"))
  object <- switch(metric_type,
    all = object,
    self = object %>% filter(marker_1 == marker_2),
    co = object %>% filter(marker_1 != marker_2)
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
  object <- object %>%
    select(any_of(c("marker_1", "marker_2", proximity_metric, "component", contrast_column, group_vars)))

  # Group data by marker pair, contrast column and optional group variables
  test_groups <- object %>%
    group_by(pick(all_of(c("marker_1", "marker_2", contrast_column, group_vars))))

  if (min_cells_per_group > 0) {
    if (backend == "data.table") {
      test_groups <- dtplyr::lazy_dt(test_groups)
    }

    # Compress the massive table into a tiny count table
    valid_groups <- test_groups %>%
      count(name = "cell_count") %>%
      filter(cell_count >= min_cells_per_group) %>%
      ungroup(!!sym(contrast_column)) %>%
      mutate(n_groups = n(), ref_present = sum(!!sym(contrast_column) == reference)) %>%
      filter(n_groups > 1, ref_present == 1) %>%
      select(-n_groups, -ref_present)

    # Filter the massive original table using the valid combinations found above
    test_groups <- test_groups %>%
      ungroup() %>%
      semi_join(valid_groups, by = c("marker_1", "marker_2", contrast_column, group_vars)) %>%
      as_tibble() %>%
      group_by(pick(all_of(c("marker_1", "marker_2", contrast_column, group_vars))))
    if (nrow(test_groups) == 0) {
      cli::cli_abort(glue(
        "Found no groups with at least {min_cells_per_group} observations."
      ))
    }
  }

  # Compute group keys (drop marker columns; this method groups by marker pair
  # internally but iterates over contrast/group-var combinations only)
  test_groups_keys <- test_groups %>%
    group_keys() %>%
    select(-marker_1, -marker_2) %>%
    filter(!!sym(contrast_column) %in% targets) %>%
    distinct() %>%
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
    verbose = verbose
  )

  # Create a list to store the results in
  keys <- apply(test_groups_keys, 1, paste, collapse = "_")
  diff_prox_res <- rep(list(NULL), length(keys)) %>%
    set_names(keys)

  # iterate over group pairs
  for (i in seq_len(nrow(test_groups_keys))) {
    if (verbose && check_global_verbosity()) {
      cli_progress_update()
    }

    # Fetch the current group keys
    # If group_vars=NULL, the key is simply one of targets.
    # Otherwise, the key is target + additional group_vars.
    cur_keys <- test_groups_keys[i, ]
    key <- as.character(cur_keys) %>% paste(collapse = "_")
    groups_keep <- c(reference, cur_keys[, contrast_column, drop = TRUE])

    # Filter data to the current comparison (contrast + group_vars)
    cur_test_groups <- .filter_cur_group(
      test_groups = test_groups,
      cur_keys = cur_keys,
      contrast_column = contrast_column,
      reference = reference,
      group_vars = group_vars
    )

    # Group by marker pair and compute ranks
    cur_test_groups_rank <- cur_test_groups %>%
      {
        if (backend == "data.table") {
          dtplyr::lazy_dt(.) %>%
            group_by(marker_1, marker_2) %>%
            mutate(r = data.table::frank(-!!sym(proximity_metric), ties.method = "average"))
        } else {
          group_by(., marker_1, marker_2) %>%
            mutate(r = rank(-!!sym(proximity_metric), ties.method = "average"))
        }
      }

    # Compute the number of ties per marker pair
    cur_test_groups_ties <- cur_test_groups_rank %>%
      group_by(marker_1, marker_2, r) %>%
      summarize(nties = n(), .groups = "drop") %>%
      group_by(marker_1, marker_2) %>%
      summarize(nties_const = sum(nties^3 - nties), .groups = "drop") %>%
      collect()

    # Compute the U statistic, Z score and p-value
    cur_test_groups_u <- cur_test_groups_rank %>%
      group_by(!!sym(contrast_column), .add = TRUE) %>%
      summarize(rank_sum = sum(r), n = n(), med = median(!!sym(proximity_metric)), .groups = "drop") %>%
      pivot_wider(names_from = !!sym(contrast_column), values_from = c("rank_sum", "n", "med")) %>%
      select(all_of(c(
        "marker_1", "marker_2",
        paste0("rank_sum_", groups_keep),
        paste0("n_", groups_keep),
        paste0("med_", groups_keep)
      ))) %>%
      rename(rs_ref = 3L, rs_tgt = 4L, n_ref = 5L, n_tgt = 6L, med_ref = 7L, med_tgt = 8L) %>%
      mutate(
        med_diff = med_tgt - med_ref,
        u = rs_tgt - n_tgt * (n_tgt + 1) / 2,
        auc = u / (n_ref * n_tgt)
      ) %>%
      # Include the number of ties in the computation which
      # were computed in the previous step
      {
        if (nrow(cur_test_groups_ties) > 0) {
          left_join(., cur_test_groups_ties, by = c("marker_1", "marker_2"))
        } else {
          mutate(., nties_const = 0)
        }
      } %>%
      mutate(z = u - (n_ref * n_tgt) / 2) %>%
      mutate(
        sigma = sqrt(
          (n_ref * n_tgt / 12) *
            ((n_ref + n_tgt + 1) - nties_const /
              ((n_ref + n_tgt) * (n_ref + n_tgt - 1))
            )
        )
      ) %>%
      # Formula for two.sided test which is currently the only option
      mutate(z = (z - sign(z) * 0.5) / sigma) %>%
      mutate(p_val = 2 * pmin(pnorm(z), pnorm(z, lower.tail = FALSE))) %>%
      collect() %>%
      na.omit()

    # Tidy up the results
    diff_prox_res[[key]] <- tibble(
      data_type = proximity_metric,
      target = cur_keys[, contrast_column, drop = TRUE],
      reference = reference,
      n_tgt = cur_test_groups_u$n_tgt,
      n_ref = cur_test_groups_u$n_ref,
      diff_median = cur_test_groups_u$med_diff,
      statistic = cur_test_groups_u$u,
      auc = cur_test_groups_u$auc,
      p = cur_test_groups_u$p_val,
      alternative = "two.sided",
      marker_1 = cur_test_groups_u$marker_1,
      marker_2 = cur_test_groups_u$marker_2
    ) %>%
      filter(!is.na(n_ref), !is.na(n_tgt)) %>%
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
  verbose = TRUE,
  ...
) {
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
  group_data <- group_data %>%
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
    verbose = verbose
  )

  # Create a list to store the results in
  keys <- apply(test_groups_keys, 1, paste, collapse = "_")
  diff_prox_res <- rep(list(NULL), length(keys)) %>%
    set_names(keys)

  # iterate over group pairs
  for (i in seq_len(nrow(test_groups_keys))) {
    if (verbose && check_global_verbosity()) {
      cli_progress_update()
    }

    # Fetch the current group keys
    cur_keys <- test_groups_keys[i, ]
    key <- as.character(cur_keys) %>% paste(collapse = "_")

    # Filter data to the current comparison (contrast + group_vars)
    cur_test_groups <- .filter_cur_group(
      test_groups = test_groups,
      cur_keys = cur_keys,
      contrast_column = contrast_column,
      reference = reference,
      group_vars = group_vars
    )

    components_reference <- cur_test_groups %>%
      filter(!!sym(contrast_column) == reference) %>%
      pull(component)
    components_target <- cur_test_groups %>%
      filter(!!sym(contrast_column) != reference) %>%
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

    de_results <- .wilcox_de_test(
      object,
      cells_1 = components_target,
      cells_2 = components_reference,
      min_diff = diff_threshold,
      min_pct = min_pct,
      min_diff_pct = min_diff_pct
    )

    # Filters in .wilcox_de_test() can leave no pairs. Skip the comparison
    # instead of building a zero-row tibble with recycled columns.
    if (nrow(de_results) == 0) {
      next
    }

    # Tidy up the results
    diff_prox_res[[key]] <- tibble(
      data_type = proximity_metric,
      target = cur_keys[, contrast_column, drop = TRUE],
      reference = reference,
      pct_tgt = de_results$pct_1,
      pct_ref = de_results$pct_2,
      diff_median = de_results$difference,
      p = de_results$p_val,
      alternative = "two.sided",
      pair = rownames(de_results)
    ) %>%
      .append_group_vars(cur_keys = cur_keys, group_vars = group_vars)
  }

  if (verbose && check_global_verbosity()) {
    cli_progress_done()
  }

  valid_results <- sapply(diff_prox_res, function(x) !is.null(x))
  if (!any(valid_results)) {
    cli::cli_abort("No valid results were generated. Please check your input data and parameters.")
  }

  .finalize_da_results(diff_prox_res, p_adjust_method, group_vars) %>%
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

  # Use default assay if assay = NULL
  assay <- .validate_or_set_assay(object, assay)

  method <- match.arg(method, choices = c("seurat", "legacy"))

  # Fetch proximity scores
  proximity_data <- .load_proximity_for_da(
    object = object,
    assay = assay,
    lazy = lazy,
    min_exp_join_count = min_exp_join_count,
    proximity_metric = proximity_metric,
    metric_type = metric_type,
    method = method
  )

  # Add group data
  group_data <- .select_group_data(object[[]], contrast_column, group_vars)
  if (method == "seurat") {
    group_data <- group_data %>%
      column_to_rownames("component")
  }

  # Add group data to colocalization table
  if (method == "legacy") {
    proximity_data <- proximity_data %>%
      left_join(y = group_data, copy = TRUE, by = "component") %>%
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

#' Utility function to compute the median difference and percentage
#' of cells expressing a feature in two groups.
#'
#' @param object A matrix where rows are features (e.g. marker pairs) and columns are cells.
#' @param cells_1 A vector of cell names for group 1
#' @param cells_2 A vector of cell names for group 2
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
  data_1 <- sparseMatrixStats::rowMedians(object[, cells_1, drop = FALSE])
  data_2 <- sparseMatrixStats::rowMedians(object[, cells_2, drop = FALSE])
  med_diff <- (data_1 - data_2)
  fc_results <- as.data.frame(x = cbind(med_diff, pct_1, pct_2))
  colnames(fc_results) <- c("difference", "pct_1", "pct_2")
  return(fc_results)
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
#'
#' @noRd
.wilcox_de_test <- function(
  data_use,
  cells_1,
  cells_2,
  min_diff = 0.01,
  min_pct = 0.01,
  min_diff_pct = -Inf
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
  features_diff <- names(which(abs(total_diff) >= min_diff))

  features <- intersect(features, features_diff)
  if (length(features) == 0) {
    cli::cli_warn("No features pass min_diff threshold; returning empty data.frame")
    return(fc_results[features, ])
  }

  data_use <- data_use[features, c(cells_1, cells_2), drop = FALSE]

  group_info <- data.frame(
    row.names = c(cells_1, cells_2),
    group = factor(c(
      rep("Group1", length(cells_1)),
      rep("Group2", length(cells_2))
    ))
  )

  data_use <- data_use[, rownames(group_info), drop = FALSE]

  if (requireNamespace("limma", quietly = TRUE)) {
    j <- seq_along(cells_1)
    p_val <- sapply(seq_len(nrow(data_use)), function(x) {
      min(2 * min(limma::rankSumTestWithCorrelation(
        index = j,
        statistics = data_use[x, ]
      )), 1)
    })
  } else {
    rlang::inform(
      c(
        "i" = "For a faster implementation of the Wilcoxon Rank Sum Test,
      please install the limma package: pak::pak('limma')",
        "v" = "This message will only appear once per R session."
      ),
      .frequency = "once",
      .frequency_id = "limma_not_installed"
    )
    p_val <- sapply(seq_len(nrow(data_use)), function(x) {
      x1 <- as.numeric(data_use[x, cells_1])
      x2 <- as.numeric(data_use[x, cells_2])
      suppressWarnings(stats::wilcox.test(x1, x2, exact = FALSE)$p.value)
    })
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
#' be included. Pairs below the threshold are treated as missing. With
#' \code{method = "seurat"}, missing scores are set to 0.
#' @param diff_threshold Minimum difference in the proximity metric required to
#' test a marker pair. Used when \code{method = "seurat"}.
#' @param min_pct Minimum fraction of cells in either group with a non-zero
#' score. Used when \code{method = "seurat"}.
#' @param min_diff_pct Minimum difference in the fraction of cells with a
#' non-zero score. Used when \code{method = "seurat"}.
#' @param method One of \code{"seurat"} or \code{"legacy"}. Passed through to
#' \code{\link{DifferentialProximityAnalysis}}.
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
  proximity_metric = "log2_ratio",
  metric_type = c("all", "self", "co"),
  backend = c("dplyr", "data.table"),
  method = c("seurat", "legacy"),
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
    proximity_metric = proximity_metric,
    verbose = verbose,
    call = call
  )

  assay <- .validate_or_set_assay(object, assay, call = call)
  metric_type <- match.arg(metric_type, choices = c("all", "self", "co"))
  backend <- match.arg(backend, choices = c("dplyr", "data.table"))
  method <- match.arg(method, choices = c("seurat", "legacy"))
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

  proximity_data <- .load_proximity_for_da(
    object = object,
    assay = assay,
    lazy = lazy,
    min_exp_join_count = min_exp_join_count,
    proximity_metric = proximity_metric,
    metric_type = metric_type,
    method = method
  )
  dots <- list(...)

  if (method == "seurat") {
    cell_labels <- .labels_for_cells(label_chr, colnames(proximity_data), group_by)
    contrast_col <- ".pxl_one_vs_rest"
    results <- .run_one_vs_rest_loop(
      label_chr = cell_labels,
      idents = idents,
      rest_label = rest_label,
      group_by = group_by,
      min_cells_per_group = min_cells_per_group,
      verbose = verbose,
      call = call,
      run_one = function(ident) {
        group_data <- data.frame(
          .pxl_one_vs_rest = .assign_one_vs_rest(
            cell_labels,
            ident,
            rest_label
          ),
          row.names = colnames(proximity_data),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        names(group_data) <- contrast_col
        args <- list(
          object = proximity_data,
          group_data = group_data,
          contrast_column = contrast_col,
          reference = rest_label,
          targets = ident,
          proximity_metric = proximity_metric,
          p_adjust_method = p_adjust_method,
          diff_threshold = diff_threshold,
          min_pct = min_pct,
          min_diff_pct = min_diff_pct,
          min_cells_per_group = min_cells_per_group,
          verbose = verbose
        )
        do.call(DifferentialProximityAnalysis, c(args, dots))
      }
    )
  } else {
    prepared <- .prepare_legacy_group_table(proximity_data, label_chr)
    cell_labels <- .labels_for_cells(
      label_chr,
      unique(as.character(prepared$data$component)),
      group_by
    )
    contrast_col <- .safe_colname(colnames(prepared$data), ".pxl_one_vs_rest")
    results <- .run_one_vs_rest_loop(
      label_chr = cell_labels,
      idents = idents,
      rest_label = rest_label,
      group_by = group_by,
      min_cells_per_group = min_cells_per_group,
      verbose = verbose,
      call = call,
      run_one = function(ident) {
        tab <- prepared$data
        tab[[contrast_col]] <- .assign_one_vs_rest(
          tab[[prepared$group_col]],
          ident,
          rest_label
        )
        tab <- tab %>% filter(!is.na(.data[[contrast_col]]))
        args <- list(
          object = tab,
          contrast_column = contrast_col,
          reference = rest_label,
          targets = ident,
          proximity_metric = proximity_metric,
          metric_type = metric_type,
          backend = backend,
          p_adjust_method = p_adjust_method,
          min_cells_per_group = min_cells_per_group,
          verbose = verbose
        )
        do.call(DifferentialProximityAnalysis, c(args, dots))
      }
    )
  }

  results
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
#' ) %>%
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
  proximity_metric = "log2_ratio",
  metric_type = c("all", "self", "co"),
  backend = c("dplyr", "data.table"),
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
) {
  call <- caller_env()
  .reject_find_all_proximity_args(..., call = call)
  assert_single_value(group_by, type = "string", call = call)
  .validate_find_all_common_args(
    min_cells_per_group = min_cells_per_group,
    proximity_metric = proximity_metric,
    verbose = verbose,
    call = call
  )
  metric_type <- match.arg(metric_type, choices = c("all", "self", "co"))
  backend <- match.arg(backend, choices = c("dplyr", "data.table"))
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
  object <- as_tibble(object) %>% ungroup()
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
    object <- object %>% filter(!is.na(component))
  }

  labels_by_cell <- object %>%
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

  group_col <- .safe_colname(colnames(object), ".pxl_group_label")
  contrast_col <- .safe_colname(colnames(object), ".pxl_one_vs_rest")
  object[[group_col]] <- unname(label_chr[as.character(object$component)])
  dots <- list(...)

  .run_one_vs_rest_loop(
    label_chr = unname(label_chr),
    idents = idents,
    rest_label = rest_label,
    group_by = group_by,
    min_cells_per_group = min_cells_per_group,
    verbose = verbose,
    call = call,
    run_one = function(ident) {
      tab <- object
      tab[[contrast_col]] <- .assign_one_vs_rest(
        tab[[group_col]],
        ident,
        rest_label
      )
      tab <- tab %>% filter(!is.na(.data[[contrast_col]]))
      args <- list(
        object = tab,
        contrast_column = contrast_col,
        reference = rest_label,
        targets = ident,
        proximity_metric = proximity_metric,
        metric_type = metric_type,
        backend = backend,
        p_adjust_method = p_adjust_method,
        min_cells_per_group = min_cells_per_group,
        verbose = verbose
      )
      do.call(DifferentialProximityAnalysis, c(args, dots))
    }
  )
}


#' Load proximity scores once for differential testing.
#'
#' @return For \code{method = "seurat"}, a sparse matrix of marker pairs by
#'   cells. For \code{method = "legacy"}, a proximity table without group
#'   columns. \code{metric_type} filtering matches
#'   \code{DifferentialProximityAnalysis.Seurat()}.
#'
#' @noRd
.load_proximity_for_da <- function(
  object,
  assay,
  lazy,
  min_exp_join_count,
  proximity_metric,
  metric_type,
  method
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
    self = proximity_data %>% filter(marker_1 == marker_2),
    co = proximity_data %>% filter(marker_1 != marker_2)
  ) %>%
    filter(join_count_expected_mean >= min_exp_join_count)
  proximity_data <- proximity_data %>% compute()
  if (method == "seurat") {
    proximity_data <- proximity_data %>%
      compute() %>%
      ProximityScoresToAssay(values_from = proximity_metric)
    missing_components <- setdiff(colnames(object), colnames(proximity_data))
    if (length(missing_components) > 0) {
      m_missing <- Matrix::rsparsematrix(
        nrow = nrow(proximity_data),
        ncol = length(missing_components),
        density = 0
      )
      rownames(m_missing) <- rownames(proximity_data)
      colnames(m_missing) <- missing_components
      proximity_data <- cbind(proximity_data, m_missing)
      proximity_data <- proximity_data[, colnames(object)]
    }
  }
  proximity_data
}


#' @noRd
.reject_find_all_proximity_args <- function(..., call = caller_env()) {
  blocked <- c(
    "contrast_column", "reference", "targets", "group_vars", "group_data"
  )
  supplied <- intersect(names(list(...)), blocked)
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
}


#' @noRd
.validate_find_all_common_args <- function(
  min_cells_per_group,
  proximity_metric,
  verbose,
  call = caller_env()
) {
  assert_single_value(min_cells_per_group, type = "numeric", call = call)
  if (is.na(min_cells_per_group) || min_cells_per_group < 0) {
    cli::cli_abort(
      c("x" = "{.arg min_cells_per_group} must be a non-negative number."),
      call = call
    )
  }
  assert_single_value(proximity_metric, type = "string", call = call)
  assert_single_value(verbose, type = "bool", call = call)
}


#' @noRd
.coerce_idents <- function(idents) {
  if (is.factor(idents) || is.numeric(idents)) {
    return(as.character(idents))
  }
  idents
}


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


#' @noRd
.assign_one_vs_rest <- function(labels, ident, rest_label) {
  labels <- as.character(labels)
  contrast <- rep(NA_character_, length(labels))
  known <- !is.na(labels)
  contrast[known & labels == ident] <- ident
  contrast[known & labels != ident] <- rest_label
  contrast
}


#' @noRd
.safe_colname <- function(existing, base) {
  candidate <- base
  while (candidate %in% existing) {
    candidate <- paste0(candidate, "_")
  }
  candidate
}


#' @noRd
.warn_missing_group_labels <- function(label_chr, group_by) {
  n_missing <- sum(is.na(label_chr))
  if (n_missing > 0) {
    cli::cli_warn(
      "Excluding {n_missing} cell{?s} with missing {.val {group_by}} labels."
    )
  }
}


#' Align named cell labels to \code{cells}, warning when metadata is missing.
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


#' @noRd
.prepare_legacy_group_table <- function(proximity_data, label_chr) {
  proximity_data <- proximity_data %>%
    collect() %>%
    as_tibble() %>%
    ungroup()
  group_col <- .safe_colname(colnames(proximity_data), ".pxl_group_label")
  proximity_data[[group_col]] <- unname(
    label_chr[as.character(proximity_data$component)]
  )
  list(data = proximity_data, group_col = group_col)
}


#' @noRd
.one_vs_rest_counts <- function(label_chr, ident) {
  list(
    n_tgt = sum(label_chr == ident, na.rm = TRUE),
    n_ref = sum(!is.na(label_chr) & label_chr != ident)
  )
}


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


#' @noRd
.is_skippable_da_error <- function(err) {
  msg <- conditionMessage(err)
  grepl(
    paste(
      "Found no groups with at least",
      "Found no valid target data",
      "No valid results were generated",
      sep = "|"
    ),
    msg
  )
}


#' Run one \code{DifferentialProximityAnalysis()} call per group level.
#'
#' @param run_one A function of one argument, the level to test.
#'
#' @noRd
.run_one_vs_rest_loop <- function(
  label_chr,
  idents,
  rest_label,
  group_by,
  min_cells_per_group,
  verbose,
  run_one,
  call = caller_env()
) {
  if (verbose && check_global_verbosity()) {
    n_tests <- length(idents)
    cli_alert_info(
      paste0(
        "Running one-versus-rest proximity tests for {n_tests} level{?s} ",
        "of {.val {group_by}}. The pooled remainder is labeled ",
        "{.val {rest_label}}."
      )
    )
  }

  results <- vector("list", length(idents))
  for (i in seq_along(idents)) {
    ident <- idents[[i]]
    counts <- .one_vs_rest_counts(label_chr, ident)
    if (
      counts$n_tgt < min_cells_per_group ||
        counts$n_ref < min_cells_per_group
    ) {
      .warn_skip_one_vs_rest(
        ident,
        rest_label,
        counts$n_tgt,
        counts$n_ref,
        min_cells_per_group
      )
      next
    }

    n_tgt <- counts$n_tgt
    n_ref <- counts$n_ref
    if (verbose && check_global_verbosity()) {
      cli_alert_info(
        paste0(
          "Testing {.val {ident}} vs {.val {rest_label}} ",
          "({n_tgt} vs {n_ref} cells)."
        )
      )
    }

    result <- tryCatch(run_one(ident), error = function(e) e)
    if (inherits(result, "error")) {
      if (.is_skippable_da_error(result)) {
        cli::cli_warn(
          paste0(
            "Skipping {.val {ident}} vs {.val {rest_label}}: ",
            "{conditionMessage(result)}"
          )
        )
        next
      }
      rlang::cnd_signal(result)
    }
    results[[i]] <- result
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
  out
}
