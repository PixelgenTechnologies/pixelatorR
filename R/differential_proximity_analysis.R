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

  .finalize_da_results(diff_prox_res, p_adjust_method, group_vars) %>%
    tidyr::separate(pair, into = c("marker_1", "marker_2"), sep = ":")
}


#' @param lazy If TRUE, the proximity scores will be loaded lazily and filtered using the
#' `duckdb` backend.
#' @param min_exp_join_count Minimum number of join counts required for a marker pair to be
#' included in the analysis.
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
  proximity_data <- ProximityScores(object[[assay]], assay = assay, lazy = lazy)
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
      ProximityScoresToAssay()
    # Add missing columns
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
      min_n_obs = min_n_obs,
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
    Matrix::rowSums(object[, cells_1, drop = FALSE] > 0) /
      length(x = cells_1),
    digits = 3
  )
  pct_2 <- round(
    Matrix::rowSums(object[, cells_2, drop = FALSE] > 0) /
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
    warning("No features pass min_diff threshold; returning empty data.frame")
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

  presto_check <- requireNamespace("presto", quietly = TRUE)
  if (presto_check[1]) {
    data_use <- data_use[, rownames(group_info), drop = FALSE]
    res <- presto::wilcoxauc(
      X = data_use,
      y = group_info[
        ,
        "group"
      ]
    )
    res <- res[1:(nrow(x = res) / 2), ]
    p_val <- res$pval
  } else {
    rlang::inform(
      c(
        "i" = "For a faster implementation of the Wilcoxon Rank Sum Test,
      please install the presto package: pak::pak('immunogenomics/presto'",
        "v" = "This message will only appear once per R session."
      ),
      .frequency = "once",
      .frequency_id = "presto_not_installed"
    )
    j <- seq_along(cells_1)
    data_use <- data_use[, rownames(group_info), drop = FALSE]
    p_val <- sapply(seq_len(nrow(data_use)), function(x) {
      return(min(2 * min(limma::rankSumTestWithCorrelation(
        index = j,
        statistics = data_use[x, ]
      )), 1))
    })
  }
  de_results <- cbind(
    data.frame(p_val, row.names = rownames(data_use)),
    fc_results[rownames(data_use), ]
  )
  return(de_results)
}
