#' @include generics.R
NULL

.datatable.aware <- TRUE


#' @rdname DifferentialProximityAnalysis
#' @method DifferentialProximityAnalysis data.frame
#'
#' @examples
#' # TODO: Update examples with real data
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
  proximity_metric = "join_count_z",
  metric_type = c("all", "self", "co"),
  backend = c("dplyr", "data.table"),
  min_n_obs = 0,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
) {
  .validate_da_input(
    object, contrast_column, reference, targets, group_vars,
    proximity_metric, min_n_obs, FALSE,
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
    abort("No data found for the specified metric type.")
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

  if (min_n_obs > 0) {
    # Filter test groups by minimum number of observations allowed
    test_groups <- test_groups %>%
      mutate(n = n()) %>%
      filter(n > min_n_obs) %>%
      ungroup(!!sym(contrast_column)) %>%
      mutate(
        ref_n = sum(!!sym(contrast_column) == reference),
        target_n = sum(!!sym(contrast_column) %in% targets)
      ) %>%
      filter(ref_n > 0, target_n > 0) %>%
      group_by(!!sym(contrast_column), .add = TRUE)
    if (nrow(test_groups) == 0) {
      abort(glue(
        "Found no groups with at least {min_n_obs} observations."
      ))
    }
  }


  # Get group keys
  test_groups_keys <- test_groups %>%
    group_keys() %>%
    select(-marker_1, -marker_2) %>%
    filter(!!sym(contrast_column) %in% targets) %>%
    distinct() %>%
    mutate_all(as.character)

  # Print information about the test set up
  if (verbose && check_global_verbosity()) {
    cli_alert_info("Computing Running Wilcoxon rank-sum test for each marker pair across the following comparisons:")
    cat_line()
    cli_div(theme = list(ul = list(`margin-left` = 2, before = "")))
    if (!is.null(group_vars)) {
      ctrst <- test_groups_keys %>% pull(!!sym(contrast_column))
      grps <- test_groups_keys %>%
        select(-!!sym(contrast_column)) %>%
        apply(1, paste, collapse = ", ")
      grps_tgts <- glue::glue("{ctrst} ({grps}) vs {reference} ({grps})")
      cli_ul(grps_tgts, id = "ul")
    } else {
      grps <- targets
      cli_ul(glue::glue("{targets} vs {reference}"), id = "ul")
    }
    cli_end(id = "ul")
    # Set up a progress bar using cli
    cli_progress_bar(
      glue("Running Wilcoxon rank-sum test:"),
      total = length(grps),
      clear = FALSE
    )
  }

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
    # Filter data to include a single comparison
    groups_keep <- c(reference, cur_keys[, contrast_column, drop = TRUE])
    cur_test_groups <- test_groups %>%
      filter(!!sym(contrast_column) %in% groups_keep)

    # Filter by group_vars if specified
    if (!is.null(group_vars)) {
      for (nm in group_vars) {
        cur_test_groups <- cur_test_groups %>%
          filter(!!sym(nm) %in% cur_keys[, nm, drop = TRUE])
      }
    }

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
      collect()

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
      method = "Wilcoxon rank-sum test test",
      alternative = "two.sided",
      marker_1 = cur_test_groups_u$marker_1,
      marker_2 = cur_test_groups_u$marker_2
    )

    # Add additional group columns (stored in cur_keys) if needed
    if (!is.null(group_vars)) {
      for (nm in group_vars) {
        diff_prox_res[[key]] <- diff_prox_res[[key]] %>%
          mutate(!!sym(nm) := cur_keys[, nm, drop = TRUE])
      }
    }
  }

  if (verbose && check_global_verbosity()) {
    cli_progress_done()
  }

  diff_prox_res <- bind_rows(diff_prox_res)

  # Adjust p-values
  diff_prox_res <- diff_prox_res %>%
    mutate(p_adj = p.adjust(p, p_adjust_method)) %>%
    relocate(p_adj, .after = "p")

  return(diff_prox_res)
}


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
  proximity_metric = "join_count_z",
  metric_type = c("all", "self", "co"),
  min_n_obs = 0,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
) {
  # Validate input parameters
  stopifnot(
    "'contrast_column' must be available in Seurat object meta.data" =
      contrast_column %in% colnames(object[[]])
  )

  # Use default assay if assay = NULL
  assay <- .validate_or_set_assay(object, assay)

  # Fetch proximity scores
  proximity_data <- object[[assay]]@proximity

  # Add group data
  group_data <- .select_group_data(object[[]], contrast_column, group_vars)

  # Add group data to colocalization table
  proximity_data <- proximity_data %>%
    left_join(y = group_data, by = "component")

  # Run differential proximity analysis
  proximity_test_results <- DifferentialProximityAnalysis(
    proximity_data,
    targets = targets,
    reference = reference,
    contrast_column = contrast_column,
    group_vars = group_vars,
    proximity_metric = proximity_metric,
    metric_type = metric_type,
    min_n_obs = min_n_obs,
    p_adjust_method = p_adjust_method,
    verbose = verbose
  )

  return(proximity_test_results)
}
