#' Select columns from a \code{tbl_df} for grouping
#'
#' @returns A \code{tbl_df} object
#'
#' @noRd
.select_group_data <- function(object, contrast_column, group_vars = NULL) {
  group_data <- object |>
    (\(.) {
      if (!is.null(group_vars)) {
        assert_vector(group_vars, type = "character", n = 1)
        if (!all(group_vars %in% colnames(.))) {
          cli::cli_abort(
            c("x" = "All {.var group_vars} must be valid column names")
          )
        }
        select(., all_of(c(contrast_column, group_vars)))
      } else {
        select(., all_of(contrast_column))
      }
    })() |>
    as_tibble(rownames = "component")
}

#' Prints the number of groups and comparisons for a differential analysis
#'
#' @noRd
.print_da_group_strategy <- function(test_groups_keys, targets, reference) {
  if (length(targets) > 1) {
    cli_alert_info(glue(
      "Running {nrow(test_groups_keys)} tests for each ",
      "of the following comparisons:\n"
    ))
  } else {
    cli_alert_info(glue(
      "Running {nrow(test_groups_keys)} tests for the following comparison:\n"
    ))
  }
  print(glue("  - {paste(targets, reference, sep = ' vs ')}"))
}

# ---------------------------------------------------------------------------
# Shared helpers used by `DifferentialProximityAnalysis()` (and friends)
# to avoid duplicated logic between the data.frame and Matrix methods.
# ---------------------------------------------------------------------------

#' Compute grouped test data and the corresponding group keys
#'
#' Groups \code{data} by \code{group_cols} and returns both the grouped tibble
#' and a tidy keys tibble filtered to the requested \code{targets}.
#'
#' @param data A data.frame-like object.
#' @param group_cols Character vector of column names to group by.
#' @param contrast_column The contrast column name.
#' @param targets Target levels of the contrast column to keep.
#' @param drop_cols Optional character vector of columns to drop from the
#'   keys tibble (e.g. \code{c("marker_1", "marker_2")} for the proximity
#'   data.frame method).
#'
#' @return A list with elements \code{test_groups} (grouped tibble) and
#'   \code{keys} (a tidy tibble of unique group key combinations).
#'
#' @noRd
.compute_test_groups_keys <- function(
  data,
  group_cols,
  contrast_column,
  targets,
  drop_cols = NULL
) {
  test_groups <- data |> group_by(pick(all_of(group_cols)))

  keys <- test_groups |> group_keys()
  if (!is.null(drop_cols)) {
    drop_cols <- intersect(drop_cols, colnames(keys))
    if (length(drop_cols) > 0) {
      keys <- keys |> select(-all_of(drop_cols))
    }
  }
  keys <- keys |>
    filter(!!sym(contrast_column) %in% targets) |>
    distinct() |>
    mutate_all(as.character)

  list(test_groups = test_groups, keys = keys)
}

#' Print the test setup and start a cli progress bar
#'
#' Emits a \pkg{cli} message describing the comparisons that will be run and
#' starts a progress bar tracking iterations across \code{keys}.
#'
#' @param keys A tibble of group keys (e.g. from
#'   \code{.compute_test_groups_keys()}).
#' @param contrast_column The contrast column name.
#' @param reference The reference level.
#' @param targets Target levels of the contrast column.
#' @param group_vars Optional grouping variables.
#' @param verbose Logical. Whether to print messages.
#' @param test_label Character. Name of the test, used in messages.
#'
#' @noRd
#' Print the test setup and start a cli progress bar
#'
#' Emits a formatted \pkg{cli} summary table of the comparisons that will be run
#' (including per-group cell counts and skipped status when \code{test_groups}
#' is supplied) and starts a transient progress bar tracking iterations across
#' \code{keys}.
#'
#' @param keys A tibble of group keys (e.g. from
#'   \code{.compute_test_groups_keys()}).
#' @param contrast_column The contrast column name.
#' @param reference The reference level.
#' @param targets Target levels of the contrast column.
#' @param group_vars Optional grouping variables.
#' @param verbose Logical. Whether to print messages.
#' @param test_label Character. Name of the test, used in messages.
#' @param test_groups Optional grouped tibble used to compute per-comparison
#'   target and reference cell counts.
#' @param min_cells_per_group Optional minimum number of cells per group used to
#'   mark skipped comparisons in the overview table.
#'
#' @noRd
.print_da_setup <- function(
  keys,
  contrast_column,
  reference,
  targets,
  group_vars,
  verbose,
  test_label = "Wilcoxon rank-sum",
  test_groups = NULL,
  min_cells_per_group = 0,
  env = rlang::caller_env()
) {
  if (!(verbose && check_global_verbosity())) {
    return(invisible(NULL))
  }

  n_rows <- nrow(keys)
  has_group_vars <- !is.null(group_vars) && length(group_vars) > 0L

  contrasts <- paste0(keys[[contrast_column]], " vs ", reference)
  grp_labels <- if (has_group_vars) {
    apply(keys[, group_vars, drop = FALSE], 1, paste, collapse = ", ")
  } else {
    rep("-", n_rows)
  }

  has_counts <- !is.null(test_groups) && "component" %in% colnames(test_groups)
  if (has_counts) {
    counts_mat <- matrix(
      vapply(seq_len(n_rows), function(i) {
        cur_keys <- keys[i, ]
        cur_grp <- .filter_cur_group(
          test_groups = test_groups,
          cur_keys = cur_keys,
          contrast_column = contrast_column,
          reference = reference,
          group_vars = group_vars
        )
        c_col <- cur_grp[[contrast_column]]
        comp_col <- cur_grp[["component"]]
        n_ref <- length(unique(comp_col[c_col == reference]))
        n_tgt <- length(unique(comp_col[c_col != reference]))
        c(n_tgt = n_tgt, n_ref = n_ref)
      }, integer(2)),
      nrow = 2L,
      dimnames = list(c("n_tgt", "n_ref"), NULL)
    )
    n_tgt_vec <- counts_mat["n_tgt", ]
    n_ref_vec <- counts_mat["n_ref", ]
    tgt_str <- format(n_tgt_vec, big.mark = ",", trim = TRUE)
    ref_str <- format(n_ref_vec, big.mark = ",", trim = TRUE)
    status_str <- ifelse(
      n_tgt_vec < min_cells_per_group | n_ref_vec < min_cells_per_group,
      paste0(cli::col_yellow(cli::symbol$cross), " skipped (< ", min_cells_per_group, " cells)"),
      cli::col_green(cli::symbol$tick)
    )
  }

  cli::cli_rule(left = paste0("Differential Proximity Analysis (", test_label, ")"))

  # Display contrast header only on first row of each consecutive contrast block
  contrast_disp <- ifelse(
    seq_len(n_rows) == 1L | contrasts != dplyr::lag(contrasts, default = ""),
    contrasts,
    ""
  )

  w_contrast <- max(nchar( c("Contrast", contrast_disp) ))
  w_group <- if (has_group_vars) max(nchar(c("Group", grp_labels))) else 0L

  if (has_counts) {
    w_tgt <- max(nchar(c("Target", tgt_str)))
    w_ref <- max(nchar(c("Reference", ref_str)))
    if (has_group_vars) {
      hdr <- paste0(
        "  ",
        cli::style_bold(format("Contrast", width = w_contrast)), "  ",
        cli::style_bold(format("Group", width = w_group)), "  ",
        cli::style_bold(format("Target", width = w_tgt, justify = "right")), "  ",
        cli::style_bold(format("Reference", width = w_ref, justify = "right")), "  ",
        cli::style_bold("Status")
      )
    } else {
      hdr <- paste0(
        "  ",
        cli::style_bold(format("Contrast", width = w_contrast)), "  ",
        cli::style_bold(format("Target", width = w_tgt, justify = "right")), "  ",
        cli::style_bold(format("Reference", width = w_ref, justify = "right")), "  ",
        cli::style_bold("Status")
      )
    }
    cli::cli_verbatim(hdr)
    for (i in seq_len(n_rows)) {
      if (has_group_vars) {
        line <- paste0(
          "  ",
          format(contrast_disp[i], width = w_contrast), "  ",
          format(grp_labels[i], width = w_group), "  ",
          format(tgt_str[i], width = w_tgt, justify = "right"), "  ",
          format(ref_str[i], width = w_ref, justify = "right"), "  ",
          status_str[i]
        )
      } else {
        line <- paste0(
          "  ",
          format(contrast_disp[i], width = w_contrast), "  ",
          format(tgt_str[i], width = w_tgt, justify = "right"), "  ",
          format(ref_str[i], width = w_ref, justify = "right"), "  ",
          status_str[i]
        )
      }
      cli::cli_verbatim(line)
    }
  } else {
    if (has_group_vars) {
      hdr <- paste0(
        "  ",
        cli::style_bold(format("Contrast", width = w_contrast)), "  ",
        cli::style_bold(format("Group", width = w_group))
      )
      cli::cli_verbatim(hdr)
      for (i in seq_len(n_rows)) {
        cli::cli_verbatim(paste0(
          "  ",
          format(contrast_disp[i], width = w_contrast), "  ",
          format(grp_labels[i], width = w_group)
        ))
      }
    } else {
      hdr <- paste0("  ", cli::style_bold(format("Contrast", width = w_contrast)))
      cli::cli_verbatim(hdr)
      for (i in seq_len(n_rows)) {
        cli::cli_verbatim(paste0("  ", format(contrast_disp[i], width = w_contrast)))
      }
    }
  }

  cli_progress_bar(
    format = paste0(
      "[{cli::pb_bar}] {cli::pb_percent} | ",
      "{cli::pb_current}/{cli::pb_total} | {cli::pb_status}"
    ),
    total = n_rows,
    clear = TRUE,
    status = "starting",
    .envir = env
  )
}

#' Format a single DifferentialProximityAnalysis comparison status string
#'
#' @noRd
.format_da_comparison_status <- function(
  cur_keys,
  contrast_column,
  reference,
  group_vars = NULL,
  n_tgt = NULL,
  n_ref = NULL
) {
  target <- as.character(cur_keys[, contrast_column, drop = TRUE])
  grp_str <- if (!is.null(group_vars) && length(group_vars) > 0L) {
    paste(as.character(cur_keys[, group_vars, drop = TRUE]), collapse = ", ")
  } else {
    NULL
  }
  cells_str <- if (!is.null(n_tgt) && !is.null(n_ref)) {
    paste0(n_tgt, " vs ", n_ref, " cells")
  } else {
    NULL
  }
  details <- paste(c(grp_str, cells_str), collapse = ": ")
  if (nzchar(details)) {
    paste0(target, " vs ", reference, " (", details, ")")
  } else {
    paste0(target, " vs ", reference)
  }
}

#' Abort with an informative message when no DifferentialProximityAnalysis comparisons produce results
#'
#' @noRd
.abort_no_valid_da_results <- function(
  n_total,
  n_passed_cells,
  min_cells_per_group,
  contrast_column,
  group_vars = NULL,
  diff_threshold = NULL,
  min_pct = NULL,
  min_diff_pct = NULL,
  only_pos = FALSE,
  call = rlang::caller_env()
) {
  if (n_total == 0L) {
    cli::cli_abort(
      c(
        "x" = "No comparisons found to test.",
        "i" = "Check {.arg contrast_column} ({.val {contrast_column}}), {.arg targets}, and {.arg reference}."
      ),
      call = call
    )
  }

  if (n_passed_cells == 0L) {
    grp_hint <- if (!is.null(group_vars) && length(group_vars) > 0L) {
      paste0(
        "Consider using fewer {.arg group_vars} (currently {.val {group_vars}}) ",
        "or lowering {.arg min_cells_per_group}."
      )
    } else {
      "Consider lowering {.arg min_cells_per_group} or checking group sizes in {.arg contrast_column} ({.val {contrast_column}})."
    }
    cli::cli_abort(
      c(
        "x" = paste0(
          "{?The/All {n_total}} comparison{?s} {?was/were} skipped because every target or ",
          "reference group had fewer than {.arg min_cells_per_group} = ",
          "{.val {min_cells_per_group}} cells."
        ),
        "i" = grp_hint
      ),
      call = call
    )
  }

  n_skipped_cells <- n_total - n_passed_cells
  msg_bullets <- c(
    "x" = paste0(
      "No marker pairs passed the pre-test filtering thresholds across the ",
      "{n_passed_cells} comparison{?s} that met {.arg min_cells_per_group} ",
      "({.val {min_cells_per_group}} cells)."
    )
  )
  if (n_skipped_cells > 0L) {
    msg_bullets <- c(
      msg_bullets,
      "!" = "{n_skipped_cells} of {n_total} comparison{?s} {?was/were} also skipped due to insufficient cells."
    )
  }
  msg_bullets <- c(
    msg_bullets,
    "i" = paste0(
      "Consider relaxing {.arg diff_threshold} ({.val {diff_threshold}}), ",
      "{.arg min_pct} ({.val {min_pct}}), {.arg min_diff_pct} ({.val {min_diff_pct}})",
      if (isTRUE(only_pos)) ", or {.code only_pos = TRUE}." else "."
    )
  )
  cli::cli_abort(msg_bullets, call = call)
}

#' Filter a grouped test tibble to a single comparison
#'
#' Restricts \code{test_groups} to the rows for the current comparison
#' defined by \code{cur_keys}: keeps only the reference and current target
#' levels of \code{contrast_column}, and (if provided) filters to the current
#' \code{group_vars} levels.
#'
#' @param test_groups A grouped tibble.
#' @param cur_keys A one-row tibble of the current group keys.
#' @param contrast_column The contrast column name.
#' @param reference The reference level.
#' @param group_vars Optional grouping variables.
#'
#' @return A filtered (still grouped) tibble.
#'
#' @noRd
.filter_cur_group <- function(
  test_groups,
  cur_keys,
  contrast_column,
  reference,
  group_vars = NULL
) {
  groups_keep <- c(reference, cur_keys[, contrast_column, drop = TRUE])
  out <- test_groups |>
    filter(!!sym(contrast_column) %in% groups_keep)

  if (!is.null(group_vars)) {
    for (nm in group_vars) {
      out <- out |>
        filter(!!sym(nm) %in% cur_keys[, nm, drop = TRUE])
    }
  }
  out
}

#' Append the current group_vars values as columns to a result tibble
#'
#' @param res A tibble of per-comparison results.
#' @param cur_keys A one-row tibble of the current group keys.
#' @param group_vars Optional grouping variables.
#'
#' @return The input tibble with one extra column per \code{group_vars}.
#'
#' @noRd
.append_group_vars <- function(res, cur_keys, group_vars) {
  if (is.null(group_vars)) {
    return(res)
  }
  for (nm in group_vars) {
    res <- res |>
      mutate(!!sym(nm) := cur_keys[, nm, drop = TRUE])
  }
  res
}

#' Combine per-comparison results and adjust p-values
#'
#' @param res_list A list of per-comparison result tibbles.
#' @param p_adjust_method A valid \code{p.adjust} method.
#'
#' @return A single tibble with a \code{p_adj} column inserted after \code{p}.
#'
#' @noRd
.finalize_da_results <- function(res_list, p_adjust_method, group_vars = NULL) {
  bind_rows(res_list) |>
    group_by(pick(all_of(group_vars))) |>
    mutate(p_adj = p.adjust(p, p_adjust_method)) |>
    relocate(p_adj, .after = "p") |>
    ungroup()
}

#' Validate \code{group_data} for matrix-based differential analysis methods
#'
#' Ensures that \code{group_data} is a \code{data.frame}, that
#' \code{contrast_column} and optional \code{group_vars} are columns of it,
#' that \code{rownames(group_data)} match \code{colnames(object)} exactly,
#' and that \code{reference} (plus any \code{targets}) are present in
#' \code{group_data[[contrast_column]]}.
#'
#' @param object A matrix-like object with sample identifiers as column names.
#' @param group_data A \code{data.frame} with rownames matching
#'   \code{colnames(object)}.
#' @param contrast_column,reference,targets,group_vars See
#'   \code{DifferentialProximityAnalysis()}.
#' @param call Environment to use for the error call.
#'
#' @return Invisibly returns \code{NULL}; called for its side effects.
#'
#' @noRd
.validate_matrix_group_data <- function(
  object,
  group_data,
  contrast_column,
  reference,
  targets = NULL,
  group_vars = NULL,
  call = caller_env()
) {
  # group_data must be a data.frame
  if (!inherits(group_data, "data.frame")) {
    cli::cli_abort(
      c("x" = "{.arg group_data} must be a {.cls data.frame}."),
      call = call
    )
  }

  # contrast_column must be a single string present in group_data
  assert_single_value(contrast_column, type = "string", call = call)
  assert_col_in_data(contrast_column, group_data, call = call)
  assert_col_class(
    contrast_column,
    group_data,
    classes = c("character", "factor"),
    call = call
  )

  # group_vars (if specified) must exist in group_data
  if (!is.null(group_vars)) {
    assert_vector(group_vars, type = "character", n = 1, call = call)
    assert_x_in_y(group_vars, colnames(group_data), call = call)
    if (contrast_column %in% group_vars) {
      cli::cli_abort(
        c(
          "x" = "{.arg contrast_column} ({.val {contrast_column}}) cannot be one of {.arg group_vars}."
        ),
        call = call
      )
    }
  }

  # rownames(group_data) must exist and match colnames(object) exactly
  if (is.null(rownames(group_data))) {
    cli::cli_abort(
      c(
        "x" = "{.arg group_data} must have rownames matching {.code colnames(object)}."
      ),
      call = call
    )
  }
  if (is.null(colnames(object))) {
    cli::cli_abort(
      c("x" = "{.arg object} must have column names."),
      call = call
    )
  }
  if (!setequal(rownames(group_data), colnames(object))) {
    n_missing <- length(setdiff(colnames(object), rownames(group_data)))
    n_extra <- length(setdiff(rownames(group_data), colnames(object)))
    cli::cli_abort(
      c(
        "x" = "{.code rownames(group_data)} must match {.code colnames(object)}.",
        "i" = "{n_missing} column(s) of {.arg object} not present in {.code rownames(group_data)}.",
        "i" = "{n_extra} rowname(s) of {.arg group_data} not present in {.code colnames(object)}."
      ),
      call = call
    )
  }

  # reference must be present in the contrast column
  assert_single_value(reference, type = "string", call = call)
  ctrst_vals <- unique(group_data[[contrast_column]])
  if (!reference %in% ctrst_vals) {
    cli::cli_abort(
      c(
        "x" = "Reference group {.val {reference}} must be present
               in the {.val {contrast_column}} column of {.arg group_data}."
      ),
      call = call
    )
  }

  # targets (if specified) must be present in the contrast column
  if (!is.null(targets)) {
    assert_vector(targets, type = "character", n = 1, call = call)
    if (reference %in% targets) {
      cli::cli_abort(
        c(
          "x" = "All {.arg targets} ({.val {targets}}) must be different
                 from {.arg reference} ({.val {reference}})."
        ),
        call = call
      )
    }
    missing_tgt <- setdiff(targets, ctrst_vals)
    if (length(missing_tgt) > 0) {
      cli::cli_abort(
        c(
          "x" = "Target group(s) {.val {missing_tgt}} not found in the
                 {.val {contrast_column}} column of {.arg group_data}."
        ),
        call = call
      )
    }
  }

  invisible(NULL)
}
