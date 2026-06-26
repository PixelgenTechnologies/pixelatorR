#' Select columns from a \code{tbl_df} for grouping
#'
#' @returns A \code{tbl_df} object
#'
#' @noRd
.select_group_data <- function(object, contrast_column, group_vars = NULL) {
  group_data <- object %>%
    {
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
    } %>%
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
  test_groups <- data %>% group_by(pick(all_of(group_cols)))

  keys <- test_groups %>% group_keys()
  if (!is.null(drop_cols)) {
    drop_cols <- intersect(drop_cols, colnames(keys))
    if (length(drop_cols) > 0) {
      keys <- keys %>% select(-all_of(drop_cols))
    }
  }
  keys <- keys %>%
    filter(!!sym(contrast_column) %in% targets) %>%
    distinct() %>%
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
#' @param verbose Logical. If \code{FALSE}, this function is a no-op.
#' @param test_label Character. Name of the test, used in messages.
#'
#' @noRd
.print_da_setup <- function(
  keys,
  contrast_column,
  reference,
  targets,
  group_vars,
  verbose,
  test_label = "Wilcoxon rank-sum test",
  env = rlang::caller_env()
) {
  if (!(verbose && check_global_verbosity())) {
    return(invisible(NULL))
  }

  cli_alert_info("Running {test_label} across the following comparisons:")
  cat_line()
  cli_div(theme = list(ul = list(`margin-left` = 2, before = "")))
  if (!is.null(group_vars)) {
    ctrst <- keys %>% pull(!!sym(contrast_column))
    grps <- keys %>%
      select(-!!sym(contrast_column)) %>%
      apply(1, paste, collapse = ", ")
    cli_ul(glue::glue("{ctrst} ({grps}) vs {reference} ({grps})"), id = "ul")
  } else {
    cli_ul(glue::glue("{targets} vs {reference}"), id = "ul")
  }
  cli_end(id = "ul")

  cli_progress_bar(
    glue("Running {test_label}:"),
    total = nrow(keys),
    clear = FALSE,
    .envir = env
  )
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
  out <- test_groups %>%
    filter(!!sym(contrast_column) %in% groups_keep)

  if (!is.null(group_vars)) {
    for (nm in group_vars) {
      out <- out %>%
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
    res <- res %>%
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
.finalize_da_results <- function(res_list, p_adjust_method) {
  bind_rows(res_list) %>%
    mutate(p_adj = p.adjust(p, p_adjust_method)) %>%
    relocate(p_adj, .after = "p")
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
