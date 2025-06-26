#' @include generics.R
NULL

#' @param name Name of the temporary table in the database containing the filtered proximity scores.
#' @method FilterProximityScores tbl_lazy
#' @rdname FilterProximityScores
#'
#' @examples
#' library(pixelatorR)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' se <- ReadPNA_Seurat(pxl_file)
#' proximity_table <- ProximityScores(se, add_marker_proportions = TRUE)
#'
#' # Filter scores
#' proximity_table_filtered <- FilterProximityScores(
#'   proximity_table,
#'   background_threshold_pct = 0.001
#' )
#'
#' # Rows kept
#' pct_rows_kept <- round(nrow(proximity_table_filtered) / nrow(proximity_table) * 100, digits = 2)
#' glue::glue("Fraction of rows kept: {pct_rows_kept}%")
#'
#' @export
#'
FilterProximityScores.tbl_lazy <- function(
  object,
  background_threshold_pct = NULL,
  background_threshold_count = NULL,
  min_cells_count = NULL,
  name = "filtered_proximity",
  ...
) {
  assert_single_value(background_threshold_pct, "numeric", allow_null = TRUE)
  if (!is.null(background_threshold_pct)) {
    assert_within_limits(background_threshold_pct, limits = c(0, 1))
    if (!all(c("p1", "p2") %in% colnames(object))) {
      cli::cli_abort(
        c(
          "i" = "The columns {.var p1} and {.var p2} are required to filter the proximity scores.",
          "i" = "Make sure to set {.var add_marker_proportions = TRUE} when calling {.fn ProximityScores}.",
          "x" = "Missing columns: {.var p1} or {.var p2}."
        )
      )
    }
  }
  assert_single_value(background_threshold_count, "numeric", allow_null = TRUE)
  if (!is.null(background_threshold_count)) {
    assert_within_limits(background_threshold_count, limits = c(0, Inf))
    if (!all(c("count_1", "count_2") %in% colnames(object))) {
      cli::cli_abort(
        c(
          "i" = "The columns {.var count_1} and {.var count_2} are required to filter the proximity scores.",
          "i" = "Make sure to set {.var add_marker_proportions = TRUE} when calling {.fn ProximityScores}.",
          "x" = "Missing columns: {.var count_1} or {.var count_2}."
        )
      )
    }
  }
  assert_single_value(min_cells_count, "numeric", allow_null = TRUE)
  if (!is.null(min_cells_count)) {
    assert_within_limits(min_cells_count, limits = c(0, Inf))
  }

  if (is.null(c(background_threshold_pct, background_threshold_count, min_cells_count))) {
    cli::cli_abort(
      c(
        "i" = "You must provide at least one of the following arguments: ",
        "{.arg background_threshold_pct}, {.arg background_threshold_count}, or {.arg min_cells_count}.",
        "x" = "No filters provided."
      )
    )
  }
  assert_single_value(name, "string")

  object <- object %>%
    {
      if (!is.null(background_threshold_pct)) {
        filter(., pmin(p1, p2) >= background_threshold_pct)
      } else {
        .
      }
    } %>%
    {
      if (!is.null(background_threshold_count)) {
        filter(., pmin(count_1, count_2) >= background_threshold_count)
      } else {
        .
      }
    } %>%
    {
      if (!is.null(min_cells_count)) {
        group_by(., marker_1, marker_2) %>%
          mutate(n_cells = n()) %>%
          filter(n_cells >= min_cells_count) %>%
          select(-n_cells) %>%
          ungroup()
      } else {
        .
      }
    } %>%
    # Here we can set the name of the filtered table in the database
    compute(name = name, overwrite = TRUE)

  return(object)
}

#' @method FilterProximityScores data.frame
#' @rdname FilterProximityScores
#'
#' @export
#'
FilterProximityScores.data.frame <- function(
  object,
  background_threshold_pct = NULL,
  background_threshold_count = NULL,
  min_cells_count = NULL,
  ...
) {
  object <- object %>%
    arrow::to_duckdb() %>%
    FilterProximityScores(
      background_threshold_pct = background_threshold_pct,
      background_threshold_count = background_threshold_count,
      min_cells_count = min_cells_count
    ) %>%
    collect()

  return(object)
}


#' @method SummarizeProximityScores tbl_lazy
#' @rdname SummarizeProximityScores
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' se <- ReadPNA_Seurat(pxl_file)
#' proximity_table <- ProximityScores(se)
#'
#' # Default method uses mean
#' SummarizeProximityScores(proximity_table) %>% head()
#'
#' # Switch to median
#' SummarizeProximityScores(proximity_table, summary_stat = "median") %>% head()
#'
#' # Ignore missing values
#' SummarizeProximityScores(proximity_table, include_missing_obs = FALSE) %>% head()
#'
#' # Return lists which can be used to compute custom summary statistics
#' SummarizeProximityScores(proximity_table, detailed = TRUE) %>%
#'   # It's important to do rowwise computations
#'   rowwise() %>%
#'   mutate(
#'     sd = sd(unlist(log2_ratio_list)),
#'     iqr = IQR(unlist(log2_ratio_list)),
#'     mad = mad(unlist(log2_ratio_list)),
#'     q90 = quantile(unlist(log2_ratio_list), 0.9)
#'   ) %>%
#'   select(marker_1, marker_2, sd, iqr, mad, q90) %>%
#'   ungroup()
#'
#' @export
#'
SummarizeProximityScores.tbl_lazy <- function(
  object,
  proximity_metric = c("log2_ratio", "join_count_z"),
  group_vars = NULL,
  include_missing_obs = TRUE,
  summary_stat = c("mean", "median"),
  detailed = FALSE,
  ...
) {
  proximity_metric <- match.arg(proximity_metric, c("log2_ratio", "join_count_z"))
  summary_stat <- match.arg(summary_stat, c("mean", "median"))
  assert_col_in_data(proximity_metric, object)
  assert_vector(group_vars, "character", n = 1, allow_null = TRUE)
  assert_col_in_data("marker_1", object)
  assert_col_in_data("marker_2", object)
  assert_col_in_data("component", object)
  assert_single_value(include_missing_obs, "bool")

  # Validate proximity_metric class
  proximity_metric_slice <- object %>%
    head() %>%
    pull(all_of(proximity_metric))
  if (!inherits(proximity_metric_slice, c("numeric", "integer"))) {
    cli::cli_abort(
      c(
        "i" = "The {.var proximity_metric} column must be of type {.cls {c('numeric', 'integer')}}.",
        "x" = "Column {.str {proximity_metric}} is a {.cls {class(proximity_metric_slice)}}."
      )
    )
  }

  if (!is.null(group_vars)) {
    for (group_var in group_vars) {
      assert_col_in_data(group_var, object)
      group_var_slice <- object %>%
        head() %>%
        pull(all_of(group_var))
      if (!inherits(group_var_slice, c("character", "factor"))) {
        cli::cli_abort(
          c(
            "i" = "The {.var group_vars} column(s) must be of type {.cls {c('character', 'factor')}}.",
            "x" = "Column {.str {group_var}} is a {.cls {class(group_var_slice)}}."
          )
        )
      }
    }
  }

  # Set the summary function
  prefix <- ifelse(summary_stat == "median", "median_", "mean_")
  summary_fkn <- switch(summary_stat,
    "median" = median,
    "mean" = mean
  )

  # Count number of components
  n_cells <- object %>%
    {
      if (!is.null(group_vars)) {
        group_by(., !!!syms(group_vars))
      } else {
        .
      }
    } %>%
    summarize(n_cells = as.integer(n_distinct(component)), .groups = "drop") %>%
    compute("n_cells", overwrite = TRUE) %>%
    collect()

  if (is.null(group_vars)) {
    n_cells <- n_cells %>% pull(n_cells)
  }

  # Calculate the median proximity score per group
  object <- object %>%
    {
      if (!is.null(group_vars)) {
        group_by(., !!!syms(group_vars))
      } else {
        .
      }
    } %>%
    group_by(marker_1, marker_2, .add = TRUE) %>%
    {
      if (detailed) {
        summarize(.,
          n_cells_detected = as.integer(n()),
          join_count_list = list(join_count),
          join_count_expected_mean_list = list(join_count_expected_mean),
          !!sym(paste0(proximity_metric, "_list")) := list(!!sym(proximity_metric)),
          .groups = "drop"
        )
      } else {
        summarize(.,
          n_cells_detected = as.integer(n()),
          !!sym(paste0(proximity_metric, "_list")) := list(!!sym(proximity_metric)),
          .groups = "drop"
        )
      }
    } %>%
    collect() %>%
    {
      if (!is.null(group_vars)) {
        left_join(., n_cells, by = group_vars)
      } else {
        mutate(., n_cells = n_cells)
      }
    } %>%
    rowwise() %>%
    mutate(n_cells_missing = n_cells - n_cells_detected) %>%
    mutate(pct_detected = n_cells_detected / n_cells) %>%
    {
      if (include_missing_obs && detailed) {
        mutate(
          .,
          !!sym(paste0(proximity_metric, "_list")) := list(
            c(unlist(!!sym(paste0(proximity_metric, "_list"))), rep(0, n_cells_missing))
          ),
          join_count_list = list(
            c(unlist(join_count_list), rep(0, n_cells_missing))
          ),
          join_count_expected_mean_list = list(
            c(unlist(join_count_expected_mean_list), rep(0, n_cells_missing))
          )
        )
      } else {
        mutate(
          .,
          !!sym(paste0(proximity_metric, "_list")) := list(
            c(unlist(!!sym(paste0(proximity_metric, "_list"))), rep(0, n_cells_missing))
          )
        )
      }
    } %>%
    mutate(
      .,
      !!sym(paste0(prefix, proximity_metric)) :=
        summary_fkn(unlist(!!sym(paste0(proximity_metric, "_list"))))
    ) %>%
    {
      if (!detailed) {
        select(., -all_of(paste0(proximity_metric, "_list")))
      } else {
        (
          .
        )
      }
    } %>%
    ungroup()

  return(object)
}

#' @method SummarizeProximityScores data.frame
#' @rdname SummarizeProximityScores
#'
#' @export
#'
SummarizeProximityScores.data.frame <- function(
  object,
  proximity_metric = "log2_ratio",
  group_vars = NULL,
  include_missing_obs = TRUE,
  summary_stat = c("mean", "median"),
  detailed = FALSE,
  ...
) {
  object <- object %>%
    arrow::to_duckdb() %>%
    SummarizeProximityScores(
      proximity_metric = proximity_metric,
      group_vars = group_vars,
      include_missing_obs = include_missing_obs,
      summary_stat = summary_stat,
      detailed = detailed
    ) %>%
    collect()

  return(object)
}
