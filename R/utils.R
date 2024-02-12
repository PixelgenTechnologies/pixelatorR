.file_ext <- function (x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

.generate_random_string <- function (
  n = 5
) {
  paste(sample(c(0:9, letters, LETTERS), n, replace = TRUE), collapse = "")
}

.tidy <- function (
  test_result
) {
  stopifnot(inherits(test_result, what = "htest"))
  tibble(estimate = test_result$estimate,
         statistic = test_result$statistic,
         p.value = test_result$p.value,
         conf.low = test_result$conf.int[1],
         conf.high = test_result$conf.int[2],
         method = test_result$method,
         alternative = test_result$alternative)
}


#' Validate polarization \code{tbl_df}
#'
#' @param polarization A \code{tbl_df} with polarization scores
#' @param cell_ids A character vector with cell IDs
#' @param markers A character vector with marker names
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @noRd
.validate_polarization <- function (
    polarization,
    cell_ids,
    markers,
    verbose = FALSE
) {
  # Set polarization to empty tibble if NULL
  polarization <- polarization %||% tibble()
  stopifnot("'polarization' must be a non-empty 'tbl_df' object" = inherits(polarization, "tbl_df"))

  # Validate polarization and colocalization
  if (length(polarization) > 0) {
    # Check column names
    stopifnot("'polarization' names are invalid" =
                all(names(polarization) ==
                      c("morans_i","morans_p_value","morans_p_adjusted","morans_z","marker","component")))
    # Check component names
    cells_in_polarization <- cell_ids %in% (polarization$component %>% unique())
    if (!all(cells_in_polarization)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning("Cells {paste(which(!cells_in_polarization), collapse=', ')} in 'counts' are missing from 'polarization' table")
      return(polarization)
    }
    # Check marker names
    markers_in_polarization <- markers %in% (polarization$marker %>% unique())
    if (!all(markers_in_polarization)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning("Markers {paste(which(!markers_in_polarization), collapse=', ')} in 'counts' are missing from 'polarization' table")
    }
  }
  return(polarization)
}

#' Validate colocalization \code{tbl_df}
#'
#' @param colocalization A \code{tbl_df} with colocalization scores
#' @param cell_ids A character vector with cell IDs
#' @param markers A character vector with marker names
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @noRd
.validate_colocalization <- function (
    colocalization,
    cell_ids,
    markers,
    verbose = FALSE
) {

  # Set colocalization to empty tibble if NULL
  colocalization <- colocalization %||% tibble()
  stopifnot("'colocalization' must be a non-empty 'tbl_df' object" = inherits(colocalization, "tbl_df"))

  if (length(colocalization) > 0) {
    # Check column names
    stopifnot("'colocalization' names are invalid" =
                all(names(colocalization) ==
                      c("marker_1","marker_2","pearson","pearson_mean","pearson_stdev","pearson_z","pearson_p_value",
                        "pearson_p_value_adjusted","jaccard","jaccard_mean","jaccard_stdev","jaccard_z","jaccard_p_value",
                        "jaccard_p_value_adjusted","component")))
    # Check component names
    cells_in_colocalization <- cell_ids %in% (colocalization$component %>% unique())
    if (!all(cells_in_colocalization)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning("Cells {paste(which(!cells_in_colocalization), collapse=', ')} in 'counts' are missing from 'colocalization' table")
      return(colocalization)
    }
    # Check marker names
    all_markers <- c(colocalization$marker_1 %>% unique(), colocalization$marker_2 %>% unique()) %>% unique()
    if (!all(markers %in% all_markers)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning("Markers {paste(which(!markers %in% all_markers), collapse=', ')} in 'counts' are missing from 'colocalization' table")
    }
  }
  return(colocalization)
}
