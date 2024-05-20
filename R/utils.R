#' Fetch file extension from file name string
#' @noRd
.file_ext <- function (x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

#' Generate a random string
#' @noRd
.generate_random_string <- function (
  n = 5
) {
  paste(sample(c(0:9, letters, LETTERS), n, replace = TRUE), collapse = "")
}

#' Tidy results from wilcoxon test
#' @noRd
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


#' Check if a path is absolute
#' @noRd
.is_absolute_path <- function (
  x
) {
  if (.Platform$OS.type == 'unix') {
    str_detect(x, '^[/~]')
  } else {
    str_detect(x, "^(~|.:)(/|\\\\)")
  }
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
  stopifnot(
    "'polarization' must be a non-empty 'tbl_df' object" =
      inherits(polarization, "tbl_df")
  )

  # Validate polarization and colocalization
  if (length(polarization) > 0) {
    # Check column names
    stopifnot(
      "'polarization' names are invalid" =
        all(sort(names(polarization)) ==
              sort(c("morans_i", "morans_p_value", "morans_p_adjusted",
                     "morans_z", "marker", "component")))
    )
    # Check component names
    cells_in_polarization <- cell_ids %in% (polarization$component %>% unique())
    if (!all(cells_in_polarization)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning(glue("Cells {paste(which(!cells_in_polarization), collapse=', ')} ",
                               "in 'counts' are missing from 'polarization' table"))
      return(polarization)
    }
    # Check marker names
    markers_in_polarization <- markers %in% (polarization$marker %>% unique())
    if (!all(markers_in_polarization)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning(glue("Markers {paste(which(!markers_in_polarization), collapse=', ')} ",
                               "in 'counts' are missing from 'polarization' table"))
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
                      c("marker_1", "marker_2", "pearson", "pearson_mean",
                        "pearson_stdev", "pearson_z", "pearson_p_value",
                        "pearson_p_value_adjusted", "jaccard", "jaccard_mean",
                        "jaccard_stdev", "jaccard_z", "jaccard_p_value",
                        "jaccard_p_value_adjusted", "component")))
    # Check component names
    cells_in_colocalization <- cell_ids %in% (colocalization$component %>% unique())
    if (!all(cells_in_colocalization)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning(glue("Cells {paste(which(!cells_in_colocalization), collapse=', ')} ",
                               "in 'counts' are missing from 'colocalization' table"))
      return(colocalization)
    }
    # Check marker names
    all_markers <- c(colocalization$marker_1 %>% unique(), colocalization$marker_2 %>% unique()) %>% unique()
    if (!all(markers %in% all_markers)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning(glue("Markers {paste(which(!markers %in% all_markers), collapse=', ')} ",
                               "in 'counts' are missing from 'colocalization' table"))
    }
  }
  return(colocalization)
}


#' Validate an fs_map tibble
#'
#' @noRd
.validate_fs_map <- function (
  fs_map
) {

  if (!inherits(fs_map, what = "tbl_df"))
    abort("'fs_map' must be a 'tbl_df'")
  if (length(fs_map) == 0)
    abort("'fs_map' must be a non-empty 'tbl_df'")
  if (!all(c("id_map", "sample", "pxl_file") == colnames(fs_map)))
    abort("'fs_map' must have columns 'id_map', 'sample', and 'pxl_file'")

  # Validate columns
  fs_map_classes <- sapply(fs_map, class)
  if (fs_map_classes["id_map"] != "list")
    abort(glue("Column 'id_map' of `fs_map` must be a 'list'. Got '{fs_map_classes['id_map']}'"))
  if (fs_map_classes["sample"] != "integer")
    abort(glue("Column 'sample' of `fs_map` must be an 'integer' vector. Got '{fs_map_classes['sample']}'"))
  if (fs_map_classes["pxl_file"] != "character")
    abort(glue("Column 'sample' of `fs_map` must be an 'character' vector. Got '{fs_map_classes['pxl_file']}'"))
  id_map_check1 <- sapply(fs_map$id_map, function(x) inherits(x, what = "tbl_df"))
  if (!all(id_map_check1))
    abort("All elements of 'id_map' must be 'tbl_df' objects")
  id_map_check2 <- sapply(fs_map$id_map, function(x) all(colnames(x) == c("current_id", "original_id")))
  if (!all(id_map_check2))
    abort("All elements of 'id_map' must be 'tbl_df' objects with columns 'current_id' and 'original_id'")
  id_map_check3 <- sapply(fs_map$id_map, function(x) {
    all(sapply(x, class) == c("character", "character"))
  })
  if (!all(id_map_check3))
    abort("Invalid classes found in 'id_map' columns 'current_id' and 'original_id'")

  # Validate pxl files
  for (i in seq_len(nrow(fs_map))) {
    f <- fs_map$pxl_file[i]
    if (!fs::file_exists(f)) {
      abort(glue("The pxl file '{col_br_blue(f)}' linked to sample {i} does not exist. ",
                 "Make sure that the path is correct and that the file has not been moved/deleted."))
    }
    # Check .pxl file for content
    pxl_files <- unzip(f, list = TRUE)$Name
    if (!all(c("adata.h5ad", "colocalization.parquet", "edgelist.parquet",
               "metadata.json", "polarization.parquet") %in% pxl_files)) {
      abort(glue("The pxl file '{col_br_blue(f)}' is invalid. "))
    }
  }

  return(invisible(NULL))
}

