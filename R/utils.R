#' Fetch file extension from file name string
#' @noRd
.file_ext <- function(x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

#' Generate a random string
#' @noRd
.generate_random_string <- function(
  n = 5
) {
  paste(sample(c(0:9, letters, LETTERS), n, replace = TRUE), collapse = "")
}

#' Tidy results from wilcoxon test
#' @noRd
.tidy <- function(
  test_result
) {
  stopifnot(inherits(test_result, what = "htest"))
  tibble(
    estimate = test_result$estimate,
    statistic = test_result$statistic,
    p.value = test_result$p.value,
    conf.low = test_result$conf.int[1],
    conf.high = test_result$conf.int[2],
    method = test_result$method,
    alternative = test_result$alternative
  )
}

#' Create a temporary directory
#'
#' This function will attempt to create a clean directory in
#' temp dir with the given name \code{dir_name}. If the directory
#' already exists, it will first be deleted. If the deletion fails,
#' a random string will be appended to the directory name to make
#' sure that the directory is unique.
#'
#' @param dir_name A character string with the directory name
#'
#' @noRd
.create_unique_temp_dir <- function(dir_name) {
  temp_layout_dir <- file.path(fs::path_temp(), dir_name)
  while (fs::dir_exists(temp_layout_dir)) {
    err <- try(fs::dir_delete(temp_layout_dir))
    if (inherits(err, "try-error")) {
      warn(glue("Failed to delete temporary directory '{col_br_blue(temp_layout_dir)}'."))
      temp_layout_dir <- file.path(fs::path_temp(), dir_name, .generate_random_string())
    }
  }
  fs::dir_create(temp_layout_dir)
  return(temp_layout_dir)
}


#' Check if a path is absolute
#' @noRd
.is_absolute_path <- function(
  x
) {
  if (.Platform$OS.type == "unix") {
    str_detect(x, "^[/~]")
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
.validate_polarization <- function(
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
    name_check <- all(c("marker", "component") %in% names(polarization))
    if (!name_check) {
      abort("Columns 'marker', and 'component' are required in the 'polarization' score table")
    }
    # Check component names
    cells_in_polarization <- cell_ids %in% (polarization$component %>% unique())
    if (!all(cells_in_polarization)) {
      if (verbose && check_global_verbosity()) {
        cli_alert_warning(glue(
          "Cells {paste(which(!cells_in_polarization), collapse=', ')} ",
          "in 'counts' are missing from 'polarization' table"
        ))
      }
      return(polarization)
    }
    # Check marker names
    markers_in_polarization <- markers %in% (polarization$marker %>% unique())
    if (!all(markers_in_polarization)) {
      if (verbose && check_global_verbosity()) {
        cli_alert_warning(glue(
          "Markers {paste(which(!markers_in_polarization), collapse=', ')} ",
          "in 'counts' are missing from 'polarization' table"
        ))
      }
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
.validate_colocalization <- function(
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
    name_check <- all(c("marker_1", "marker_2", "component") %in% names(colocalization))
    if (!name_check) {
      abort("Columns 'marker_1', 'marker_2', and 'component' are required in the 'colocalization' score table")
    }
    # Check component names
    cells_in_colocalization <- cell_ids %in% (colocalization$component %>% unique())
    if (!all(cells_in_colocalization)) {
      if (verbose && check_global_verbosity()) {
        cli_alert_warning(glue(
          "Cells {paste(which(!cells_in_colocalization), collapse=', ')} ",
          "in 'counts' are missing from 'colocalization' table"
        ))
      }
      return(colocalization)
    }
    # Check marker names
    all_markers <- c(colocalization$marker_1 %>% unique(), colocalization$marker_2 %>% unique()) %>% unique()
    if (!all(markers %in% all_markers)) {
      if (verbose && check_global_verbosity()) {
        cli_alert_warning(glue(
          "Markers {paste(which(!markers %in% all_markers), collapse=', ')} ",
          "in 'counts' are missing from 'colocalization' table"
        ))
      }
    }
  }
  return(colocalization)
}


#' Validate an fs_map tibble
#'
#' @noRd
.validate_fs_map <- function(
  fs_map
) {
  if (!inherits(fs_map, what = "tbl_df")) {
    abort("'fs_map' must be a 'tbl_df'")
  }
  if (length(fs_map) == 0) {
    abort("'fs_map' must be a non-empty 'tbl_df'")
  }
  if (!all(c("id_map", "sample", "pxl_file") == colnames(fs_map))) {
    abort("'fs_map' must have columns 'id_map', 'sample', and 'pxl_file'")
  }

  # Validate columns
  fs_map_classes <- sapply(fs_map, class)
  if (fs_map_classes["id_map"] != "list") {
    abort(glue("Column 'id_map' of `fs_map` must be a 'list'. Got '{fs_map_classes['id_map']}'"))
  }
  if (fs_map_classes["sample"] != "integer") {
    abort(glue("Column 'sample' of `fs_map` must be an 'integer' vector. Got '{fs_map_classes['sample']}'"))
  }
  if (fs_map_classes["pxl_file"] != "character") {
    abort(glue("Column 'sample' of `fs_map` must be an 'character' vector. Got '{fs_map_classes['pxl_file']}'"))
  }
  id_map_check1 <- sapply(fs_map$id_map, function(x) inherits(x, what = "tbl_df"))
  if (!all(id_map_check1)) {
    abort("All elements of 'id_map' must be 'tbl_df' objects")
  }
  id_map_check2 <- sapply(fs_map$id_map, function(x) all(colnames(x) == c("current_id", "original_id")))
  if (!all(id_map_check2)) {
    abort("All elements of 'id_map' must be 'tbl_df' objects with columns 'current_id' and 'original_id'")
  }
  id_map_check3 <- sapply(fs_map$id_map, function(x) {
    all(sapply(x, class) == c("character", "character"))
  })
  if (!all(id_map_check3)) {
    abort("Invalid classes found in 'id_map' columns 'current_id' and 'original_id'")
  }

  # Validate pxl files
  for (i in seq_len(nrow(fs_map))) {
    f <- fs_map$pxl_file[i]
    if (!fs::file_exists(f)) {
      abort(glue(
        "The pxl file '{col_br_blue(f)}' linked to sample {i} does not exist. ",
        "Make sure that the path is correct and that the file has not been moved/deleted."
      ))
    }
    # Check .pxl file for content
    pxl_files <- unzip(f, list = TRUE)$Name


    required_files <- c("adata.h5ad", "edgelist.parquet", "metadata.json")
    expected_files <- c(
      "adata.h5ad", "colocalization.parquet", "edgelist.parquet",
      "metadata.json", "polarization.parquet"
    )

    if (!all(expected_files %in% pxl_files)) {
      warn(glue("The pxl file '{col_br_blue(f)}' is missing: \n
                {paste(setdiff(expected_files, pxl_files), collapse = '\n')}"))
    }
    if (!all(required_files %in% pxl_files)) {
      abort(glue("The pxl file '{col_br_blue(f)}' is invalid. "))
    }
  }

  return(invisible(NULL))
}

#' Utility function to capture warnings and errors
#'
#' @param expr Expression to evaluate
#'
#' @noRd
evaluate_with_catch <- function(expr) {
  result <- NULL
  warning_message <- NULL
  error_message <- NULL

  tryCatch(
    {
      # Capture warnings
      withCallingHandlers(
        {
          result <- eval(expr)
        },
        warning = function(w) {
          warning_message <<- w$message
          invokeRestart("muffleWarning")
        }
      )
    },
    error = function(e) {
      error_message <<- e$message
    }
  )

  list(
    result = result,
    error = error_message,
    warning = warning_message
  )
}


#' Utility function to capture warnings and errors
#'
#' @param ... Any number of R expressions, which should each evaluate
#' to (a logical vector of all) TRUE
#' @param message Message to display if the condition is not met
#' @param class Class of the error to throw
#' @param not Logical indicating whether to negate the condition
#' @param call Environment to use for the error call
#'
#' @noRd
abort_if_not <- function(
  ...,
  message = NULL,
  class = NULL,
  not = TRUE,
  call = caller_env()
) {
  params <- list2(...)

  if (length(params) > 1) {
    for (x in seq(params)) {
      abort_if_not(
        params[[x]],
        message = names(params)[[x]],
        class = class,
        call = call
      )
    }

    invisible()
  }

  condition <- params[[1]]

  if (is.logical(condition) && ((!condition && not) || (condition && !not))) {
    if (is_named(params)) {
      message <- message %||% names(params)
    }

    abort(
      message = glue(message, .envir = call),
      class = class,
      call = call
    )
  }

  invisible()
}


#' @param object A tibble
#' @param contrast_column A character vector of length 1
#' @param reference A character vector of length 1
#' @param targets A character vector of length >= 1 or NULL
#' @param group_vars A character vector of length >= 1 or NULL
#' @param spatial_metric A character vector of length 1
#' @param min_n_obs An integer of length 1
#' @param conf_int Either TRUE or FALSE
#' @param cl An integer, a cluster object or NULL
#' @param data_type A character vector of length 1
#'
#' @noRd
.validate_dpa_dca_input <- function(
  object,
  contrast_column,
  reference,
  targets,
  group_vars,
  spatial_metric,
  min_n_obs,
  conf_int,
  cl,
  data_type
) {
  metric_name <- paste0(data_type, "_metric")

  # Validate contrast column
  abort_if_not(
    "contrast_column = '{contrast_column}' is invalid
    contrast_column must be a valid column name" =
      inherits(contrast_column, what = "character") &&
        (length(contrast_column) == 1) &&
        (contrast_column %in% colnames(object))
  )
  group_vector <- object[, contrast_column, drop = TRUE]
  abort_if_not(
    "contrast_column = '{contrast_column}' is invalid
    contrast_column must be a character vector or a factor" =
      inherits(group_vector, what = c("character", "factor")),
    "contrast_column = '{contrast_column}' is invalid
    contrast_column must have at least 2 groups" =
      length(unique(group_vector)) > 1
  )

  # Validate reference
  abort_if_not(
    "reference = '{reference}' is invalid
    reference must be present in the '{contrast_column}' column" =
      inherits(reference, what = "character") &&
        (length(reference) == 1) &&
        (reference %in% group_vector)
  )

  # Validate targets
  if (!is.null(targets)) {
    abort_if_not(
      "targets must be a character vector" =
        inherits(targets, what = "character")
    )

    abort_if_not(
      "targets is invalid
      all targets must be different from reference = '{reference}'" =
        !(reference %in% targets)
    )
    for (target in targets) {
      abort_if_not(
        "'{target}' in targets is invalid
        '{target}' must be present in the '{contrast_column}' column" =
          target %in% group_vector
      )
    }
  }

  # Check for component and marker or marker_1/marker_2 columns
  if (data_type == "polarity") {
    abort_if_not(
      "'component' and 'marker' must be present in the {data_type} score table" =
        all(c("marker", "component") %in% colnames(object))
    )
  }
  if (data_type == "colocalization") {
    abort_if_not(
      "'component', 'marker_1' and 'marker_2' must be present in the {data_type} score table" =
        all(c("marker_1", "marker_2", "component") %in% colnames(object))
    )
  }

  # Validate spatial metric
  abort_if_not(
    "{metric_name} = '{spatial_metric}' is missing from the {data_type} score table." =
      spatial_metric %in% colnames(object),
    "conf_int = '{conf_int}'
    must be TRUE or FALSE" =
      inherits(conf_int, what = "logical") & (length(conf_int) == 1)
  )
  if (!inherits(object[, spatial_metric, drop = TRUE], what = "numeric")) {
    abort(
      glue(
        "Column '{spatial_metric}' (polarity_metric) is a '{class(object[, spatial_metric, drop = TRUE])}' ",
        "vector but must be a 'numeric' vector.\n"
      )
    )
  }

  # Validate group_vars
  if (!is.null(group_vars)) {
    abort_if_not(
      "group_vars is invalid
      group_vars must be a character vector with valid column names" =
        inherits(group_vars, what = "character") &&
          (length(group_vars) >= 1) &&
          all(group_vars %in% colnames(object))
    )
    abort_if_not(
      "contrast_column = '{contrast_column}' cannot be one of group_vars" =
        !contrast_column %in% group_vars
    )
    for (group_var in group_vars) {
      abort_if_not(
        "group variable '{group_var}' must be a character vector or a factor" =
          inherits(object[, group_var, drop = TRUE],
            what = c("character", "factor")
          )
      )
    }
  }

  # Validate min_n_obs
  abort_if_not(
    "min_n_obs = '{min_n_obs}' is invalid
    min_n_obs must be an integer" =
      inherits(min_n_obs, what = "numeric") &&
        (length(min_n_obs) == 1) &&
        (min_n_obs >= 0)
  )

  # Validate cl
  if (!is.null(cl)) {
    abort_if_not(
      "'cl' must be a cluster object or an integer" =
        inherits(cl, what = c("cluster", "numeric"))
    )
  }
}


#' @param object A tibble
#' @param contrast_column A character vector of length 1
#' @param reference A character vector of length 1
#' @param targets A character vector of length >= 1 or NULL
#' @param group_vars A character vector of length >= 1 or NULL
#'
#' @noRd
.validate_daa_input <- function(
  object,
  contrast_column,
  reference,
  targets,
  group_vars
) {

  # Validate contrast column
  abort_if_not(
    "contrast_column = '{contrast_column}' is invalid
    contrast_column must be a valid column name" =
      inherits(contrast_column, what = "character") &&
      (length(contrast_column) == 1) &&
      (contrast_column %in% colnames(object))
  )
  group_vector <- object[, contrast_column, drop = TRUE]
  abort_if_not(
    "contrast_column = '{contrast_column}' is invalid
    contrast_column must be a character vector or a factor" =
      inherits(group_vector, what = c("character", "factor")),
    "contrast_column = '{contrast_column}' is invalid
    contrast_column must have at least 2 groups" =
      length(unique(group_vector)) > 1
  )

  # Validate reference
  abort_if_not(
    "reference = '{reference}' is invalid
    reference must be present in the '{contrast_column}' column" =
      inherits(reference, what = "character") &&
      (length(reference) == 1) &&
      (reference %in% group_vector)
  )

  # Validate targets
  if (!is.null(targets)) {
    abort_if_not(
      "targets must be a character vector" =
        inherits(targets, what = "character")
    )

    abort_if_not(
      "targets is invalid
      all targets must be different from reference = '{reference}'" =
        !(reference %in% targets)
    )
    for (target in targets) {
      abort_if_not(
        "'{target}' in targets is invalid
        '{target}' must be present in the '{contrast_column}' column" =
          target %in% group_vector
      )
    }
  }

  # Validate group_vars
  if (!is.null(group_vars)) {
    abort_if_not(
      "group_vars is invalid
      group_vars must be a character vector with valid column names" =
        inherits(group_vars, what = "character") &&
        (length(group_vars) >= 1) &&
        all(group_vars %in% colnames(object))
    )
    abort_if_not(
      "contrast_column = '{contrast_column}' cannot be one of group_vars" =
        !contrast_column %in% group_vars
    )
    for (group_var in group_vars) {
      abort_if_not(
        "group variable '{group_var}' must be a character vector or a factor" =
          inherits(object[, group_var, drop = TRUE],
                   what = c("character", "factor")
          )
      )
    }
  }
}
