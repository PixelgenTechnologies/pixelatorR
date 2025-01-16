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
  assert_class(test_result, "htest")
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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# VALIDATION FUNCTIONS
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
#' @param call Environment to use for the error call
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @noRd
.validate_polarization <- function(
  polarization,
  cell_ids,
  markers,
  call = caller_env(),
  verbose = FALSE
) {
  # Set polarization to empty tibble if NULL
  polarization <- polarization %||% tibble()
  assert_class(polarization, "tbl_df", call = call)

  # Validate polarization and colocalization
  if (length(polarization) > 0) {
    # Check column names
    name_check <- all(c("marker", "component") %in% names(polarization))
    if (!name_check) {
      cli::cli_abort(
        c(
          "i" = "Columns {.str marker} and {.str component} must be present in the 'polarization' score table"
        ),
        call = call
      )
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
#' @param call Environment to use for the error call
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @noRd
.validate_colocalization <- function(
  colocalization,
  cell_ids,
  markers,
  call = caller_env(),
  verbose = FALSE
) {
  # Set colocalization to empty tibble if NULL
  colocalization <- colocalization %||% tibble()
  assert_class(colocalization, "tbl_df", call = call)

  if (length(colocalization) > 0) {
    # Check column names
    name_check <- all(c("marker_1", "marker_2", "component") %in% names(colocalization))
    if (!name_check) {
      cli::cli_abort(
        c(
          "i" = "Columns {.str marker_1}, {.str marker_2} and {.str component} must",
          " " = "be present in the 'colocalization' score table"
        ),
        call = call
      )
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
  fs_map,
  allow_empty_df = FALSE,
  call = caller_env()
) {
  if (inherits(fs_map, "data.frame") && length(fs_map) == 0 && allow_empty_df) {
    return(invisible(NULL))
  }
  assert_non_empty_object(fs_map, "tbl_df", call = call)
  if (!all(c("id_map", "sample", "pxl_file") == colnames(fs_map))) {
    cli::cli_abort(
      c(
        "i" = "{.var fs_map} must have columns {.str id_map}, {.str sample}, and {.str pxl_file}",
        "x" = "Found column(s) {.val {colnames(fs_map)}} in {.var fs_map}"
      ),
      call = call
    )
  }

  # Validate columns
  assert_col_class("id_map", fs_map, classes = "list", call = call)
  assert_col_class("sample", fs_map, classes = "integer", call = call)
  assert_col_class("pxl_file", fs_map, classes = "character", call = call)
  id_map_check1 <- sapply(fs_map$id_map, function(x) inherits(x, what = "tbl_df"))
  if (!all(id_map_check1)) {
    cli::cli_abort(
      c("x" = "All elements of {.var id_map} must be {.cls tbl_df} objects"),
      call = call
    )
  }
  id_map_check2 <- sapply(fs_map$id_map, function(x) all(colnames(x) == c("current_id", "original_id")))
  if (!all(id_map_check2)) {
    cli::cli_abort(
      c("x" = "All elements of {.var id_map} must be {.cls tbl_df} objects with",
        " " = "columns {.str current_id} and {.str original_id}"),
      call = call
    )
  }
  id_map_check3 <- sapply(fs_map$id_map, function(x) {
    all(sapply(x, class) == c("character", "character"))
  })
  if (!all(id_map_check3)) {
    cli::cli_abort(
      c("x" = "Invalid classes found in {.var id_map} columns 'current_id' and 'original_id'"),
      call = call
    )
  }

  # Validate pxl files
  for (i in seq_len(nrow(fs_map))) {
    f <- fs_map$pxl_file[i]
    if (!fs::file_exists(f)) {
      cli::cli_abort(
        c(
          "x" = "The pxl file {.file {f}} linked to sample {.val {i}} does not exist. ",
          "i" = "Make sure that the path is correct and that the file has not been moved/deleted."
        ),
        call = call
      )
    }
    # Check .pxl file for content
    pxl_files <- utils::unzip(f, list = TRUE)$Name


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
      cli::cli_abort(
        c("x" = "The pxl file {.file {f}} is missing required data."),
        call = call
      )
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
#' @param call Environment to use for the error call
#'
#' @noRd
.validate_da_input <- function(
  object,
  contrast_column,
  reference,
  targets,
  group_vars,
  spatial_metric = NULL,
  min_n_obs = 0,
  conf_int = FALSE,
  cl = NULL,
  data_type = NULL,
  call = caller_env()
) {
  metric_name <- data_type %||% paste0(data_type, "_metric")

  # Validate contrast column
  assert_single_value(contrast_column, type = "string", call = call)
  assert_col_in_data(contrast_column, object, call = call)
  assert_col_class(contrast_column, object, classes = c("character", "factor"), call = call)
  group_vector <- object[, contrast_column, drop = TRUE]
  if (!length(unique(group_vector)) > 1) {
    cli::cli_abort(
      c(
        "i" = "Group variable {.val {contrast_column}} must have at least 2 groups",
        "x" = "Group variable {.val {contrast_column}} only has 1 unique group"
      ),
      call = call
    )
  }

  # Validate reference
  assert_single_value(reference, type = "string", call = call)
  if (!(reference %in% group_vector)) {
    cli::cli_abort(
      c("x" = "Reference group {.val {reference}} must be present in the {.val {contrast_column}} column"),
      call = call
    )
  }

  # Validate targets
  if (!is.null(targets)) {
    assert_vector(targets, type = "character", n = 1, call = call)
    if (reference %in% targets) {
      cli::cli_abort(
        c("x" = "All {.var targets} ({.val {targets}}) must be different from {.var reference} ({.val {reference}})"),
        call = call
      )
    }

    for (target in targets) {
      if (!target %in% group_vector) {
        cli::cli_abort(
          c("x" = "Target group {.val {target}} must be present in the {.val {contrast_column}} column"),
          call = call
        )
      }
    }
  }

  # Validate spatial metrics if available
  if (!is.null(data_type)) {
    # Check for component and marker or marker_1/marker_2 columns
    if (data_type == "polarity") {
      assert_col_in_data("marker", object, call = call)
      assert_col_in_data("component", object, call = call)
    }
    if (data_type == "colocalization") {
      assert_col_in_data("marker_1", object, call = call)
      assert_col_in_data("marker_2", object, call = call)
      assert_col_in_data("component", object, call = call)
    }

    # Validate spatial metric
    assert_col_in_data(spatial_metric, object, call = call)
    assert_single_value(conf_int, type = "bool", call = call)
    assert_col_class(spatial_metric, object, classes = "numeric", call = call)
  }

  # Validate group_vars
  if (!is.null(group_vars)) {
    assert_vector(group_vars, type = "character", n = 1, call = call)
    assert_x_in_y(group_vars, colnames(object), call = call)
    if (contrast_column %in% group_vars) {
      cli::cli_abort(
        c("x" = "Group variable {.val {contrast_column}} cannot be one of {.var group_vars}"),
        call = call
      )
    }
    for (group_var in group_vars) {
      assert_col_class(group_var, object, classes = c("character", "factor"), call = call)
    }
  }

  # Validate min_n_obs
  assert_single_value(min_n_obs, type = "numeric", call = call)
  if (min_n_obs < 0) {
    cli::cli_abort(
      c("x" = "{.var min_n_obs} must be a non-negative integer"),
      call = call
    )
  }

  # Validate cl
  assert_class(cl, classes = c("numeric", "cluster"), allow_null = TRUE, call = call)
}


#' Validate assay argument
#'
#' @return assay
#'
#' @noRd
.validate_or_set_assay <- function(object, assay = NULL, call = caller_env()) {
  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    assert_single_value(assay, type = "character", call = call)
  } else {
    assay <- DefaultAssay(object)
  }
  return(assay)
}
