# Declarations used in package check
globalVariables(
<<<<<<< HEAD
  names = c(
    "id_map", "component_new", "tau", "tau_type", "umi_per_upia",
    "upia1", "upia2", "component", "rn", "x", "y", "z", "name",
    "type", "g", "from", "to", "node_type", "id", "layout",
    "pearson_z", "p", "p.value", ".", "original_id", "current_id",
    "graph_projection", "label", "in_gate", "dens", "xmax",
    "xmin", "ymax", "ymin", "marker_x", "marker_y", "val", ".x",
    "value_x", "value_y", "bi_prob", "pxl_file", "value",
    "marker_1", "marker_2", "graph_projection", "modality",
    "mixture_component", "morans_z", "upia", "upib", "marker",
    "n", "norm_factor", "nodes", "group", "molecules"
  ),
  package = "pixelatorR",
=======
  names = c('id_map', 'component_new', 'tau', 'tau_type', 'umi_per_upia',
            'upia1', 'upia2', 'component', 'rn', 'x', 'y', 'z', 'name',
            'type', 'g', 'from', 'to', 'node_type', 'id', 'layout',
            'pearson_z', 'p', 'p.value', '.', 'original_id', 'current_id',
            'graph_projection', "label", "in_gate", "dens", "xmax",
            "xmin", "ymax", "ymin", "marker_x", "marker_y",
            "value_x", "value_y", "bi_prob", "pxl_file", "value", "marker_1", "marker_1",
            'graph_projection', "modality", "mixture_component", "frequency"),
  package = 'pixelatorR',
>>>>>>> 03c6361 (fix documentation)
  add = TRUE
)

#' Check global option for verbosity
#'
#' By setting the global option \code{options(pixelatorR.verbose = FALSE)},
#' users can turn off all verbosity in pixelatorR. This function can be used to
#' check the value of this global options
#'
#' @noRd
check_global_verbosity <- function() {
  if (is.null(getOption("pixelatorR.verbose"))) {
    return(TRUE)
  }
  getOption("pixelatorR.verbose")
}


# ***********************************
# Check if certain packages are installed. When calling these functions,
# if the package is missing, users will be asked to install the package.
# ***********************************

expect_scales <- function(...) {
  rlang::check_installed("scales", ...)
}

expect_graphlayouts <- function(...) {
  rlang::check_installed("graphlayouts", ...)
}

expect_pheatmap <- function(...) {
  rlang::check_installed("pheatmap", ...)
}

expect_duckdb <- function(...) {
  rlang::check_installed("duckdb", ...)
}

expect_MASS <- function(...) {
  rlang::check_installed("MASS", ...)
}

expect_mclust <- function(...) {
  rlang::check_installed("mclust", ...)
}

expect_limma <- function(...) {
  rlang::check_installed("limma", ...)
}

expect_styler <- function(...) {
  rlang::check_installed("styler", ...)
}
