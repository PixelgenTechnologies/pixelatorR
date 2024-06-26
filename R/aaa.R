# Declarations used in package check
globalVariables(
  names = c('id_map', 'component_new', 'tau', 'tau_type', 'umi_per_upia',
            'upia1', 'upia2', 'component', 'rn', 'x', 'y', 'z', 'name',
            'type', 'g', 'from', 'to', 'node_type', 'id', 'layout',
            'pearson_z', 'p', 'p.value', '.', 'original_id', 'current_id',
            'graph_projection', "label", "in_gate", "dens", "xmax",
            "xmin", "ymax", "ymin"),
  package = 'pixelatorR',
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
  if (is.null(getOption("pixelatorR.verbose"))) return(TRUE)
  getOption("pixelatorR.verbose")
}


# ***********************************
# Check if certain packages are installed. When calling these functions,
# if the package is missing, users will be asked to install the package.
# ***********************************

expect_scales <- function(...) {
  rlang::check_installed('scales', ...)
}

expect_graphlayouts <- function(...) {
  rlang::check_installed('graphlayouts', ...)
}

expect_pheatmap <- function(...) {
  rlang::check_installed('pheatmap', ...)
}

expect_duckdb <- function(...) {
  rlang::check_installed('duckdb', ...)
}

expect_MASS <- function(...) {
  rlang::check_installed('MASS', ...)
}

