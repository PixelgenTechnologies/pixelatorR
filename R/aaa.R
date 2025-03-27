# Declarations used in package check
globalVariables(
  names = c(
    ".", ".x", "bi_prob", "component", "component_new", "current_id",
    "dens", "frequency", "from", "g", "graph_projection", "group",
    "hjust", "id", "id_map", "in_gate", "join_count", "join_count_expected_mean",
    "label", "layout", "marker", "marker_1", "marker_2", "marker_x",
    "marker_y", "marker1", "marker2", "med_ref", "med_tgt", "mixture_component",
    "modality", "molecules", "morans_z", "n", "n_inside", "n_ref",
    "n_tgt", "name", "node_type", "nodes", "norm_factor", "nties",
    "nties_const", "original_id", "p", "p_adj", "p_val", "p_val_adj",
    "p.value", "pct.1", "pct.2", "pearson_z", "pxl_file", "quadrant",
    "r", "ref_n", "rn", "rs_ref", "rs_tgt", "sigma", "target_n",
    "tau", "tau_type", "to", "total", "type", "u", "umi_per_upia",
    "umi1", "umi2", "upia", "upia1", "upia2", "upib", "val", "value",
    "value_x", "value_y", "vjust", "x", "x_label", "xmax", "xmin",
    "y", "y_label", "ymax", "ymin", "z"
  ),
  package = "pixelatorR",
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

expect_ComplexHeatmap <- function(...) {
  rlang::check_installed("ComplexHeatmap", ...)
}

expect_Seurat <- function(...) {
  rlang::check_installed("Seurat", ...)
}

expect_zip <- function(...) {
  rlang::check_installed("zip", ...)
}

expect_dtplyr <- function(...) {
  rlang::check_installed("dtplyr", ...)
}
