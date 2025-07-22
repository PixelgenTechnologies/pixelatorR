# Declarations used in package check
globalVariables(
  names = c(".", ".x", "bi_prob", "color", "community", "comp", "component",
            "component_new", "count_1", "count_2", "current_id", "dens",
            "dev_png", "doublet_nn_rate", "doublet_nns", "doublet_p", "doublet_p_adj",
            "edge", "frequency", "from", "g", "graph_edges", "graph_projection",
            "graph_proteins", "graph_reads", "group", "hjust", "hp", "hue",
            "hup", "i1", "i2", "id", "id_map", "in_gate", "index", "item",
            "join_count", "join_count_expected_mean", "join_count_expected_mean_list",
            "join_count_list", "label", "layout", "level", "marker", "marker_1",
            "marker_2", "marker_x", "marker_y", "marker1", "marker2", "med_ref",
            "med_tgt", "mixture_component", "modality", "molecules", "morans_z",
            "n", "n_cells", "n_cells_detected", "n_cells_missing", "n_inside",
            "n_ref", "n_tgt", "name", "neighbor", "node_type", "node_val",
            "nodes", "norm_factor", "nties", "nties_const", "original_id",
            "p", "p_adj", "p_val", "p_val_adj", "p.value", "p1", "p2", "patch",
            "pct.1", "pct.2", "pearson_z", "png", "pxl_file", "quadrant",
            "r", "read_count", "receiver_freq", "receiver_unmixed_freq",
            "ref_n", "rn", "rs_ref", "rs_tgt", "sigma", "simulated", "size",
            "target_freq", "target_n", "target_unmixed_freq", "tau", "tau_type",
            "text_color", "to", "total", "tp", "trial", "tup", "type", "u",
            "uei_count", "umi_count", "umi_per_upia", "umi1", "umi2", "upia",
            "upia1", "upia2", "upib", "val", "value", "value_x", "value_y",
            "vjust", "x", "x_label", "xmax", "xmin", "y", "y_label", "ymax",
            "ymin", "z"),
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

expect_RcppAnnoy <- function(...) {
  rlang::check_installed("RcppAnnoy", ...)
}

expect_pcaMethods <- function(...) {
  rlang::check_installed("pcaMethods", ...)
}

expect_RcppML <- function(...) {
  rlang::check_installed("RcppML", ...)
}

expect_ggrepel <- function(...) {
  rlang::check_installed("ggrepel", ...)
}
