# Declarations used in package check
globalVariables(
  names = c(
    ".", ".", ".score_rank", ".x", ".x", "alias", "bi_prob", "bi_prob", 
    "cell_count", "color", "color", "community", "community", "comp", 
    "comp", "component", "component", "component_new", "component_new", 
    "componentId", "componentId", "count_1", "count_1", "count_2", 
    "count_2", "current_id", "current_id", "degree", "degree", "dens", 
    "dens", "dev_png", "dev_png", "doublet_nn_rate", "doublet_nn_rate", 
    "doublet_nns", "doublet_nns", "doublet_p", "doublet_p", "doublet_p_adj", 
    "doublet_p_adj", "edge", "edge", "edges", "edges", "f1", "f1", 
    "f2", "f2", "frac", "frac", "frequency", "frequency", "from", 
    "from", "g", "g", "graph_edges", "graph_edges", "graph_projection", 
    "graph_projection", "graph_proteins", "graph_proteins", "graph_reads", 
    "graph_reads", "group", "group", "hjust", "hjust", "hp", "hp", 
    "hue", "hue", "hup", "hup", "i", "i1", "i1", "i2", "i2", "id", 
    "id", "id_map", "id_map", "in_gate", "in_gate", "index", "index", 
    "item", "item", "j", "join_count", "join_count", "join_count_expected_mean", 
    "join_count_expected_mean", "join_count_expected_mean_list", 
    "join_count_expected_mean_list", "join_count_list", "join_count_list", 
    "label", "label", "label_text", "label_text", "layout", "layout", 
    "lbl", "lbl", "level", "level", "marker", "marker", "marker_1", 
    "marker_1", "marker_1_ordered", "marker_1_unordered", "marker_1_unordered", 
    "marker_2", "marker_2", "marker_2_ordered", "marker_2_unordered", 
    "marker_2_unordered", "marker_x", "marker_x", "marker_y", "marker_y", 
    "marker1", "marker1", "marker2", "marker2", "med_ref", "med_ref", 
    "med_tgt", "med_tgt", "min_n_obs", "mixture_component", "mixture_component", 
    "modality", "modality", "molecules", "molecules", "morans_z", 
    "morans_z", "n", "n", "n_cells", "n_cells", "n_cells_detected", 
    "n_cells_detected", "n_cells_missing", "n_cells_missing", "n_edges", 
    "n_edges", "n_groups", "n_inside", "n_inside", "n_nodes", "n_nodes", 
    "n_ref", "n_ref", "n_tgt", "n_tgt", "name", "name", "neighbor", 
    "neighbor", "network", "node_type", "node_type", "node_val", 
    "node_val", "nodes", "nodes", "norm_factor", "norm_factor", "nties", 
    "nties", "nties_const", "nties_const", "original_id", "original_id", 
    "p", "p", "p_adj", "p_adj", "p_val", "p_val", "p_val_adj", "p_val_adj", 
    "p.value", "p.value", "p1", "p1", "p2", "p2", "pair", "pair", 
    "patch", "patch", "pct.1", "pct.1", "pct.2", "pct.2", "pearson_z", 
    "pearson_z", "png", "png", "pxl_file", "pxl_file", "quadrant", 
    "quadrant", "r", "r", "read_count", "read_count", "reads", "reads", 
    "reads_in_component", "reads_in_component", "receiver_freq", 
    "receiver_freq", "receiver_unmixed_freq", "receiver_unmixed_freq", 
    "ref_n", "ref_n", "ref_present", "references", "rn", "rn", "rs_ref", 
    "rs_ref", "rs_tgt", "rs_tgt", "score", "sigma", "sigma", "simulated", 
    "simulated", "size", "size", "sources", "string_protein_id", 
    "swap", "target_freq", "target_freq", "target_n", "target_n", 
    "target_unmixed_freq", "target_unmixed_freq", "tau", "tau", "tau_type", 
    "tau_type", "text_color", "text_color", "theoretical_max_edges", 
    "theoretical_max_edges", "theoretical_max_nodes", "theoretical_max_nodes", 
    "to", "to", "total", "total", "tp", "tp", "trial", "trial", "tup", 
    "tup", "type", "type", "u", "u", "uei_count", "uei_count", "umi_category", 
    "umi_category", "umi_count", "umi_count", "umi_high", "umi_high", 
    "umi_low", "umi_low", "umi_per_upia", "umi_per_upia", "umi1", 
    "umi1", "umi2", "umi2", "uniprot", "uniprot_a", "uniprot_a_ordered", 
    "uniprot_b", "uniprot_b_ordered", "uniprot_id", "upia", "upia", 
    "upia1", "upia1", "upia2", "upia2", "upib", "upib", "val", "val", 
    "value", "value", "value_x", "value_x", "value_y", "value_y", 
    "vjust", "vjust", "weight", "weight", "x", "x", "x_label", "x_label", 
    "xmax", "xmax", "xmin", "xmin", "y", "y", "y_label", "y_label", 
    "y_position", "y_position", "ymax", "ymax", "ymin", "ymin", "z", 
    "z"
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

expect_RSpectra <- function(...) {
  rlang::check_installed("RSpectra", ...)
}

expect_pls <- function(...) {
  rlang::check_installed("pls", ...)
}

expect_matrixStats <- function(...) {
  rlang::check_installed("matrixStats", ...)
}

expect_FNN <- function(...) {
  rlang::check_installed("FNN", ...)
}

expect_sparseMatrixStats <- function(...) {
  rlang::check_installed("sparseMatrixStats", ...)
}
