#' Filter markers by isotype-relative and/or absolute abundance
#'
#' Identifies markers that are positive in a sufficient fraction of cells.
#' A cell is positive for a marker if that cell's counts-per-million (CPM)
#' clears **every** cutoff that is supplied:
#'
#' - If `isotype_ratio` is set, CPM must be greater than
#'   `isotype_ratio` times the median CPM of `isotype_markers`.
#' - If `abundance_threshold` is set, CPM must be greater than
#'   `abundance_threshold`.
#' - If both are set, both conditions must hold (AND). Neither argument
#'   overrides the other. At least one cutoff is required.
#'
#' A marker is kept if at least `min_cell_fraction` of cells are positive.
#' Isotype controls are always dropped from the result.
#'
#' When `group_column` is set, the isotype median and whether each cell is
#' positive are computed independently in each group.
#'
#' @param object A `Seurat` object with a `counts` layer.
#' @param isotype_markers Character vector of isotype control marker names
#'   (for example `c("mIgG1", "mIgG2a", "mIgG2b")`.
#' @param isotype_ratio Numeric relative cutoff versus the median isotype CPM,
#'   or `NULL` to skip this cutoff. Default is `1.5`.
#' @param abundance_threshold Numeric absolute CPM cutoff, or `NULL` to skip
#'   this cutoff. Default is `NULL`.
#' @param min_cell_fraction Minimum fraction of cells that must be positive
#'   for a marker to be kept. Default is `0.05`.
#' @param group_column Optional metadata column name. If provided, the
#'   function returns a named list with one result per group.
#' @param return_stats Logical; if `TRUE`, return a tibble of per-marker
#'   statistics (including a `kept` column) instead of marker names.
#'
#' @return If `group_column` is `NULL`: a character vector of kept marker names,
#'   or a tibble if `return_stats` is `TRUE`. If `group_column` is set: a
#'   named list of the same, one element per group.
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' seur <-
#'   ReadPNA_Seurat(
#'     minimal_pna_pxl_file(),
#'     overwrite = TRUE,
#'     load_proximity_scores = FALSE,
#'     verbose = FALSE
#'   )
#'
#' # Keep markers above 1.5 times the isotype median CPM in at least 5% of cells
#' kept <- FilterMarkers(
#'   object = seur,
#'   isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b")
#' )
#'
#' # Require both the isotype ratio and an absolute CPM floor
#' kept_and <- FilterMarkers(
#'   object = seur,
#'   isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
#'   isotype_ratio = 1.5,
#'   abundance_threshold = 2000,
#'   min_cell_fraction = 0.05
#' )
#'
#' # Per-group filtering
#' seur$sample <- c("S1", "S1", "S2", "S2", "S2")
#' kept_by_sample <- FilterMarkers(
#'   object = seur,
#'   isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
#'   group_column = "sample"
#' )
#'
#' @export
#'
FilterMarkers <- function(
  object,
  isotype_markers,
  isotype_ratio = 1.5,
  abundance_threshold = NULL,
  min_cell_fraction = 0.05,
  group_column = NULL,
  return_stats = FALSE
) {
  assert_class(object, "Seurat")
  assert_vector(isotype_markers, type = "character", n = 1)
  assert_x_in_y(isotype_markers, rownames(object))
  assert_single_value(isotype_ratio, type = "numeric", allow_null = TRUE)
  assert_single_value(abundance_threshold, type = "numeric", allow_null = TRUE)
  assert_single_value(min_cell_fraction, type = "numeric")
  assert_within_limits(min_cell_fraction, limits = c(0, 1))
  assert_single_value(group_column, type = "string", allow_null = TRUE)
  assert_col_in_data(group_column, object[[]], allow_null = TRUE)
  assert_single_value(return_stats, type = "bool")

  if (is.null(isotype_ratio) && is.null(abundance_threshold)) {
    cli::cli_abort(
      c(
        "i" = "At least one of {.arg isotype_ratio} or {.arg abundance_threshold} must be set.",
        "x" = "Both cutoffs were {.val NULL}."
      )
    )
  }

  assert_x_in_y("counts", Layers(object))
  counts_mat <- LayerData(object, layer = "counts")

  filter_one <- function(mat) {
    lib_sizes <- Matrix::colSums(mat)
    if (any(lib_sizes <= 0)) {
      cli::cli_abort(
        c(
          "i" = "All cells must have a positive library size in the counts layer.",
          "x" = "Found {sum(lib_sizes <= 0)} cell(s) with zero counts."
        )
      )
    }
    cpm <- Matrix::t(Matrix::t(mat) / lib_sizes) * 1e6

    isotype_cpm <- as.vector(as.matrix(cpm[isotype_markers, , drop = FALSE]))
    isotype_median_cpm <- stats::median(isotype_cpm)

    ratio_cut <- if (!is.null(isotype_ratio)) {
      isotype_ratio * isotype_median_cpm
    } else {
      -Inf
    }
    abs_cut <- if (!is.null(abundance_threshold)) {
      abundance_threshold
    } else {
      -Inf
    }
    cutoff <- max(ratio_cut, abs_cut)

    n_cells <- ncol(cpm)
    positive_fraction <- as.numeric(Matrix::rowSums(cpm > cutoff) / n_cells)
    kept <- (positive_fraction >= min_cell_fraction) &
      !(rownames(cpm) %in% isotype_markers)

    markers_stats <- tibble(
      marker = rownames(cpm),
      positive_fraction = unname(positive_fraction),
      isotype_median_cpm = isotype_median_cpm,
      isotype_ratio = if (is.null(isotype_ratio)) NA_real_ else isotype_ratio,
      abundance_threshold = if (is.null(abundance_threshold)) {
        NA_real_
      } else {
        abundance_threshold
      },
      kept = unname(kept)
    )

    if (isTRUE(return_stats)) {
      return(markers_stats)
    }
    return(markers_stats$marker[markers_stats$kept])
  }

  if (is.null(group_column)) {
    return(filter_one(counts_mat))
  }

  group_ids <- object[[]][[group_column]]
  unique_groups <- unique(group_ids)
  unique_groups <- unique_groups[!is.na(unique_groups)]
  if (!is.null(levels(group_ids))) {
    unique_groups <- levels(group_ids)[levels(group_ids) %in% unique_groups]
  } else {
    unique_groups <- sort(unique_groups)
  }

  results <- list()
  for (grp in unique_groups) {
    cells_in_group <- which(group_ids == grp)
    if (length(cells_in_group) == 0) {
      cli::cli_warn("Skipping empty group {.val {grp}}.")
      next
    }
    results[[as.character(grp)]] <- filter_one(
      counts_mat[, cells_in_group, drop = FALSE]
    )
  }

  return(results)
}
