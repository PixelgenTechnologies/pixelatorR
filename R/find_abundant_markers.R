#' Find abundant markers by isotype-relative and/or absolute abundance
#'
#' Identifies markers that are positive in a sufficient fraction of cells.
#' A cell is positive for a marker if that cell's counts-per-million (CPM)
#' clears **every** cutoff that is supplied. CPM is computed from counts
#' plus a pseudocount of 1, so no entry is zero:
#' \eqn{\mathrm{CPM}_{g,c} = 10^6 (x_{g,c} + 1) / \sum_{g'}(x_{g',c} + 1)},
#' where \eqn{x_{g,c}} is the count for marker \eqn{g} in cell \eqn{c}.
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
#' The isotype median is computed from the cells that enter a given call to
#' the filter. Without `group_column`, that is all cells in `object`. With
#' `group_column` (for example cell type), the median, positivity, and
#' `min_cell_fraction` rule are computed independently in each group so a
#' marker that is absent in one type does not pull the cutoff down for
#' another. Unused factor levels are omitted.
#'
#' `isotype_ratio` defaults to `1.5`, matching the relative cutoff used by
#' the original marker-filter helper. `min_cell_fraction` defaults to `0.05`
#' so markers that are abundant only in a small population can still be
#' kept when the function is applied to a whole sample.
#'
#' @param object A `Seurat` object with counts. Split Assay5 layers from
#'   `merge()` (for example `counts.1`, `counts.2`) are joined in a local
#'   copy via `JoinLayers()`; `object` is not modified.
#' @param isotype_markers Character vector of isotype control marker names
#'   (for example `c("mIgG1", "mIgG2a", "mIgG2b")`).
#' @param isotype_ratio Numeric relative cutoff versus the median isotype CPM,
#'   or `NULL` to skip this cutoff. Default is `1.5`.
#' @param abundance_threshold Numeric absolute CPM cutoff, or `NULL` to skip
#'   this cutoff. Default is `NULL`.
#' @param min_cell_fraction Minimum fraction of cells that must be positive
#'   for a marker to be kept. Default is `0.05`, small enough to keep markers
#'   that are abundant only in a minority population of a mixed sample.
#' @param group_column Optional metadata column name (typically cell type).
#'   If provided, the isotype median is computed within each group and the
#'   function returns a named list with one result per observed group.
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
#' kept <- FindAbundantMarkers(
#'   object = seur,
#'   isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b")
#' )
#'
#' # Require both the isotype ratio and an absolute CPM floor
#' kept_and <- FindAbundantMarkers(
#'   object = seur,
#'   isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
#'   isotype_ratio = 1.5,
#'   abundance_threshold = 2000,
#'   min_cell_fraction = 0.05
#' )
#'
#' # Per-group filtering (isotype median computed within each group)
#' seur$sample <- c("S1", "S1", "S2", "S2", "S2")
#' kept_by_sample <- FindAbundantMarkers(
#'   object = seur,
#'   isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
#'   group_column = "sample"
#' )
#'
#' @export
#'
FindAbundantMarkers <- function(
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

  count_layers <- grep("^counts", Layers(object), value = TRUE)
  if (length(count_layers) == 0L) {
    cli::cli_abort(
      c(
        "i" = "A counts layer is required.",
        "x" = "No layers matching {.val counts} were found in {.arg object}."
      )
    )
  }
  if (length(count_layers) > 1L || !identical(count_layers, "counts")) {
    object <- JoinLayers(object)
  }
  counts_mat <- LayerData(object, layer = "counts")

  filter_one <- function(mat) {
    raw_lib_sizes <- Matrix::colSums(mat)
    if (any(raw_lib_sizes <= 0)) {
      cli::cli_abort(
        c(
          "i" = "All cells must have a positive library size in the counts layer.",
          "x" = "Found {sum(raw_lib_sizes <= 0)} cell(s) with zero counts."
        )
      )
    }
    lib_sizes <- raw_lib_sizes + nrow(mat)
    cpm <- Matrix::t((Matrix::t(mat) + 1) / lib_sizes) * 1e6

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
    results[[as.character(grp)]] <- filter_one(
      counts_mat[, cells_in_group, drop = FALSE]
    )
  }

  return(results)
}
