#' @include generics.R
NULL


#' @rdname CellGraphs
#' @method CellGraphs Seurat
#' @export
#' @concept cellgraphs
#'
#' @examples
#'
#' # CellGraphs getter Seurat
#' # ---------------------------------
#' se <- ReadMPX_Seurat(pxl_file)
#'
#' # Get cellgraphs from a Seurat object
#' CellGraphs(se)
#'
CellGraphs.Seurat <- function(
  object,
  ...
) {
  assay <- DefaultAssay(object = object)
  return(CellGraphs(object = object[[assay]]))
}

#' @method FSMap Seurat
#'
#' @rdname FSMap
#'
#' @examples
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' seur_obj <- ReadMPX_Seurat(pxl_file)
#'
#' # Check PXL file paths in a Seurat object
#' FSMap(seur_obj)
#'
#' @export
#'
FSMap.Seurat <- function(
  object,
  ...
) {
  assay <- DefaultAssay(object = object)
  FSMap(object[[assay]])
  return(FSMap(object[[assay]]))
}

#' @method FSMap<- Seurat
#'
#' @rdname FSMap
#'
#' @examples
#' \dontrun{
#' # Update PXL file paths in a Seurat object
#' new_pxl_file_paths <- c("<path_to_pxl_file>")
#' FSMap(seur_obj) <- FSMap(seur_obj) %>%
#'   mutate(pxl_file = new_pxl_file_paths)
#' }
#'
#' @export
#'
"FSMap<-.Seurat" <- function(
  object,
  ...,
  value
) {
  # Validate value
  .validate_fs_map(value)
  assay <- DefaultAssay(object = object)
  cg_assay <- object[[assay]]
  assert_mpx_assay(cg_assay)
  slot(cg_assay, name = "fs_map") <- value
  object[[assay]] <- cg_assay
  return(object)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @export
#' @method CellGraphs<- Seurat
#' @rdname CellGraphs
#' @concept cellgraphs
#'
#' @examples
#' # CellGraphs setter Seurat
#' # ---------------------------------
#'
#' # Set cellgraphs in a Seurat object
#' CellGraphs(se) <- cg_assay@cellgraphs
#'
"CellGraphs<-.Seurat" <- function(
  object,
  ...,
  value
) {
  assay <- DefaultAssay(object = object)
  CellGraphs(object = object[[assay]]) <- value
  return(object)
}


#' @param assay Name of a \code{CellGraphAssay}
#' @param meta_data_columns A character vector with meta.data column names.
#' This option can be useful to join meta.data columns with the polarization
#' score table.
#' @param add_marker_counts A logical value indicating whether to add marker
#' counts to the polarization score table.
#'
#' @method PolarizationScores Seurat
#'
#' @rdname PolarizationScores
#'
#' @export
#'
PolarizationScores.Seurat <- function(
  object,
  assay = NULL,
  meta_data_columns = NULL,
  add_marker_counts = FALSE,
  ...
) {
  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)
  cg_assay <- object[[assay]]
  assert_mpx_assay(cg_assay)

  # Get polarizaation scores from CellGraphAssay
  pol_scores <- PolarizationScores(cg_assay, add_marker_counts)

  # Handle adding meta data columns
  if (!is.null(meta_data_columns)) {
    assert_vector(meta_data_columns, type = "character", n = 1)
    meta_data_columns_valid <- meta_data_columns %in% colnames(object[[]])
    if (any(!meta_data_columns_valid)) {
      cli::cli_abort(
        c(
          "x" = "The following meta data columns were not found in the meta.data slot: ",
          " " = "{.val {meta_data_columns[!meta_data_columns_valid]}}"
        )
      )
    }

    # Add additional meta.data slots
    pol_scores <- pol_scores %>%
      left_join(
        y = object[[]] %>%
          as_tibble(rownames = "component") %>%
          select(component, all_of(meta_data_columns)),
        by = "component"
      )
  }

  return(pol_scores)
}


#' @method PolarizationScores<- Seurat
#'
#' @rdname PolarizationScores
#'
#' @export
#'
"PolarizationScores<-.Seurat" <- function(
  object,
  assay = NULL,
  ...,
  value
) {
  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  # Fetch polarization scores from assay
  cg_assay <- object[[assay]]
  PolarizationScores(cg_assay) <- value
  object[[assay]] <- cg_assay

  return(object)
}


#' @param assay Name of a \code{CellGraphAssay}
#' @param meta_data_columns A character vector with meta.data column names.
#' This option can be useful to join meta.data columns with the colocalization
#' score table.
#'
#' @method ColocalizationScores Seurat
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
ColocalizationScores.Seurat <- function(
  object,
  assay = NULL,
  meta_data_columns = NULL,
  add_marker_counts = FALSE,
  ...
) {
  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)
  cg_assay <- object[[assay]]
  assert_mpx_assay(cg_assay)

  # Get colocalization scores
  coloc_scores <- ColocalizationScores(cg_assay, add_marker_counts)

  # Handle adding meta data columns
  if (!is.null(meta_data_columns)) {
    assert_vector(meta_data_columns, type = "character", n = 1)
    meta_data_columns_valid <- meta_data_columns %in% colnames(object[[]])
    if (any(!meta_data_columns_valid)) {
      cli::cli_abort(
        c(
          "x" = "The following meta data columns were not found in the meta.data slot: ",
          " " = "{.val {meta_data_columns[!meta_data_columns_valid]}}"
        )
      )
    }

    # Add additional meta.data slots
    coloc_scores <- coloc_scores %>%
      left_join(
        y = object[[]] %>%
          as_tibble(rownames = "component") %>%
          select(component, all_of(meta_data_columns)),
        by = "component"
      )
  }

  return(coloc_scores)
}


#' @method ColocalizationScores<- Seurat
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
"ColocalizationScores<-.Seurat" <- function(
  object,
  assay = NULL,
  ...,
  value
) {
  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  # Fetch polarization scores from assay
  cg_assay <- object[[assay]]
  ColocalizationScores(cg_assay) <- value
  object[[assay]] <- cg_assay

  return(object)
}
