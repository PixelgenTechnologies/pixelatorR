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
#' pxl_file <- minimal_mpx_pxl_file()
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
#' pxl_file <- minimal_mpx_pxl_file()
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
#' library(pixelatorR)
#'
#' # Create example Seurat object
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#'
#' # Replace FSMap in Seurat object
#' FSMap(seur_obj) <- FSMap(seur_obj) %>%
#'   mutate(pxl_file = fs::path_rel(pxl_file))
#'
#' @export
#'
"FSMap<-.Seurat" <- function(
  object,
  ...,
  value
) {
  # Validate value
  assay <- DefaultAssay(object = object)
  pixel_assay <- object[[assay]]
  assert_class(pixel_assay, classes = c("CellGraphAssay", "CellGraphAssay5", "PNAAssay", "PNAAssay5"))
  FSMap(pixel_assay) <- value
  object[[assay]] <- pixel_assay
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


#' @param assay Name of a \code{PNAAssay}
#' @param meta_data_columns A character vector with meta.data column names.
#' This option can be useful to join meta.data columns with the proximity
#' score table.
#'
#' @method ProximityScores Seurat
#'
#' @rdname ProximityScores
#'
#' @examples
#' library(pixelatorR)
#'
#' # Create example Seurat object
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#'
#' # Get proximity scores
#' proximity <- ProximityScores(seur_obj)
#' proximity
#'
#' # Get proximity scores with additional meta data
#' proximity <-
#'   ProximityScores(seur_obj, meta_data_columns = c("n_umi", "n_edges"))
#' proximity
#'
#' # Get proximity scores with marker_1 and marker_2 counts
#' proximity <-
#'   ProximityScores(seur_obj, add_marker_counts = TRUE)
#' proximity
#'
#' @export
#'
ProximityScores.Seurat <- function(
  object,
  assay = NULL,
  meta_data_columns = NULL,
  add_marker_counts = FALSE,
  add_marker_proportions = FALSE,
  lazy = FALSE,
  calc_log2ratio = TRUE,
  ...
) {
  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)
  pixel_assay <- object[[assay]]
  assert_class(pixel_assay, classes = c("CellGraphAssay", "CellGraphAssay5", "PNAAssay", "PNAAssay5"))

  # Get proximity scores
  proximity_scores <- ProximityScores(pixel_assay, add_marker_counts, add_marker_proportions, lazy, calc_log2ratio, ...)

  if (inherits(proximity_scores, "tbl_lazy")) {
    con <- proximity_scores[[1]]$con
  } else {
    con <- NULL
  }

  # Add optional meta data columns to proximity scores table
  if (!is.null(meta_data_columns)) {
    mdata <- .meta_data_columns(object, con, lazy, meta_data_columns)

    # Add additional meta.data slots
    proximity_scores <- proximity_scores %>%
      left_join(
        y = mdata,
        by = "component"
      )
  }

  return(proximity_scores)
}


#' @method ProximityScores<- Seurat
#'
#' @rdname ProximityScores
#'
#' @examples
#' # Update proximity scores in Seurat object
#' ProximityScores(seur_obj) <- ProximityScores(seur_obj) %>%
#'   mutate(ratio = join_count / join_count_expected_mean)
#'
#' @export
#'
"ProximityScores<-.Seurat" <- function(
  object,
  assay = NULL,
  ...,
  value
) {
  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  # Update proximity scores in the assay
  pixel_assay <- object[[assay]]
  ProximityScores(pixel_assay) <- value
  object[[assay]] <- pixel_assay

  return(object)
}


#' @param assay Name of a \code{CellGraphAssay}
#' @param meta_data_columns A character vector with meta.data column names.
#' This option can be useful to join meta.data columns with the proximity
#' score table.
#'
#' @method Edgelists Seurat
#' @rdname Edgelists
#'
#' @export
#'
Edgelists.Seurat <- function(
  object,
  assay = NULL,
  meta_data_columns = NULL,
  lazy = TRUE,
  ...
) {
  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)
  pixel_assay <- object[[assay]]
  assert_class(pixel_assay, classes = c("CellGraphAssay", "CellGraphAssay5", "PNAAssay", "PNAAssay5"))

  # Get edgelist
  edgelists <- Edgelists(object[[assay]], lazy, ...)

  if (inherits(edgelists, "tbl_lazy")) {
    con <- edgelists$src$con
  } else {
    con <- NULL
  }

  # Add optional meta data columns to edgelist table
  if (!is.null(meta_data_columns)) {
    mdata <- .meta_data_columns(object, con, lazy, meta_data_columns)

    # Add additional meta.data slots
    edgelists <- edgelists %>%
      left_join(
        y = mdata,
        by = "component"
      )
  }

  return(edgelists)
}


#' Utility method to select meta data from a Suurat object
#'
#' @param object A Seurat object
#' @param con A database connection
#' @param lazy A logical indicating whether to lazy load the table
#' @param meta_data_columns A character vector with meta.data column names
#' @param call The calling environment
#'
#' @return A tibble with selected meta data columns
#'
#' @noRd
#'
.meta_data_columns <- function(
  object,
  con = NULL,
  lazy = FALSE,
  meta_data_columns,
  call = caller_env()
) {
  assert_class(object, "Seurat")
  assert_vector(meta_data_columns, type = "character", n = 1, call = call)
  assert_class(con, "duckdb_connection", allow_null = TRUE, call = call)
  assert_single_value(lazy, "bool")
  meta_data_columns_valid <- meta_data_columns %in% colnames(object[[]])
  if (any(!meta_data_columns_valid)) {
    cli::cli_abort(
      c(
        "x" = "The following meta data columns were not found in the meta.data slot: ",
        " " = "{.val {meta_data_columns[!meta_data_columns_valid]}}"
      )
    )
  }

  # Get selected meta data
  mdata <- object[[]] %>%
    select(all_of(meta_data_columns)) %>%
    as_tibble(rownames = "component")

  if (lazy) {
    copy_to(con, mdata)
    mdata <- tbl(con, "mdata")
  }

  return(mdata)
}
