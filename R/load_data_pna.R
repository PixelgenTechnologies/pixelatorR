#' Read a count matrix from a PXL file with PNA data
#'
#' @param pxl_file Path to a PXL file
#'
#' @family PXL-data-loaders
#'
#' @return A count matrix or a list if \code{return_list = TRUE}
#'
#' @examples
#' library(pixelatorR)
#'
#' # Load example counts
#' pxl_file <- minimal_pna_pxl_file()
#' counts <- ReadPNA_counts(pxl_file)
#' counts[1:5, 1:5]
#'
#' @export
#'
ReadPNA_counts <- function(
  pxl_file
) {
  assert_single_value(pxl_file, type = "string")
  assert_file_ext(pxl_file, ext = "pxl")
  assert_file_exists(pxl_file)

  # Fetch counts from data base
  db <- PixelDB$new(pxl_file)
  on.exit(db$close())
  X <- db$counts()

  return(X)
}


#' Load data from PNA PXL file into a \code{Seurat} object
#'
#' This function can be used to load data from a PXL file, and returns
#' a \code{Seurat} object.
#'
#' By default, the count matrix is returned in a \code{PNAAssay} object.
#' If the global option \code{Seurat.object.assay.version} is set to \code{"v5"},
#' the function will return a \code{PNAAssay5} object instead.
#'
#' @param assay Assay name
#' @param return_pna_assay Logical specifying whether the count data should be
#' stored in a \code{PNAAssay}/\code{PNAAssay5}. If set to \code{FALSE}, the count
#' data will be stored in a \code{Assay}/\code{Assay5} instead. For the latter case,
#' many features provided in \code{pixelatorR} will be unavailable. This can be useful
#' if you only intend to analyze the abundance data.
#' @param load_proximity_scores Logical specifying whether the proximity scores should
#' be loaded into the \code{PNAAssay}/\code{PNAAssay5}. If you only intend to analyze
#' abundance data or PNA graphs, you can set this parameter to \code{FALSE} to use less
#' memory. This parameter only have an effect if \code{return_pna_assay = TRUE}.
#' @param ... Additional parameters passed to \code{\link[SeuratObject]{CreateSeuratObject}}
#' @inheritParams ReadPNA_counts
#' @inheritParams ReadPNA_proximity
#'
#' @family PXL-data-loaders
#'
#' @examples
#' library(pixelatorR)
#'
#' # Crete example Seurat object
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' seur_obj
#'
#' @return An object of class \code{Seurat}
#'
#' @export
#'
ReadPNA_Seurat <- function(
  pxl_file,
  assay = "PNA",
  return_pna_assay = TRUE,
  load_proximity_scores = TRUE,
  calc_log2_ratio = TRUE,
  verbose = TRUE,
  ...
) {
  # Validate input parameters
  assert_single_value(pxl_file, type = "string")
  assert_file_exists(pxl_file)
  assert_file_ext(pxl_file, ext = "pxl")
  assert_single_value(assay, type = "string")
  assert_single_value(return_pna_assay, type = "bool")
  assert_single_value(load_proximity_scores, type = "bool")
  assert_single_value(calc_log2_ratio, type = "bool")
  assert_single_value(verbose, type = "bool")

  # Setup connection to PXL file
  db <- PixelDB$new(pxl_file)
  on.exit(db$close())

  # Load count matrix from database
  X <- db$counts()

  # Create an empty cellgraphs list
  empty_graphs <- rep(list(NULL), ncol(X)) %>% set_names(nm = colnames(X))

  if (return_pna_assay) {
    # Create CellGraphAssay(5)
    pna_assay_create_func <- switch(getOption("Seurat.object.assay.version", "v3"),
      "v3" = CreatePNAAssay,
      "v5" = CreatePNAAssay5
    )
    pna_assay <- pna_assay_create_func(
      counts = X,
      cellgraphs = empty_graphs,
      fs_map = tibble(
        id_map = list(tibble(
          current_id = colnames(X),
          original_id = colnames(X)
        )),
        sample = 1L,
        pxl_file = ifelse(.is_absolute_path(pxl_file), normalizePath(pxl_file), pxl_file)
      )
    )
    Key(pna_assay) <- paste0(assay, "_")

    # Load proximity scores
    if (load_proximity_scores) {
      proximity <- db$proximity(calc_log2_ratio)
      if (all(c("marker1", "marker2") %in% names(proximity))) {
        proximity <- proximity %>%
          rename(marker_1 = !!sym("marker1"), marker_2 = !!sym("marker2"))
      }
      pna_assay@proximity <- proximity
    }
  } else {
    pna_assay <- switch(getOption("Seurat.object.assay.version", "v3"),
      "v3" = CreateAssayObject(counts = X),
      "v5" = CreateAssay5Object(counts = X)
    )
    Key(pna_assay) <- paste0(assay, "_")
  }

  # Create Seurat object
  seur_obj <- CreateSeuratObject(counts = pna_assay, assay = assay, ...)
  if (verbose && check_global_verbosity()) {
    cli_alert_success(
      c(
        "Created a {.cls Seurat} object with {.val {ncol(X)}} cells ",
        "and {.val {nrow(X)}} targeted surface proteins"
      )
    )
  }

  # Extract meta data
  meta_data <- db$cell_meta()

  if (!all(rownames(meta_data) == colnames(seur_obj))) {
    cli::cli_abort(
      c(
        "i" = "The cell IDs in the {.cls Seurat} object and the metadata",
        " " = "collected from the adata.h5ad file do not match. ",
        "x" = "Failed to create {.cls Seurat} object. The PXL file appears to be corrupt."
      )
    )
  }
  seur_obj@meta.data <- meta_data

  # Extract feature meta data (var)
  feature_meta_data <- db$protein_meta()

  if (getOption("Seurat.object.assay.version", "v3") == "v3") {
    if (!"marker" %in% colnames(feature_meta_data)) {
      feature_meta_data <- feature_meta_data %>%
        mutate(marker = rownames(.))
    }

    # Replace underscore with dash
    feature_meta_data$marker <- gsub("_", "-", feature_meta_data$marker)
    rownames(feature_meta_data) <- feature_meta_data$marker
    check <- try({
      seur_obj[[assay]]@meta.features <- feature_meta_data
    })
  } else {
    check <- try({
      seur_obj[[assay]]@meta.data <- feature_meta_data
    })
  }
  if (inherits(check, "try-error")) {
    warn(glue(
      "Failed to set feature meta data."
    ))
  }

  # Set variable features to all features
  VariableFeatures(seur_obj) <- rownames(seur_obj)

  # Remove redundant meta data columns
  seur_obj@meta.data <- seur_obj@meta.data %>% select(-contains(c("nFeature", "nCount", "_index")))

  return(seur_obj)
}


#' Read metadata from a PNA PXL file
#'
#' @param pxl_file Path to a PXL file
#'
#' @family PXL-data-loaders
#'
#' @return A \code{tbl_df} with sample meta data
#'
#' @examples
#' library(pixelatorR)
#'
#' # Create example Seurat object
#' pxl_file <- minimal_pna_pxl_file()
#' ReadPNA_metadata(pxl_file)
#'
#' @export
#'
ReadPNA_metadata <- function(
  pxl_file
) {
  assert_file_exists(pxl_file)
  assert_file_ext(pxl_file, ext = "pxl")

  # Read meta data from database
  db <- PixelDB$new(pxl_file)
  on.exit(db$close())
  meta_data <- db$run_meta()
  if (nrow(meta_data) == 0) {
    cli::cli_abort(
      c("x" = "No metadata found in the PXL file")
    )
  }

  return(meta_data)
}

#' Load the edge list from a PNA PXL file
#'
#' Note that the umi1 and umi2 sequences are encoded as \code{integer64} which
#' is not natively supported by R. The \code{bit64} package is required to
#' read these IDs and conversion to \code{integer} or \code{numeric} should
#' be avoided. It is safer to convert these to character vectors if needed
#' for downstream processing.
#'
#' @param pxl_file Path to a PXL file with a PNA data set
#' @param cells A character vector with component IDs. If NULL, all components are loaded.
#' @param umi_data_type One of "int64", "string" or "suffixed_string". Default is "int64".
#'  - "int64": The UMIs are encoded as int64
#'  - "string": The UMIs are encoded as character
#'  - "suffixed_string": The UMIs are encoded as character with a suffix '-umi1' or '-umi2' added
#'
#' @export
#'
ReadPNA_edgelist <- function(
  pxl_file,
  cells = NULL,
  umi_data_type = c("int64", "string", "suffixed_string")
) {
  assert_single_value(pxl_file, type = "string")
  assert_file_exists(pxl_file)
  assert_file_ext(pxl_file, ext = "pxl")
  assert_vector(cells, type = "character", allow_null = TRUE, n = 1)

  # Read edgelist from database
  db <- PixelDB$new(pxl_file)
  on.exit(db$close())

  available_components <- db$counts() %>% colnames()
  if (!any(cells %in% available_components) && !is.null(cells)) {
    cli::cli_abort(
      c(
        "i" = "All {.var cells} must be present in the PXL file",
        "x" = "The following are missing: {.val {setdiff(cells, available_components)}}"
      )
    )
  }

  el <- db$components_edgelist(cells, umi_data_type, include_all_columns = TRUE) %>%
    as_tibble()

  return(el)
}

#' Load layouts from a PNA PXL file
#'
#' @param pxl_file Path to a PXL file containing PNA data
#' @param cells A character vector with component IDs. If NULL, all components are loaded.
#' @param add_marker_counts Logical specifying if marker counts should be added to the
#' layout table(s).
#' @param verbose Logical specifying if verbose output should be printed
#'
#' @return A list with one \code{tbl_df} object for each component. Each object contains the
#' node IDs and their x, y, and z coordinates.
#'
#' @export
#'
ReadPNA_layouts <- function(
  pxl_file,
  cells = NULL,
  add_marker_counts = FALSE,
  verbose = FALSE
) {
  assert_single_value(pxl_file, type = "string")
  assert_file_exists(pxl_file)
  assert_file_ext(pxl_file, ext = "pxl")
  assert_vector(cells, type = "character", allow_null = TRUE, n = 1)

  # Read layouts from database
  db <- PixelDB$new(pxl_file)
  on.exit(db$close())

  # Check that pre-computed layouts exist
  if (!"layouts" %in% db$names()) {
    cli::cli_abort(
      c("x" = "No layouts found in the PXL file")
    )
  }

  available_components <- db$counts() %>% colnames()
  if (!any(cells %in% available_components) && !is.null(cells)) {
    cli::cli_abort(
      c(
        "i" = "All {.var cells} must be present in the PXL file",
        "x" = "The following are missing: {.val {setdiff(cells, available_components)}}"
      )
    )
  }
  layouts <- db$components_layout(cells, add_marker_counts, verbose)

  return(layouts)
}


#' Load the Proximity scores table from a PNA PXL file
#'
#' @param pxl_file Path to a PXL file containing PNA data
#' @param calc_log2_ratio A logical specifying whether to calculate and add
#' a log2ratio column to the output table. Default is \code{TRUE}
#' @param verbose Print messages
#'
#' @family data-loaders
#'
#' @rdname ReadPNA_proximity
#'
#' @return A \code{tbl_df} or a \code{Table} with PNA Proximity scores
#'
#' @examples
#' library(pixelatorR)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' proximity_tbl <- ReadPNA_proximity(pxl_file)
#' proximity_tbl
#'
#' @export
#'
ReadPNA_proximity <- function(
  pxl_file,
  calc_log2_ratio = TRUE,
  verbose = TRUE
) {
  # Check input parameters
  assert_single_value(pxl_file, type = "string")
  assert_file_ext(pxl_file, ext = "pxl")
  assert_file_exists(pxl_file)
  assert_single_value(calc_log2_ratio, type = "bool")
  assert_single_value(verbose, type = "bool")

  if (verbose && check_global_verbosity()) {
    cli_alert_info("Loading proximity scores from: {.file {pxl_file}}")
  }

  # Setup connection to PXL file
  db <- PixelDB$new(pxl_file)
  on.exit(db$close())

  # Load proximity scores from database
  proximity_scores <- db$proximity(calc_log2_ratio)

  if (verbose && check_global_verbosity()) {
    cli_alert_success("Returning a {.cls {class(proximity_scores)[1]}} with the PNA proximity scores")
  }

  return(proximity_scores)
}
