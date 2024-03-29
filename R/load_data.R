# Declarations used in package check
globalVariables(
  names = c('.'),
  package = 'pixelatorR',
  add = TRUE
)

#' Read a count matrix from a pxl file
#'
#' @param filename Path to a .pxl file
#' @param return_list If TRUE, returns a list with the expression matrix and
#' the output from \code{h5read}
#' @param verbose Print messages
#'
#' @import rlang
#' @importFrom rhdf5 h5read
#' @importFrom utils unzip
#'
#' @return A count matrix or a list if \code{return_list = TRUE}
#'
#' @examples
#' library(pixelatorR)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                          package = "pixelatorR")
#' counts <- ReadMPX_counts(pxl_file)
#' counts[1:5, 1:5]
#'
#' @export
#'
ReadMPX_counts <- function (
  filename,
  return_list = FALSE,
  verbose = TRUE
) {
  stopifnot("filename must be a character of length 1" = is.character(filename) & (length(filename) == 1))
  if (!file.exists(filename)) abort(glue("{filename} doesn't exist"))

  # Reads anndata from a .pxl file by unzipping it into a temporary folder and reading the anndata object within.
  if (verbose & check_global_verbosity())
    cli_alert_info("Loading count data from {filename}")

  # Unzip pxl file
  original_filename <- filename
  if (endsWith(filename, ".pxl")) {
    # LL: Returns error message
    filename <- tryCatch(unzip(filename, "adata.h5ad", exdir = tempdir()),
                         error = function(e) e,
                         warning = function(w) w)
    if (inherits(x = filename, what = "simpleWarning")) abort("Failed to unzip data")
  } else {
    abort(glue("Invalid file format .{.file_ext(filename)}. Expected a .pxl file."))
  }

  # Read temporary file
  anndata_hier <- h5read(filename, "/")

  # Remove temporary file
  if (endsWith(original_filename, ".pxl")) file.remove(filename)

  # Extract contents
  X <- anndata_hier$X
  colnames(X) <- anndata_hier$obs$component
  rownames(X) <- anndata_hier$var$marker

  if (return_list) {
    return(list(X = X, anndata_hier = anndata_hier))
  } else {
    return(X)
  }
}

#' Load data from PXL file into a \code{Seurat} object
#'
#' This wrapper function can be used to load data from a PXL file, and returns
#' a \code{Seurat} object with the desired data. By default, the MPX count matrix
#' is returned as an \code{Assay} and the graph data is excluded.
#'
#' The graph data is stored as an edgelist in the PXL file, and this edgelist
#' can be inconvenient to work with in memory for larger datasets. If you want to
#' have the edgelist available in your \code{Seurat} object, you can set \code{return_cellgraphassay = TRUE}
#' to return a \code{\link{CellGraphAssay}} class object instead. The \code{\link{CellGraphAssay}}
#' object extends the \code{Assay} class but can keep additional slots with graph-related data.
#' The edgelist is not loaded into memory, but is provided as an arrow Dataset stored in the
#' \code{\link{CellGraphAssay}}. This arrow Dataset makes it possible manipulate and fetch
#' information from the edgelist whenever necessary, without having to load it into memory.
#' See \code{\link{CellGraphAssay}} for more details.
#'
#' @param assay Assay name
#' @param return_cellgraphassay Should data be loaded as a \code{CellGraphAssay} object?
#' @param load_cell_graphs Should the cellgraphs be loaded into the \code{CellGraphAssay} object?
#' @param load_polarity_scores,load_colocalization_scores Logical specifying if the polarity
#' and colocalization scores should be loaded. These parameters only have an effect if
#' \code{return_cellgraphassay = TRUE}.
#' @param add_additional_assays If other matrix representations are stored in the PXL file,
#' for instance CLR-normalized counts or denoised, set this parameter to {TRUE} to load these
#' in separate Assays.
#' @param edgelist_outdir A directory where the edgelist should be stored
#' @param overwrite Should \code{edgelist_outdir} be overwritten?
#' @param ... Additional parameters passed to \code{\link{CreateSeuratObject}}
#' @inheritParams ReadMPX_counts
#'
#' @import rlang
#' @importFrom SeuratObject CreateSeuratObject CreateAssayObject `VariableFeatures<-`
#' @importFrom stats setNames
#'
#' @family data-loaders
#'
#' @return An object of class \code{Seurat}
#'
#' @examples
#'
#' library(pixelatorR)
#' # Set arrow data output directory to temp for tests
#' options(pixelatorR.arrow_outdir = tempdir())
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
#' seur_obj
#'
#' @export
#'
ReadMPX_Seurat <- function (
  filename,
  assay = "mpxCells",
  return_cellgraphassay = TRUE,
  load_cell_graphs = FALSE,
  load_polarity_scores = TRUE,
  load_colocalization_scores = TRUE,
  add_additional_assays = FALSE,
  edgelist_outdir = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {

  stopifnot(
    "assay must be a character of length 1" =
      is.character(assay) &&
      (length(assay) == 1)
  )

  # Load count matrix
  data <- ReadMPX_counts(filename = filename, return_list = TRUE, verbose = FALSE)
  X <- data$X; anndata_hier <- data$anndata_hier

  # Load edgelist
  empty_graphs <- rep(list(NULL), ncol(X)) %>% setNames(nm = colnames(X))

  if (return_cellgraphassay) {

    # Create CellGraphAssay
    cg_assay <- CreateCellGraphAssay(counts = X,
                                     cellgraphs = empty_graphs,
                                     arrow_dir = filename,
                                     outdir = edgelist_outdir,
                                     overwrite = overwrite)
    Key(cg_assay) <- paste0(assay, "_")

    if (load_cell_graphs) {
      if (verbose && check_global_verbosity())
        cli_alert_warning(glue(col_br_red("Loading cell graphs into the Seurat object may take time and uses a lot of memory.",
                                          " Consider using LoadCellGraphs to load selected cell graphs at a later stage.")))
      cg_assay <- LoadCellGraphs(cg_assay, verbose = verbose)
    }
    # Load polarity scores
    if (load_polarity_scores) {
      polarization <- ReadMPX_item(filename = filename, items = "polarization", verbose = FALSE)
      cg_assay@polarization <- polarization
    }
    # Load colocalization scores
    if (load_polarity_scores) {
      colocalization <- ReadMPX_item(filename = filename, items = "colocalization", verbose = FALSE)
      cg_assay@colocalization <- colocalization
    }
  } else {
    cg_assay <- CreateAssayObject(counts = X)
    Key(cg_assay) <- paste0(assay, "_")
  }

  # Create Seurat object
  seur_obj <- CreateSeuratObject(counts = cg_assay, assay = assay, ...)
  if (verbose && check_global_verbosity())
    cli_alert_success("Created a 'Seurat' object with {col_br_blue(ncol(X))} cells and {col_br_blue(nrow(X))} targeted surface proteins")

  if (add_additional_assays) {
    for (layr in names(anndata_hier$obsm)) {

      # Add only normalized counts as individual assays
      if (!layr %in% c("denoised", "normalized_clr", "normalized_rel")) next

      X_assay <-
        anndata_hier$obsm[[layr]] %>%
        as_tibble() %>%
        column_to_rownames("component") %>%
        t()
      layer_assay <- CreateAssayObject(counts = X_assay)
      suppressWarnings(seur_obj[[layr]] <- layer_assay)
    }
  }


  # Extract meta data
  seur_obj@meta.data <-
    anndata_hier$obs %>%
    {
      newdata <- .
      for (name in names(newdata)) {
        if (length(newdata[[name]]) == dim(anndata_hier$X)[2]) next
        indicies <- newdata[[name]]$codes + 1
        categories <- newdata[[name]]$categories
        newdata[[name]] <- factor(categories[indicies], categories)
      }
      newdata
    } %>%
    as.data.frame() %>%
    column_to_rownames("component")

  # Extract variable meta data (var)
  var_meta <- anndata_hier$var %>%
    as.data.frame()
  rownames(var_meta) <- var_meta$marker

  seur_obj[[assay]]@meta.features <- var_meta

  # Set variable features to all features
  VariableFeatures(seur_obj) <- rownames(seur_obj)

  return(seur_obj)
}

#' Read a pixel data item
#'
#' Function that unzips a .pxl file into a temporary folder and reads any
#' number of items (not anndata) as specified in items. If more than one
#' item is fetched, the output is a list of data frames. Otherwise the output
#' is a single \code{tbl_df}.
#'
#' @param filename Path to a .pxl file
#' @param items One of "colocalization", "polarization", "edgelist"
#' @param verbose Print messages
#'
#' @import rlang
#' @importFrom arrow read_parquet
#' @importFrom utils read.csv
#'
#' @family data-loaders
#'
#' @rdname ReadMPX_item
#'
#' @return Either an object of class \code{tbl_df} or a list of \code{tbl_df}
#' objects if multiple \code{items} are selected
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' polarization <- ReadMPX_item(pxl_file, items = "polarization")
#' polarization
#'
#' # Alternative 2
#' polarization <- ReadMPX_polarization(pxl_file)
#'
#' @export
#'
ReadMPX_item <- function (
  filename,
  items = c("colocalization", "polarization", "edgelist"),
  verbose = TRUE
) {

  # Check input parameters
  stopifnot(
    "filename must be a character of length 1" =
      is.character(filename) &&
      (length(filename) == 1),
    "Invalid items" =
      all(items %in% c("colocalization", "polarization", "edgelist")),
    "Expected a .pxl file" =
      .file_ext(filename) == "pxl"
  )
  if (!file.exists(filename)) abort(glue("{filename} doesn't exist"))

  if (verbose && check_global_verbosity())
    cli_alert_info("Loading item(s) from: {filename}")

  suppressWarnings(
    read_items <-
      items %>%
      set_names(., .) %>%
      lapply(function(item) {

        if (verbose && check_global_verbosity())
          cli_alert("  Loading {col_br_magenta(item)} data")

        # generate a file name for a new temporary directory.
        # When running R CMD check on windows latest release, we
        # do not have permissions to modify files in TEMPDIR.
        # For this reason, we create a new directory with a unique
        # name instead.
        exdir_temp <- fs::file_temp()

        # Unzip item to temporary directory
        unzipped_filename <- unzip(filename,  paste0(item, ".parquet"), exdir = exdir_temp)

        # Read contents of parquet file
        outdata <- read_parquet(unzipped_filename)

        # Try to delete temporary directory or throw a warning if it fails.
        check <- tryCatch({
          fs::dir_delete(exdir_temp)
        }, error = function(e)
          e,
        warning = function(w)
          w)
        if (inherits(check, what = "error"))
          cli_alert_warning("Failed to remove temporary directory {exdir_temp}")

        return(outdata)
      })
  )

  if (length(read_items) == 1) {
    read_items <- read_items[[1]]
  }

  if (verbose && check_global_verbosity()) {
    if (inherits(read_items, what = "list")) {
      cli_alert_success("Returning a list of tibbles with items: {paste(items, collapse = ', ')}")
    } else {
      cli_alert_success("Returning a '{class(read_items)[1]}' object")
    }
  }
  return(read_items)
}


#' @rdname ReadMPX_item
#'
#' @export
#'
ReadMPX_polarization <- function (
  filename,
  verbose = TRUE
) {
  ReadMPX_item(filename, items = "polarization", verbose = verbose)
}

#' @rdname ReadMPX_item
#'
#' @export
#'
ReadMPX_colocalization <- function (
  filename,
  verbose = TRUE
) {
  ReadMPX_item(filename, items = "colocalization", verbose = verbose)
}

#' @rdname ReadMPX_item
#'
#' @export
#'
ReadMPX_edgelist <- function (
  filename,
  verbose = TRUE
) {
  ReadMPX_item(filename, items = "edgelist", verbose = verbose)
}
