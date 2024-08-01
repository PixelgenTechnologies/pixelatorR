#' @include generics.R
NULL

#' @param pxl_files_dir The directory where the PXL files are stored
#' @param verbose Print messages
#' @param ... Additional arguments. Currently not used.
#'
#' @rdname RestorePaths
#' @method RestorePaths MPXAssay
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#'
#' # Copy PXL file to tempdir
#' tmp_pxl_file <- file.path(fs::path_temp(), "five_cells.pxl")
#' fs::file_copy(pxl_file, tmp_pxl_file)
#' seur_obj <- ReadMPX_Seurat(tmp_pxl_file)
#'
#' # Now we can load graphs
#' seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])
#'
#' # Removing or moving PXL file will make graphs inaccessible
#' fs::file_delete(tmp_pxl_file)
#' \dontrun{
#' # This will now fail since the PXL file is missing
#' seur_obj <- LoadCellGraphs(seur_obj, force = TRUE)
#' }
#'
#' # We can restore paths to the PXL file using RestorePaths.
#' # All we need is to provide the path to a directory where the PXL file is located.
#' pxl_files_dir <- system.file("extdata/five_cells",
#'                             package = "pixelatorR")
#' seur_obj <- RestorePaths(seur_obj, pxl_files_dir = pxl_files_dir)
#'
#' # Now LoadCellGraphs should work as expected
#' seur_obj <- LoadCellGraphs(seur_obj, force = TRUE)
#'
#' @export
#'
RestorePaths.MPXAssay <- function (
  object,
  pxl_files_dir,
  verbose = TRUE,
  ...
) {

  # Check that pxl_files_dir exists
  if (!fs::dir_exists(pxl_files_dir)) {
    abort(glue(
      "The directory '{cli::col_blue(pxl_files_dir)}' does not exist.",
      " Please provide a valid directory path."
    ))
  }

  fs_map <- FSMap(object)

  # Get PXL file names
  pxl_files <- fs_map %>% pull(all_of("pxl_file")) %>% basename()

  # Update pxl_files
  pxl_files_updated <- fs::path(pxl_files_dir, pxl_files)
  file_checks <- fs::file_exists(pxl_files_updated)
  if (!all(file_checks)) {
    abort(glue(
      "The directory '{cli::col_blue(pxl_files_dir)}' is missing the following files:\n",
      "{paste(pxl_files[!file_checks], collapse = '\n')}"
    ))
  }

  # Update pxl_files
  fs_map <- fs_map %>%
    mutate(pxl_file = pxl_files_updated %>% as.character())

  FSMap(object) <- fs_map

  if (verbose && check_global_verbosity())
    cli_alert_success(glue(
      "Successfully updated PXL file paths."
    ))

  return(object)
}

#' @rdname RestorePaths
#' @method RestorePaths CellGraphAssay
#' @docType methods
#' @export
#'
RestorePaths.CellGraphAssay <- RestorePaths.MPXAssay

#' @rdname RestorePaths
#' @method RestorePaths CellGraphAssay5
#' @docType methods
#' @export
#'
RestorePaths.CellGraphAssay5 <- RestorePaths.MPXAssay

#' @param assay Assay name
#'
#' @rdname RestorePaths
#' @method RestorePaths Seurat
#' @docType methods
#'
#' @export
#'
RestorePaths.Seurat <- function (
  object,
  pxl_files_dir,
  assay = NULL,
  verbose = TRUE,
  ...
) {

  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    stopifnot(
      "'assay' must be a character of length 1" =
        is.character(assay) &&
        (length(assay) == 1)
    )
  } else {
    # Use default assay if assay = NULL
    assay <- DefaultAssay(object)
  }

  cg_assay <- object[[assay]]
  cg_assay <- RestorePaths(cg_assay, pxl_files_dir = pxl_files_dir, verbose = verbose, ...)
  object[[assay]] <- cg_assay

  return(object)
}
