#' @param file Name of the file where the object is written to
#' @param overwrite Set to TRUE to overwrite existing directory
#'
#' @import rlang
#'
#' @rdname WriteMPX
#' @method WriteMPX CellGraphAssay
#'
#' @examples
#'
#' library(pixelatorR)
#' # Set arrow data output directory to temp for tests
#' options(pixelatorR.arrow_outdir = tempdir())
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/mock_data",
#'                         "mock_mpx_data.pxl",
#'                         package = "pixelatorR")
#' seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE, return_cellgraphassay = TRUE)
#'
#' # Save CellGraphAssay
#' outfile <- tempfile()
#' WriteMPX(object = seur_obj[["mpxCells"]], outfile, overwrite = TRUE)
#'
#' # Reload object
#' cg_assay <- readRDS(outfile)
#'
#' # The arrow_dir now points to the directory where the parquet files are stored
#' cg_assay@arrow_dir
#' list.files(cg_assay@arrow_dir)
#'
#' @export
#'
WriteMPX.CellGraphAssay <- function (
  object,
  file,
  overwrite = FALSE,
  ...
) {

  # Check input parameters
  stopifnot(
    "'file' must be a 'character' of length 1" =
      inherits(file, what = "character") && (length(file) == 1),
    "'overwrite' must be TRUE/FALSE" =
      is.logical(overwrite)
  )

  arrow_dir <- ArrowDir(object)

  # If arrow_data is dead or missing, use saveRDS
  if (!dir.exists(arrow_dir)) {
    abort(glue("ArrowDir '{arrow_dir}' is missing.\n",
               "If the R session was restarted or if the directory was moved, ",
               "you will have to rerun all steps from from the original pxl file ",
               "and re-export the CellGraphAssay object.\n\n",
               "Use saveRDS(object, file = <path to rds file>) to save this object without the edgelist data."))
  } else {
    cellgraphassay_dir <- dirname(path = file)
    if (!dir.exists(cellgraphassay_dir)) {
      abort(glue("Invalid output directory {cellgraphassay_dir}"))
    }

    # Create new directory
    file.copy(from = ArrowDir(object),
              to = cellgraphassay_dir,
              recursive = TRUE,
              overwrite = overwrite)

    # Update arrow_dir
    slot(object, name = "arrow_dir") <- file.path(cellgraphassay_dir, basename(ArrowDir(object)))

    # Export Seurat object
    saveRDS(object, file = file, ...)
  }
}


#' @param file Name of the file where the object is written to
#' @import rlang
#'
#' @rdname WriteMPX
#' @method WriteMPX Seurat
#'
#' @examples
#'
#' # Save Seurat object
#' outfile <- tempfile()
#' WriteMPX(object = seur_obj, outfile, overwrite = TRUE)
#'
#' # Reload object
#' seur_obj <- readRDS(outfile)
#'
#' @export
#'
WriteMPX.Seurat <- function (
  object,
  file,
  overwrite = FALSE,
  ...
) {

  # Check input parameters
  stopifnot(
    "'file' must be a 'character' of length 1" =
      inherits(file, what = "character") && (length(file) == 1),
    "'overwrite' must be TRUE/FALSE" =
      is.logical(overwrite)
  )

  # Check for CellGraphAssays
  assay_classes <- sapply(object@assays, class)
  if (sum(assay_classes == "CellGraphAssay") > 1)
    abort("Found more than one CellGraphAssays in Seurat object. Cannot export more than one CellGraphAssay")

  # Only export Seurat object if there is no CellGraphAssay present
  if (sum(assay_classes == "CellGraphAssay") == 0) {
    cli_alert_warning("Found no CellGraphAssays. Exporitng Seurat object with saveRDS(...)")
    saveRDS(object, file = file, ...)
  } else {

    # Get CellGraphAssay
    assay <- names(assay_classes)[which(assay_classes == "CellGraphAssay")]

    # Get directory name
    seurat_dir <- dirname(path = file)
    if (!dir.exists(seurat_dir)) {
      abort(glue("Invalid output directory {seurat_dir}"))
    }

    # Fetch arrow data
    arrow_dir <- ArrowDir(object, assay = assay)

    if (!dir.exists(arrow_dir)) {
      abort(glue("ArrowDir '{arrow_dir}' is missing.\n",
                 "If the R session was restarted or if the directory was moved, ",
                 "you will have to rerun all steps from from the original pxl file ",
                 "and re-export the Seurat object.\n\n",
                 "Use saveRDS(object, file = <path to rds file>) to save this object without the edgelist data."))
    } else {
      # Export parquet file to the same directory as the Seurat object
      file.copy(from = ArrowDir(object),
                to = seurat_dir,
                recursive = TRUE,
                overwrite = overwrite)

      # Update arrow_dir
      ArrowDir(object) <- file.path(seurat_dir, basename(ArrowDir(object)))

      # Export Seurat object
      saveRDS(object, file = file, ...)
    }
  }
}
