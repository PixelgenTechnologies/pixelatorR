#' @param file Name of the file where the object is written to
#' @param overwrite Set to TRUE to overwrite existing directory
#'
#' @import rlang
#' @import glue
#'
#' @rdname WriteMPX
#' @method WriteMPX CellGraphAssay
#'
#' @examples
#'
#' library(pixelatorR)
#' library(SeuratObject)
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
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

  arrow_data <- slot(object, name = "arrow_data")

  # If arrow_data is dead or missing, use saveRDS
  msg <- tryCatch(arrow_data %>% nrow(), error = function(e) "Error")
  msg <- msg %||% "Error"
  if (msg == "Error") {
    cli_alert_danger("Arrow data could not be retrieved. The edge list will not be saved.")
    saveRDS(object, file = file, ...)
  } else {
    cellgraphassay_dir <- dirname(path = file)
    if (!dir.exists(cellgraphassay_dir)) {
      abort(glue("Invalid output directory {cellgraphassay_dir}"))
    }

    # Export parquet file to the same directory as the CellGraphAssay object
    arrow_dir_path <-
      export_edgelist_to_parquet(
        object = arrow_data,
        outdir = cellgraphassay_dir,
        overwrite = overwrite,
        verbose = FALSE
      )

    # Update arrow_dir
    slot(object, name = "arrow_dir") <- arrow_dir_path

    # Export Seurat object
    saveRDS(object, file = file, ...)
  }
}


#' @param file Name of the file where the object is written to
#' @import rlang
#' @import glue
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
    abort("Found two CellGraphAssays in Seurat object. Cannot export more than one CellGraphAssay")

  # Only export Seurat object if there is no CellGraphAssay present
  if (sum(assay_classes == "CellGraphAssay") == 0) {
    cli_alert_warning("Found no CellGraphAssays. Exporitng Seurat object with saveRDS(...)")
    saveRDS(object, file = file, ...)
  } else {

    # Get directory name
    seurat_dir <- dirname(path = file)
    if (!dir.exists(seurat_dir)) {
      abort(glue("Invalid output directory {seurat_dir}"))
    }

    # Fetch arrow data
    arrow_data <- ArrowData(object)

    # If arrow_data is dead or missing, use saveRDS
    msg <- tryCatch(arrow_data %>% nrow(), error = function(e) "Error")
    msg <- msg %||% "Error"
    if (msg == "Error") {
      cli_alert_danger("Arrow data could not be retrieved. The edge list will not be saved.")
      saveRDS(object, file = file, ...)
    } else {
      # Export parquet file to the same directory as the Seurat object
      arrow_dir_path <-
        export_edgelist_to_parquet(
          object = arrow_data,
          outdir = seurat_dir,
          overwrite = overwrite,
          verbose = FALSE
        )

      # Update arrow_dir
      ArrowDir(object) <- arrow_dir_path

      # Export Seurat object
      saveRDS(object, file = file, ...)
    }
  }
}
