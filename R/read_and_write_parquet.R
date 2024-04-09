#' Read edgelists from .pxl files
#'
#' This function uses arrow to read edgelists from one or several
#' .pxl files. The edgelists are stored in parquet files in
#' the \code{outdir} directory which can be modified on disk.
#'
#' @param pxl_file Path to a .pxl file
#' @param edge_list_file Path to the output edgelist.parquet file
#' @param verbose Print messages
#' @param ... Parameters passed to other methods
#'
#' @import rlang
#' @importFrom arrow open_dataset
#' @importFrom utils unzip
#'
#' @return Nothing. The edgelist is saved to a parquet file
#' set with \code{edge_list_file}
#'
#' @examples
#' library(pixelatorR)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' edgelist_arrow <- ReadMPX_arrow_edgelist(pxl_file)
#' edgelist_arrow
#'
#' @export
#'
ReadMPX_arrow_edgelist <- function (
  pxl_file,
  edge_list_file = NULL,
  verbose = TRUE,
  ...
) {

  # Check input parameters
  stopifnot(
    "'pxl_file' must be a non-empty character of length 1" =
      inherits(pxl_file, what = "character") &&
      (length(pxl_file) == 1)
  )

  # Validate path
  if (!fs::file_exists(pxl_file)) {
    abort(glue("File '{pxl_file}' does not exist"))
  }

  available_files <- unzip(pxl_file, list = TRUE)$Name
  if (!"edgelist.parquet" %in% available_files) {
    abort(glue(".pxl file {filename} does not contain an 'edgelist.parquet' file"))
  }

  # Extract the edgelist parquet file
  if (!is.null(edge_list_file)) {
    stopifnot(
      "'edge_list_file' must be a non-empty character of length 1" =
        inherits(edge_list_file, what = "character") &&
        (length(edge_list_file) == 1),
      "'edge_list_file' must be a .parquet file" =
        .file_ext(edge_list_file) == "parquet"
    )
  } else {
    edge_list_dir <- fs::path_temp()
    if (verbose && check_global_verbosity())
      cli_alert_info("Extracting edgelist.parquet file to {col_br_blue(file.path(edge_list_dir, 'edgelist.parquet'))}")
  }

  # Extract the edgelist.parquet file
  zip::unzip(pxl_file, files = "edgelist.parquet", exdir = edge_list_dir)

  # Read the parquet file
  ds <- arrow::open_dataset(fs::path(edge_list_dir, "edgelist.parquet"))

  if (verbose && check_global_verbosity()) cli_alert_success("Returning FileSystemDataset")

  # Return list with additional info if return_list = TRUE
  return(ds)
}
