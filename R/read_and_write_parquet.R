#' Read edgelists from .pxl files
#'
#' This function uses arrow to read edgelists from one or several
#' .pxl files. The edgelists are stored in parquet files in
#' the \code{outdir} directory which can be modified on disk.
#'
#' @param path Path to a .pxl file or a directory containing .parquet files
#' @param outdir A path to an existing directory
#' @param return_list If TRUE, returns a list with an \code{ArrowObject},
#' a vector with sample IDs and a the path to the parquet output directory
#' @param overwrite If the output directory already exists, set this parameter
#' to TRUE to overwrite the directory
#' @param verbose Print messages
#'
#' @import rlang
#' @import glue
#' @import cli
#' @importFrom arrow open_dataset
#' @importFrom utils unzip
#'
#' @return An \code{ArrowObject}
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
#'                         package = "pixelatorR")
#' edgelist_arrow <- ReadMPX_arrow_edgelist(pxl_file, overwrite = TRUE)
#' edgelist_arrow
#'
#' @export
#'
ReadMPX_arrow_edgelist <- function (
  path,
  outdir = NULL,
  return_list = FALSE,
  overwrite = FALSE,
  verbose = TRUE
) {

  # Check input parameters
  stopifnot(
    "'path' must be a non-emty character of length 1" =
      inherits(path, what = "character") &&
      (length(path) == 1)
  )

  # Validate path
  if (dir.exists(path)) {
    # Check for parquet files if a directory is provided
    fs <- list.files(path, recursive = TRUE)
    stopifnot("Found no files in provided directory" = length(fs) > 0)
    for (f in fs) {
      if (.file_ext(f) != "parquet")
        abort(glue("'path' must provide a valid path to directory containing .parquet files"))
    }
  } else {
    # Check for pxl files if files are provided
    if (!file.exists(path)) abort(glue("file '{path}' doesn't exist"))
    if (.file_ext(path) != "pxl") abort(glue("'path' must provide a valid path to a .pxl file"))
  }

  # Use getOption("pixelatorR.arrow_outdir") if no outdir is provided
  if (is.null(outdir)) {
    outdir <- getOption("pixelatorR.arrow_outdir")
    if (verbose && check_global_verbosity())
      cli_alert_info("Output directory set to '{outdir}'")
  } else {
    stopifnot(
      "'outdir' must be a character vector of length 1" =
        is.character(outdir) &&
        (length(outdir) == 1)
    )
    outdir <- normalizePath(outdir)
    if (!dir.exists(outdir))  abort(glue("outdir '{outdir}' doesn't exist"))
  }

  # If pxl files are provided, make a copy of the parquet files stored internally
  # Otherwise, simply load the edgelist from the directory provided
  if (!dir.exists(path)) {

    # Create a new folder names by date, hour and minute
    session_tmpdir_random <- file.path(outdir, paste0(.generate_random_string(), "-", format(Sys.time(), "%Y-%m-%d-%H%M%S")))
    if (dir.exists(session_tmpdir_random) && !overwrite)
      abort(glue("Output directory '{session_tmpdir_random}' already exists.",
                 " Set 'overwrite=TRUE' to overwrite this directory."))

    # Create empty directory
    if (verbose && check_global_verbosity())
      cli_alert_info("Copying parquet files to {session_tmpdir_random}")
    suppressWarnings({dir.create(path = session_tmpdir_random)})

    # Copy parquet file
    newdir <- file.path(session_tmpdir_random, paste0("sample=S", 1))
    suppressWarnings({dir.create(path = newdir)})
    unzipped_filename <- unzip(path, paste0("edgelist.parquet"), exdir = newdir)

    # Open copied parquet file
    ds <- open_dataset(session_tmpdir_random)

    # Update path to session_tmpdir_random
    path <- session_tmpdir_random

  } else {

    # Open data set directly if a directory is provided
    ds <- open_dataset(path)

  }

  if (verbose && check_global_verbosity()) cli_alert_success("Returning FileSystemDataset")

  # Return list with additional info if return_list = TRUE
  if (return_list) {
    return(list(ArrowObject = ds, arrow_dir = path))
  } else {
    return(ds)
  }
}

#' Export a FileSystemDataset to .parquet files
#'
#' This function can be used to export a \code{FileSystemDataset} to a
#' folder containing .parquet files readable with \code{\link{open_dataset}}
#'
#' @param object An \code{FileSystemDataset}
#' @param outdir A path to an existing directory
#' @param suffix Add a suffix to output directory name
#' @param overwrite If the output directory already exists, set this parameter
#' to TRUE to overwrite the directory
#' @param verbose Print messages
#'
#' @import rlang
#' @import glue
#' @import cli
#' @importFrom arrow write_dataset
#' @importFrom utils unzip
#'
#' @return Path to the output folder
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
#'                         package = "pixelatorR")
#' edgelist_arrow <- ReadMPX_arrow_edgelist(pxl_file, overwrite = TRUE)
#'
#' # Manipulate edgelist with dplyr verbs
#' sorted_components <- edgelist_arrow %>%
#'   group_by(component) %>%
#'   summarize(component_size = n()) %>%
#'   arrange(-component_size)
#'
#' # Inspect results
#' sorted_components %>% collect()
#'
#' # Export results to a parquet file
#' tmp_dir <- tempdir()
#' outdir <- file.path(tmp_dir, "sorted_components")
#' dir.create(outdir, showWarnings = FALSE)
#' outpath <- export_edgelist_to_parquet(object = sorted_components,
#'                                       outdir = outdir,
#'                                       overwrite = TRUE)
#' outpath
#'
#' @export
#'
export_edgelist_to_parquet <- function (
  object,
  outdir,
  suffix = "",
  overwrite = FALSE,
  verbose = TRUE
) {
  stopifnot(inherits(object, what = c("arrow_dplyr_query", "FileSystemDataset", "tbl_df")))
  stopifnot("'outdir' must be a character vector of length 1" = is.character(outdir) & (length(outdir) == 1))
  outdir <- normalizePath(outdir)
  if (!dir.exists(outdir))  abort(glue("outdir '{outdir}' doesn't exist"))
  session_tmpdir_random <- file.path(outdir, paste0(.generate_random_string(), "-", format(Sys.time(), "%Y-%m-%d-%H%M%S")), suffix)

  if ((!overwrite) && dir.exists(session_tmpdir_random))
    abort(glue("Output directory '{session_tmpdir_random}' already exists. ",
               "Set 'overwrite=TRUE' to overwrite this directory."))
  if (verbose && check_global_verbosity())
    cli_alert_info("Saving edgelist to {session_tmpdir_random}")
  write_dataset(object %>% ungroup(), path = session_tmpdir_random, existing_data_behavior = "overwrite")
  return(session_tmpdir_random)
}
