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
#'
#' @return Nothing. The edgelist is saved to a parquet file
#' set with \code{edge_list_file}
#'
#' @examples
#' library(pixelatorR)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' edgelist_arrow <- ReadMPX_arrow_edgelist(pxl_file)
#' edgelist_arrow
#'
#' @export
#'
ReadMPX_arrow_edgelist <- function(
  pxl_file,
  edge_list_file = NULL,
  verbose = TRUE,
  ...
) {
  # Check input parameters
  assert_single_value(pxl_file, type = "string")
  assert_file_ext(pxl_file, ext = "pxl")
  assert_file_exists(pxl_file)

  available_files <- utils::unzip(pxl_file, list = TRUE)$Name
  if (!"edgelist.parquet" %in% available_files) {
    cli::cli_abort(
      c("x" = "{.str edgelist.parquet} is missing from {.file pxl_file}")
    )
  }

  # Extract the edgelist parquet file
  if (!is.null(edge_list_file)) {
    assert_single_value(edge_list_file, type = "string")
    assert_file_ext(edge_list_file, ext = "parquet")
    edge_list_dir <- dirname(edge_list_file)
    fname <- basename(edge_list_file)
  } else {
    edge_list_dir <- fs::path_temp()
    fname <- "edgelist.parquet"
    if (verbose && check_global_verbosity()) {
      cli_alert_info("Extracting edgelist.parquet file to {col_br_blue(file.path(edge_list_dir, 'edgelist.parquet'))}")
    }
  }

  # Extract the edgelist.parquet file
  utils::unzip(pxl_file, files = "edgelist.parquet", exdir = edge_list_dir)
  if (!is.null(edge_list_file)) {
    fs::file_move(fs::path(edge_list_dir, "edgelist.parquet"), fs::path(edge_list_dir, fname))
  }

  # Read the parquet file
  ds <- arrow::open_dataset(fs::path(edge_list_dir, fname))

  # Convert selected columns uint64 to string
  sel_fields <- grep(pattern = "^umi", names(ds), value = TRUE)
  sel_fields_type <- sapply(sel_fields, function(f) {
    arrow::schema(ds)$GetFieldByName(f)$ToString() %>%
      stringr::str_extract("(?<= ).*")
  })
  sel_fields <- sel_fields[sel_fields_type == "uint64"]
  if (length(sel_fields) > 0) {
    ds <- ds %>%
      mutate(across(all_of(sel_fields), as.character))
  }

  if (verbose && check_global_verbosity()) cli_alert_success("Returning {class(ds)[1]}")

  # Return list with additional info if return_list = TRUE
  return(ds)
}
