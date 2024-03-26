#' Check global option for verbosity
#'
#' By setting the global option \code{options(pixelatorR.verbose = FALSE)},
#' users can turn off all verbosity in pixelatorR. This function can be used to
#' check the value of this global options
#'
#' @noRd
check_global_verbosity <- function() {
  if (is.null(getOption("pixelatorR.verbose"))) return(TRUE)
  getOption("pixelatorR.verbose")
}

#' Generates a message once when the following
#' function ReadMPX_Seurat is called
#'
#' @noRd
.initial_call_message <- function() {
  cli_alert_info("Graph data will be stored on disk until loaded with LoadCellGraphs.")
  cli_text()
  cli_alert(glue("The directory where the graph data is stored can be ",
                 "set with \n the global option 'pixelatorR.arrow_outdir'."))
  cli_text()
  cli_alert(glue("Certain operations, such as subset and merge, will create new \n",
                 "directories in this folder. If the size of the directory exceeds \n",
                 "the size limit to trigger cleanup of unused files \n",
                 "(see global option 'pixelatorR.arrowdir_maxsize'), \n",
                 "these functions will trigger a cleanup to remove directories \n",
                 "that are not linked to any global variable in the current R session."))
  cli_text()
  cli_alert(glue("You can find more information about this behavior by \n",
                 "typing ?pixelatorR_options in the console."))
  cli_text()
  cli_alert_info(col_red("This message is displayed once per session."))
  options(pixelatorR.startup_message = FALSE)
}


#' Utility function to run clean up of edge list directories
#' from other functions
#'
#' @noRd
.run_clean <- function (
  default_size_gb = 5
) {

  edgelist_size <- suppressMessages(edgelist_directories_du())
  if (is.null(edgelist_size)) return(invisible(NULL))

  # Trigger garbage cleaning if the edgelist directories exceed the
  # maximum allowed size
  if (edgelist_size > getOption("pixelatorR.arrowdir_maxsize", fs::fs_bytes(default_size_gb*1024^3))) {
    cli_alert_info(
      glue(
        "Edge list directories exceed size limit to trigger cleanup of unused files ",
        "(see getOption('pixelatorR.arrowdir_maxsize')."
      )
    )
    attempt <- try({edgelist_directories_clean()})
    if (inherits(attempt, "try-error")) {
      cli_alert_warning(
        "Failed to clean edge list directories. "
      )
    }
  }
}


# ***********************************
# Check if certain packages are installed. If the package is missing, users will
# be asked to install the package.
# ***********************************

expect_jsonlite <- function(...) {
  rlang::check_installed('jsonlite', ...)
}

expect_scales <- function(...) {
  rlang::check_installed('scales', ...)
}

expect_graphlayouts <- function(...) {
  rlang::check_installed('graphlayouts', ...)
}

expect_pheatmap <- function(...) {
  rlang::check_installed('pheatmap', ...)
}

expect_duckdb <- function(...) {
  rlang::check_installed('duckdb', ...)
}
