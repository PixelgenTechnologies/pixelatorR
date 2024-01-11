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
