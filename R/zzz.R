#' @docType package
#' @name pixelatorR-package
#' @rdname pixelatorR-package
#'
"_PACKAGE"

.onLoad <- function(libname, pkgname) {

  # set verbosity
  if (is.null(getOption("pixelatorR.verbose"))) {
    options(pixelatorR.verbose = TRUE)
  }

}
