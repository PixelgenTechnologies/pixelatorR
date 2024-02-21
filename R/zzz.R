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

  # set arrow data output directory
  if (is.null(getOption("pixelatorR.arrow_outdir"))) {

    if (Sys.getenv("PIXELATORR_ARROWDIR") == "") {
      options(pixelatorR.arrow_outdir = file.path(getwd(), "edgelists"))
    } else {
      options(pixelatorR.arrow_outdir = Sys.getenv("PIXELATORR_ARROWDIR"))
    }
  }

}
