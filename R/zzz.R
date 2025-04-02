#' @docType package
#' @name pixelatorR-package
#' @rdname pixelatorR-package
#'
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  run_on_load()
}

on_load({
  if (is.null(getOption("pixelatorR.verbose"))) {
    options(pixelatorR.verbose = TRUE)
  }
})
