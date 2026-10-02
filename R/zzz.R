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

#' Clean up cell-plot rgl resize watchers
#'
#' Cancels idle polling, deletes recorded temporary overlay textures, and
#' clears the device registry when the package namespace is unloaded.
#'
#' @param libpath Path to the package library being unloaded.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.onUnload <- function(libpath) {
  .cell_rgl_cancel_poll()
  keys <- ls(envir = .cell_rgl_chrome_registry)
  for (key in keys) {
    .cell_rgl_drop_chrome(key)
  }
  return(invisible(NULL))
}
