#' @include generics.R
NULL

#' @param verbose Print messages
#'
#' @rdname RestoreArrowConnection
#' @method RestoreArrowConnection CellGraphAssay
#'
#' @export
#'
RestoreArrowConnection.CellGraphAssay <- function (
  object,
  verbose = TRUE,
  ...
) {
  # Check for valid arrow_data
  arrow_data <- slot(object, name = "arrow_data")
  msg <- tryCatch(arrow_data$.data %>% nrow(), error = function(e) "Error")
  msg <- msg %||% "Error"
  if (msg == "Error") {
    if (length(slot(object, name = "arrow_dir")) == 0) {
      if (verbose && check_global_verbosity())
        cli_alert_info(glue("Cannot establish connection when 'arrow_dir' ",
                            "is missing. Returning unmodified object."))
      return(object)
    }

    # Read arrow data
    arrow_data <- ReadMPX_arrow_edgelist(
      path = slot(object, name = "arrow_dir"),
      outdir = slot(object, name = "arrow_dir"),
      verbose = FALSE
    )

    # Place arrow data in correct CelLGraphAssay slot
    slot(object, name = "arrow_data") <- arrow_data

    if (verbose && check_global_verbosity())
      cli_alert_success("Successfully restored connection. Returning updated object.")
    return(object)
  } else {
    if (verbose && check_global_verbosity())
      cli_alert_info("'arrow_data' is active. Returning unmodified object.")
    return(object)
  }
}

#' @param verbose Print messages
#'
#' @rdname RestoreArrowConnection
#' @method RestoreArrowConnection Seurat
#'
#' @export
#'
RestoreArrowConnection.Seurat <- function (
  object,
  verbose = TRUE,
  ...
) {
  # Check for CellGraphAssays
  assay_classes <- sapply(object@assays, class)
  if (!any(assay_classes == "CellGraphAssay")) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("Found no 'CellGraphAssays' in 'Seurat' object. Returning unmodified object.")
      return(object)
    }
  }

  # handle CellGraphAssays
  for (assay_name in names(object@assays[assay_classes == "CellGraphAssay"])) {
    if (verbose && check_global_verbosity())
      cli_alert_info("Restoring connection for CellGraphAssay '{assay_name}'")
    cg_assay <- object[[assay_name]]
    cg_assay <- RestoreArrowConnection(cg_assay, verbose = FALSE)
    object[[assay_name]] <- cg_assay
  }
  if (verbose && check_global_verbosity())
    cli_alert_success("Successfully resotered connections.")
  return(object)
}
