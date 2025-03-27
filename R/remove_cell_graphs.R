#' @rdname RemoveCellGraphs
#' @method RemoveCellGraphs MPXAssay
#'
#' @export
#'
RemoveCellGraphs.MPXAssay <- function(
  object,
  ...
) {
  slot(object, name = "cellgraphs") <- rep(list(NULL), ncol(object)) %>% set_names(nm = colnames(object))
  return(object)
}

#' Remove cell graphs from a \code{PNAssay} or \code{PNAssay5} object
#'
#' @param object An object with cell graphs
#' @param ... Additional arguments (not used)
#'
#' @rdname RemoveCellGraphs
#' @method RemoveCellGraphs PNAAssay
#'
#' @export
#'
RemoveCellGraphs.PNAAssay <- function(
    object,
    ...
) {
  slot(object, name = "cellgraphs") <- rep(list(NULL), ncol(object)) %>% set_names(nm = colnames(object))
  return(object)
}

#' @rdname RemoveCellGraphs
#' @method RemoveCellGraphs PNAAssay5
#'
#' @export
#'
RemoveCellGraphs.PNAAssay5 <- RemoveCellGraphs.PNAAssay

#' @param assay The name of the target assay
#'
#' @rdname RemoveCellGraphs
#' @method RemoveCellGraphs Seurat
#'
#' @export
#'
RemoveCellGraphs.Seurat <- function(
  object,
  assay = NULL,
  ...
) {
  if (!is.null(assay)) {
    assert_single_value(assay, type = "string")
  } else {
    # Use default assay if assay = NULL
    assay <- DefaultAssay(object)
  }

  pixel_assay <- object[[assay]]
  assert_class(
    pixel_assay,
    classes = c("CellGraphAssay", "CellGraphAssay5", "PNAAssay", "PNAAssay5")
  )
  pixel_assay <- RemoveCellGraphs(pixel_assay)

  object[[assay]] <- pixel_assay
  return(object)
}
