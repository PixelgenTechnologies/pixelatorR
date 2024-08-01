#' @rdname RemoveCellGraphs
#' @method RemoveCellGraphs MPXAssay
#'
#' @export
#'
RemoveCellGraphs.MPXAssay <- function (
  object,
  ...
) {
  slot(object, name = "cellgraphs") <- rep(list(NULL), ncol(object)) %>% set_names(nm = colnames(object))
  return(object)
}

#' @param assay The name of the target assay
#'
#' @import rlang
#'
#' @rdname RemoveCellGraphs
#' @method RemoveCellGraphs Seurat
#'
#' @export
#'
RemoveCellGraphs.Seurat <- function (
  object,
  assay = NULL,
  ...
) {

  if (!is.null(assay)) {
    stopifnot("'assay' must be a character of length 1" = is.character(assay) & (length(assay) == 1))
  } else {
    # Use default assay if assay = NULL
    assay <- DefaultAssay(object)
  }

  cg_assay <- object[[assay]]
  if (!is(cg_assay, "MPXAssay")) {
    abort(glue("Invalid assay type '{class(cg_assay)}'. Expected a 'CellGraphAssay'",
               " or a 'CellGraphAssay5' object."))
  }
  cg_assay <- RemoveCellGraphs(cg_assay)

  object[[assay]] <- cg_assay
  return(object)
}
