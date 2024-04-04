#' @rdname RemoveCellGraphs
#' @method RemoveCellGraphs CellGraphAssay
#'
#' @export
#'
RemoveCellGraphs.CellGraphAssay <- function (
  object,
  ...
){
  slot(object, name = "cellgraphs") <- rep(list(NULL), ncol(object)) %>% setNames(nm = colnames(object))
  return(object)
}

#' @rdname RemoveCellGraphs
#' @method RemoveCellGraphs CellGraphAssay5
#'
#' @export
#'
RemoveCellGraphs.CellGraphAssay5 <- function (
  object,
  ...
){
  slot(object, name = "cellgraphs") <- rep(list(NULL), ncol(object)) %>% setNames(nm = colnames(object))
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
){

  if (!is.null(assay)) {
    stopifnot("'assay' must be a character of length 1" = is.character(assay) & (length(assay) == 1))
  } else {
    # Use default assay if assay = NULL
    assay <- DefaultAssay(object)
  }

  cg_assay <- object[[assay]]
  if (!inherits(cg_assay, what = "CellGraphAssay")) {
    abort(glue("Invalid assay type '{class(cg_assay)}'. Expected a 'CellGraphAssay'"))
  }
  cg_assay <- RemoveCellGraphs(cg_assay)

  object[[assay]] <- cg_assay
  return(object)
}

