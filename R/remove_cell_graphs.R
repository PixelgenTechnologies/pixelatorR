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

#' @param assay The name of the target assay
#'
#' @import rlang
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

  cg_assay <- object[[assay]]
  assert_mpx_assay(cg_assay)
  cg_assay <- RemoveCellGraphs(cg_assay)

  object[[assay]] <- cg_assay
  return(object)
}
