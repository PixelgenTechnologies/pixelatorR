#' @include CellGraph.R
NULL

# -------------------------------------------------------
# Class definition
# -------------------------------------------------------

#' The CellGraphList class
#'
#' A simple container for a named list of \code{\link{CellGraph}} objects.
#'
#' @slot cellgraphs A named list of \code{\link{CellGraph}} objects
#'
#' @name CellGraphList-class
#' @rdname CellGraphList-class
#' @exportClass CellGraphList
#' @concept cellgraph
CellGraphList <- setClass(
  Class = "CellGraphList",
  slots = list(
    cellgraphs = "list"
  ),
  prototype = list(
    cellgraphs = list()
  )
)

# -------------------------------------------------------
# Create methods
# -------------------------------------------------------

#' Create a CellGraphList object
#'
#' @param cellgraphs A named list of \code{\link{CellGraph}} objects
#'
#' @return A \code{\link{CellGraphList}} object
#'
#' @examples
#' library(pixelatorR)
#'
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#' se <- LoadCellGraphs(se, cells = colnames(se)[1:2], verbose = FALSE)
#' cgl <- CreateCellGraphList(CellGraphs(se)[1:2])
#' cgl
#'
#' @export
#' @concept cellgraph
#'
CreateCellGraphList <- function(cellgraphs) {
  .validate_cellgraph_list(cellgraphs)
  new(Class = "CellGraphList", cellgraphs = cellgraphs)
}

# -------------------------------------------------------
# Methods
# -------------------------------------------------------

#' CellGraphList Methods
#'
#' Methods for \code{\link{CellGraphList}} objects
#'
#' @param object A \code{\link{CellGraphList}} object
#' @param x A \code{\link{CellGraphList}} object
#' @param i Index
#' @param j Unused
#' @param value Replacement value
#' @param drop Unused
#' @param ... Currently not used
#'
#' @name CellGraphList-methods
#' @rdname CellGraphList-methods
#'
#' @concept cellgraph
#'
NULL

#' Show method for \code{CellGraphList} object
#'
#' @describeIn CellGraphList-methods Show a \code{CellGraphList} object
#' @method show CellGraphList
#' @docType methods
#'
setMethod(
  f = "show",
  signature = "CellGraphList",
  definition = function(object) {
    n <- length(slot(object, "cellgraphs"))
    cat(
      "A CellGraphList with", col_br_blue(n),
      "CellGraph objects\n"
    )
    nm <- names(object)
    if (!is.null(nm) && length(nm) > 0) {
      shown <- nm[seq_len(min(length(nm), 5))]
      extra <- if (length(nm) > 5) ", ..." else ""
      cat("Names:", col_br_blue(paste(shown, collapse = ", ")), extra, "\n")
    }
  }
)

#' @describeIn CellGraphList-methods Number of \code{CellGraph} objects
#' @export
#'
setMethod(
  f = "length",
  signature = "CellGraphList",
  definition = function(x) {
    length(slot(x, "cellgraphs"))
  }
)

#' @describeIn CellGraphList-methods Names of \code{CellGraph} objects
#' @export
#'
setMethod(
  f = "names",
  signature = "CellGraphList",
  definition = function(x) {
    names(slot(x, "cellgraphs"))
  }
)

#' @describeIn CellGraphList-methods Set names of \code{CellGraph} objects
#' @export
#'
setMethod(
  f = "names<-",
  signature = c(x = "CellGraphList", value = "ANY"),
  definition = function(x, value) {
    names(slot(x, "cellgraphs")) <- value
    .validate_cellgraph_list(slot(x, "cellgraphs"))
    x
  }
)

#' @describeIn CellGraphList-methods Extract a \code{CellGraph}
#' @export
#'
setMethod(
  f = "[[",
  signature = c(x = "CellGraphList", i = "ANY", j = "missing"),
  definition = function(x, i, j, ..., drop = TRUE) {
    slot(x, "cellgraphs")[[i]]
  }
)

#' @describeIn CellGraphList-methods Replace a \code{CellGraph}
#' @export
#'
setMethod(
  f = "[[<-",
  signature = c(x = "CellGraphList", i = "ANY", j = "missing", value = "ANY"),
  definition = function(x, i, j, ..., value) {
    slot(x, "cellgraphs")[[i]] <- value
    .validate_cellgraph_list(slot(x, "cellgraphs"))
    x
  }
)

#' @describeIn CellGraphList-methods Subset a \code{CellGraphList}
#' @export
#'
setMethod(
  f = "[",
  signature = c(x = "CellGraphList", i = "ANY", j = "missing", drop = "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    CreateCellGraphList(slot(x, "cellgraphs")[i])
  }
)

#' @describeIn CellGraphList-methods Convert to a list of \code{CellGraph} objects
#' @method as.list CellGraphList
#' @export
#'
as.list.CellGraphList <- function(x, ...) {
  slot(x, "cellgraphs")
}

# -------------------------------------------------------
# Internal helpers
# -------------------------------------------------------

#' Validate a list of CellGraph objects
#'
#' @noRd
#'
.validate_cellgraph_list <- function(cellgraphs, call = caller_env()) {
  assert_class(cellgraphs, "list", call = call)
  if (length(cellgraphs) == 0) {
    return(invisible(NULL))
  }
  is_cg <- vapply(cellgraphs, function(x) is(x, "CellGraph"), logical(1))
  if (!all(is_cg)) {
    cli::cli_abort(
      c("x" = "All elements of {.arg cellgraphs} must be {.cls CellGraph} objects."),
      call = call
    )
  }
  nm <- names(cellgraphs)
  if (is.null(nm) || any(nm == "")) {
    cli::cli_abort(
      c("x" = "The {.arg cellgraphs} list must be named."),
      call = call
    )
  }
  if (anyDuplicated(nm)) {
    cli::cli_abort(
      c("x" = "The {.arg cellgraphs} list must have unique names."),
      call = call
    )
  }
  invisible(NULL)
}
