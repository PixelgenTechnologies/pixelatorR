#' @include CellGraph.R
NULL

#' Treat CellGraph objects as length-1 vctrs
#'
#' \code{vctrs::new_list_of()} requires a prototype that is a true vector.
#' An empty S3 \code{CellGraph} list is used as that prototype; S4
#' \code{CellGraph} objects are proxied as length-1 lists so they can be
#' stored and type-checked as elements of a \code{CellGraphList}.
#'
#' @noRd
#' @export
#' @importFrom vctrs vec_proxy
#' @method vec_proxy CellGraph
#'
vec_proxy.CellGraph <- function(x, ...) {
  if (isS4(x)) {
    list(.cellgraph = x)
  } else {
    unclass(x)
  }
}

#' @noRd
#' @export
#' @importFrom vctrs vec_restore
#' @method vec_restore CellGraph
#'
vec_restore.CellGraph <- function(x, to, ...) {
  if (length(x) == 0) {
    return(structure(list(), class = "CellGraph"))
  }
  if (is.list(x) && !is.null(x$.cellgraph)) {
    return(x$.cellgraph)
  }
  x
}

#' The CellGraphList class
#'
#' A named list of \code{\link{CellGraph}} objects. \code{CellGraphList}
#' is a \code{\link[vctrs:list_of]{vctrs} list_of} subclass, so subsetting,
#' concatenation, and replacement type-check elements against
#' \code{CellGraph}. Unloaded graphs may be stored as \code{NULL}.
#' Use \code{lapply.CellGraphList} to apply a function without dropping
#' the \code{CellGraphList} class (\code{base::lapply} is not an S3 generic).
#'
#' @param cellgraphs A named list of \code{\link{CellGraph}} objects.
#' Unloaded graphs may be represented as \code{NULL}.
#' @param x,X A \code{\link{CellGraphList}} object
#' @param FUN A function to apply to each element
#' @param ... Currently not used
#'
#' @return \code{CreateCellGraphList}: a \code{CellGraphList} object
#'
#' @examples
#' library(pixelatorR)
#' library(tidygraph)
#' library(dplyr)
#'
#' # Build a small dummy cell graph
#' edges <- tibble(from = c("a", "b"), to = c("b", "c"))
#' g <- as_tbl_graph(edges, directed = FALSE) %N>%
#'   mutate(node_type = c("umi1", "umi2", "umi1"))
#' attr(g, "type") <- "bipartite"
#' cg <- CreateCellGraphObject(cellgraph = g)
#'
#' # Repeat the CellGraph in a named list and convert to a CellGraphList
#' cgl <- CreateCellGraphList(list(cell_1 = cg, cell_2 = cg))
#' cgl
#'
#' @name CellGraphList
#' @rdname CellGraphList
#' @importFrom vctrs new_list_of
#' @export
#' @concept cellgraph
#'
CreateCellGraphList <- function(cellgraphs = list()) {
  .validate_cellgraph_list(cellgraphs)
  vctrs::new_list_of(
    x = cellgraphs,
    ptype = structure(list(), class = "CellGraph"),
    class = "CellGraphList"
  )
}

#' @rdname CellGraphList
#' @method print CellGraphList
#' @export
#'
print.CellGraphList <- function(x, ...) {
  n <- length(x)
  n_loaded <- sum(vapply(x, function(el) inherits(el, "CellGraph"), logical(1)))
  cat(
    "A CellGraphList with", col_br_blue(n_loaded),
    "loaded CellGraph object(s) out of", col_br_blue(n), "\n"
  )
  nm <- names(x)
  if (!is.null(nm) && length(nm) > 0) {
    shown <- nm[seq_len(min(length(nm), 5))]
    extra <- if (length(nm) > 5) ", ..." else ""
    cat("Names:", col_br_blue(paste(shown, collapse = ", ")), extra, "\n")
  }
  invisible(x)
}

#' @rdname CellGraphList
#' @method lapply CellGraphList
#' @export
#'
lapply.CellGraphList <- function(X, FUN, ...) {
  CreateCellGraphList(lapply(as.list(X), FUN, ...))
}

#' Validate a list of CellGraph objects
#'
#' @noRd
#'
.validate_cellgraph_list <- function(cellgraphs, call = caller_env()) {
  assert_class(cellgraphs, "list", call = call)
  if (length(cellgraphs) == 0) {
    return(invisible(NULL))
  }
  is_ok <- vapply(cellgraphs, function(x) {
    is.null(x) || inherits(x, "CellGraph")
  }, logical(1))
  if (!all(is_ok)) {
    cli::cli_abort(
      c("x" = "All elements of {.arg cellgraphs} must be {.cls CellGraph} objects or {.cls NULL}."),
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

#' Drop the CellGraphList class so the object can be stored in an S4 list slot
#'
#' @noRd
#'
.unclass_cellgraph_list <- function(x) {
  if (inherits(x, "CellGraphList")) {
    return(as.list(x))
  }
  x
}
