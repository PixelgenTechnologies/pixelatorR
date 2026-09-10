#' @include CellGraph.R
NULL

#' The CellGraphList class
#'
#' A named list of \code{\link{CellGraph}} objects. Subsetting,
#' concatenation, and replacement type-check elements against
#' \code{CellGraph}. Unloaded graphs may be stored as \code{NULL};
#' \code{x[[i]] <- NULL} keeps the name and stores \code{NULL}
#' rather than dropping the element.
#' Use \code{lapply.CellGraphList} to apply a function without dropping
#' the \code{CellGraphList} class (\code{base::lapply} is not an S3 generic).
#'
#' @param cellgraphs A named list of \code{\link{CellGraph}} objects.
#' Unloaded graphs may be represented as \code{NULL}.
#' @param x,X A \code{\link{CellGraphList}} object
#' @param i Index to extract or replace
#' @param value A \code{\link{CellGraph}}, \code{NULL}, or a list of those
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
#' @export
#' @concept cellgraph
#'
CreateCellGraphList <- function(cellgraphs = list()) {
  .validate_cellgraph_list(cellgraphs)
  structure(cellgraphs, class = c("CellGraphList", "list"))
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
#' @method [ CellGraphList
#' @export
#'
`[.CellGraphList` <- function(x, i, ...) {
  CreateCellGraphList(NextMethod())
}

#' @rdname CellGraphList
#' @method [<- CellGraphList
#' @export
#'
`[<-.CellGraphList` <- function(x, i, value) {
  if (inherits(value, "CellGraphList")) {
    value <- as.list.CellGraphList(value)
  }
  x <- unclass(x)
  if (is.null(value)) {
    value <- vector("list", length(x[i]))
  }
  x[i] <- value
  CreateCellGraphList(x)
}

#' @rdname CellGraphList
#' @method [[<- CellGraphList
#' @export
#'
`[[<-.CellGraphList` <- function(x, i, value) {
  if (!(is.null(value) || inherits(value, "CellGraph"))) {
    cli::cli_abort(
      c("x" = "Replacement values must be {.cls CellGraph} objects or {.cls NULL}.")
    )
  }
  x <- unclass(x)
  if (is.null(value)) {
    x[i] <- list(NULL)
  } else {
    x[[i]] <- value
  }
  CreateCellGraphList(x)
}

#' @rdname CellGraphList
#' @method names<- CellGraphList
#' @export
#'
`names<-.CellGraphList` <- function(x, value) {
  x <- unclass(x)
  names(x) <- value
  CreateCellGraphList(x)
}

#' @rdname CellGraphList
#' @method c CellGraphList
#' @export
#'
c.CellGraphList <- function(...) {
  pieces <- lapply(list(...), function(elt) {
    if (is.null(elt)) {
      return(list())
    }
    if (inherits(elt, "CellGraphList")) {
      return(unclass(elt))
    }
    if (inherits(elt, "CellGraph")) {
      return(list(elt))
    }
    if (is.list(elt)) {
      return(elt)
    }
    cli::cli_abort(
      c("x" = "Can only concatenate {.cls CellGraphList} with lists of {.cls CellGraph} or {.cls NULL}.")
    )
  })
  CreateCellGraphList(do.call(c, pieces))
}

#' @rdname CellGraphList
#' @method as.list CellGraphList
#' @export
#'
as.list.CellGraphList <- function(x, ...) {
  unclass(x)
}

#' @rdname CellGraphList
#' @method lapply CellGraphList
#' @export
#'
lapply.CellGraphList <- function(X, FUN, ...) {
  CreateCellGraphList(lapply(as.list.CellGraphList(X), FUN, ...))
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
    return(unclass(x))
  }
  x
}
