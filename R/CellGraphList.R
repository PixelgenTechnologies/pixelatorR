#' @include CellGraph.R
NULL

#' The CellGraphList class
#'
#' A named list of \code{\link{CellGraph}} objects. \code{CellGraphList}
#' extends \code{list}, so list operations such as \code{[[}, \code{lapply},
#' \code{names}, and \code{length} work as usual. The only specialized
#' behavior is printing, which shows a short summary instead of each
#' \code{CellGraph}.
#'
#' @param cellgraphs A named list of \code{\link{CellGraph}} objects
#' @param x A \code{\link{CellGraphList}} object
#' @param ... Currently not used
#'
#' @return \code{CreateCellGraphList}: a \code{CellGraphList} object
#'
#' @examples
#' library(pixelatorR)
#'
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#' se <- LoadCellGraphs(se, cells = colnames(se)[1:2], verbose = FALSE)
#' cgl <- CreateCellGraphList(CellGraphs(se)[1:2])
#' cgl
#'
#' @name CellGraphList
#' @rdname CellGraphList
#' @export
#' @concept cellgraph
#'
CreateCellGraphList <- function(cellgraphs) {
  .validate_cellgraph_list(cellgraphs)
  structure(cellgraphs, class = c("CellGraphList", "list"))
}

#' @rdname CellGraphList
#' @method print CellGraphList
#' @export
#'
print.CellGraphList <- function(x, ...) {
  n <- length(x)
  cat(
    "A CellGraphList with", col_br_blue(n),
    "CellGraph objects\n"
  )
  nm <- names(x)
  if (!is.null(nm) && length(nm) > 0) {
    shown <- nm[seq_len(min(length(nm), 5))]
    extra <- if (length(nm) > 5) ", ..." else ""
    cat("Names:", col_br_blue(paste(shown, collapse = ", ")), extra, "\n")
  }
  invisible(x)
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
