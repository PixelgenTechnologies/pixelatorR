#' @include CellGraph.R
NULL

#' The CellGraphList class
#'
#' A named list of \code{\link{CellGraph}} objects. Unloaded graphs may be
#' stored as \code{NULL}. See \code{\link{CellGraphList-methods}} for
#' subsetting, replacement, and concatenation.
#'
#' @param cellgraphs A named list of \code{\link{CellGraph}} objects.
#' Unloaded graphs may be represented as \code{NULL}.
#'
#' @return A \code{CellGraphList} object
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
#' @seealso \code{\link{CellGraphList-methods}}
#' @name CellGraphList
#' @rdname CellGraphList
#' @export
#' @concept cellgraph
#'
CreateCellGraphList <- function(cellgraphs = list()) {
  .validate_cellgraph_list(cellgraphs)
  structure(cellgraphs, class = c("CellGraphList", "list"))
}


# -------------------------------------------------------
# Methods
# -------------------------------------------------------

#' CellGraphList Methods
#'
#' Methods for \code{\link{CellGraphList}} objects. Subsetting, concatenation,
#' and replacement type-check elements against \code{\link{CellGraph}}.
#' Unloaded graphs may be stored as \code{NULL}; \code{x[[i]] <- NULL} keeps
#' the name and stores \code{NULL} rather than dropping the element.
#'
#' @param x A \code{\link{CellGraphList}} object
#' @param i Index to extract or replace
#' @param value A \code{\link{CellGraph}}, \code{NULL}, or a list of those
#' @param ... Currently not used
#'
#' @return \code{[}, \code{[<-}, \code{[[<-}, \code{names<-}, and \code{c}:
#' a \code{CellGraphList}. \code{as.list}: a named list.
#' \code{print}: \code{x}, invisibly.
#'
#' @examples
#' library(pixelatorR)
#' library(tidygraph)
#' library(dplyr)
#'
#' edges <- tibble(from = c("a", "b"), to = c("b", "c"))
#' g <- as_tbl_graph(edges, directed = FALSE) %N>%
#'   mutate(node_type = c("umi1", "umi2", "umi1"))
#' attr(g, "type") <- "bipartite"
#' cg <- CreateCellGraphObject(cellgraph = g)
#' cgl <- CreateCellGraphList(list(cell_1 = cg, cell_2 = cg))
#'
#' # Print and subset
#' print(cgl)
#' cgl[1]
#'
#' # Unload a graph without dropping its name
#' cgl[[1]] <- NULL
#'
#' # Concatenate while keeping the class
#' c(cgl[1], cgl[2])
#'
#' @name CellGraphList-methods
#' @rdname CellGraphList-methods
#' @seealso \code{\link{CellGraphList}}
#' @concept cellgraph
#'
NULL

#' @describeIn CellGraphList-methods Print a \code{CellGraphList}
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

#' @describeIn CellGraphList-methods Subset a \code{CellGraphList}. Unknown
#' character names raise an error.
#' @method [ CellGraphList
#' @export
#'
`[.CellGraphList` <- function(x, i, ...) {
  if (!missing(i) && is.character(i)) {
    unknown <- unique(i[is.na(i) | !i %in% names(x)])
    if (length(unknown) > 0) {
      cli::cli_abort(
        c("x" = "Unknown name{?s} in {.cls CellGraphList}: {.val {unknown}}.")
      )
    }
  }
  CreateCellGraphList(NextMethod())
}

#' @describeIn CellGraphList-methods Replace a subset of graphs. \code{NULL}
#' unloads the selected cells without dropping their names.
#' @method [<- CellGraphList
#' @export
#'
`[<-.CellGraphList` <- function(x, i, value) {
  if (inherits(value, "CellGraphList")) {
    value <- as.list.CellGraphList(value)
  } else if (inherits(value, "CellGraph")) {
    value <- list(value)
  }
  x <- unclass(x)
  if (is.null(value)) {
    value <- vector("list", length(x[i]))
  }
  x[i] <- value
  CreateCellGraphList(x)
}

#' @describeIn CellGraphList-methods Replace a single graph. \code{NULL}
#' unloads that cell without dropping its name.
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

#' @describeIn CellGraphList-methods Set names. Names must be unique and
#' non-missing.
#' @method names<- CellGraphList
#' @export
#'
`names<-.CellGraphList` <- function(x, value) {
  x <- unclass(x)
  names(x) <- value
  CreateCellGraphList(x)
}

#' @describeIn CellGraphList-methods Concatenate \code{CellGraphList} objects
#' with lists of \code{CellGraph} or \code{NULL}. A bare \code{CellGraph}
#' uses the argument name, or \code{CellGraph1}, \code{CellGraph2}, ...
#' when unnamed.
#' @method c CellGraphList
#' @export
#'
c.CellGraphList <- function(...) {
  dots <- list(...)
  dot_names <- names(dots) %||% rep("", length(dots))
  pieces <- lapply(seq_along(dots), function(i) {
    elt <- dots[[i]]
    nm <- dot_names[[i]]
    if (is.null(elt)) {
      return(list())
    }
    if (inherits(elt, "CellGraphList")) {
      return(unclass(elt))
    }
    if (inherits(elt, "CellGraph")) {
      if (is.null(nm) || is.na(nm) || !nzchar(nm)) {
        nm <- paste0("CellGraph", i)
      }
      return(stats::setNames(list(elt), nm))
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

#' @describeIn CellGraphList-methods Convert to a named list
#' @method as.list CellGraphList
#' @export
#'
as.list.CellGraphList <- function(x, ...) {
  unclass(x)
}

#' @describeIn CellGraph-methods Pull node-level data from each loaded
#' \code{CellGraph} in a \code{CellGraphList}. Unlike
#' \code{\link{FetchLayoutData}}, this does not require a stored layout and
#' does not reserve coordinate names, so \code{vars} may include \code{x},
#' \code{y}, or \code{z} when those columns exist on the graphs. \code{component}
#' is reserved for the source graph ID. Node IDs are used as row names and
#' must be unique across the graphs being combined. Variables missing from a
#' graph are filled with \code{NA}. Variables missing from every graph are
#' omitted, with the same warning as \code{FetchData.CellGraph}. Graphs that
#' do not have a requested \code{layer} omit only features from that layer;
#' metadata, vertex attributes, reductions, and other layer values are kept.
#' \code{clean} defaults to \code{FALSE} so those missing values are kept.
#' \code{add_protein = TRUE} adds a \code{protein} column from the one-hot
#' counts matrix of each graph.
#' @method FetchData CellGraphList
#' @export
#'
FetchData.CellGraphList <- function(
  object,
  vars,
  cells = NULL,
  layer = NULL,
  clean = FALSE,
  add_protein = FALSE,
  ...
) {
  cells <- .resolve_loaded_cellgraph_ids(object, cells, fn = "FetchData")
  assert_single_value(add_protein, type = "bool")

  if (!is.null(vars) && length(vars) > 0) {
    vars <- as.character(vars)
    if ("component" %in% vars) {
      cli::cli_abort(
        c(
          "x" = "{.arg vars} cannot include {.val component}.",
          "i" = "{.val component} identifies the source {.cls CellGraph} in the result."
        )
      )
    }
    if (isTRUE(add_protein) && "protein" %in% vars) {
      cli::cli_abort(
        c(
          "x" = "{.arg vars} cannot include reserved column name {.val protein}.",
          "i" = "Protein labels are added with {.arg add_protein} = {.code TRUE}."
        )
      )
    }
  }
  if (isTRUE(clean)) {
    clean <- "all"
  } else if (isFALSE(clean)) {
    clean <- "none"
  }
  clean <- rlang::arg_match0(clean, values = c("all", "none"))

  pieces <- lapply(cells, function(nm) {
    cg <- object[[nm]]
    node_ids <- .cg_node_map(cg)
    df <- .fetch_layout_vars(
      object = cg,
      vars = vars,
      cells = node_ids,
      layer = layer,
      fill_missing = FALSE,
      missing_layer = "omit"
    )
    row_ids <- rownames(df)
    if (is.null(row_ids) || length(row_ids) != nrow(df)) {
      row_ids <- node_ids
    }
    list(nm = nm, df = df, row_ids = row_ids, node_ids = node_ids, cg = cg)
  })
  found <- unique(unlist(lapply(pieces, function(p) names(p$df)), use.names = FALSE))
  .warn_unfound_fetch_vars(vars, found)
  keep_vars <- if (is.null(vars) || length(vars) == 0) {
    found
  } else {
    intersect(vars, found)
  }

  frames <- .match_fill_classes(lapply(pieces, function(p) {
    n <- length(p$node_ids)
    df <- data.frame(
      component = rep_len(p$nm, n),
      stringsAsFactors = FALSE,
      check.names = FALSE,
      row.names = p$node_ids
    )
    if (isTRUE(add_protein)) {
      df$protein <- .node_protein_labels(p$cg, nodes = p$node_ids)
    }
    fetched <- p$df
    src_ids <- p$row_ids
    for (v in keep_vars) {
      if (v %in% names(fetched) && nrow(fetched) > 0) {
        df[[v]] <- fetched[[v]][match(p$node_ids, src_ids)]
      } else {
        df[[v]] <- NA
      }
    }
    df
  }))
  if (length(frames) > 1) {
    all_ids <- unlist(lapply(frames, rownames), use.names = FALSE)
    dup <- unique(all_ids[duplicated(all_ids)])
    if (length(dup) > 0) {
      cli::cli_abort(
        c(
          "x" = "Node IDs are duplicated across {.cls CellGraph} objects.",
          "i" = "Example: {.val {head(dup, 3)}}",
          "i" = "Each node ID can appear in only one component when combining with {.fn FetchData}."
        )
      )
    }
  }
  fetched <- do.call(rbind, frames)
  if (is.null(fetched)) {
    fetched <- data.frame(component = character(), stringsAsFactors = FALSE)
  }

  skip_clean <- "component"
  if (isTRUE(add_protein)) {
    skip_clean <- c(skip_clean, "protein")
  }
  value_cols <- setdiff(names(fetched), skip_clean)
  if (identical(clean, "all") && length(value_cols) > 0 && nrow(fetched) > 0) {
    no_data <- which(apply(fetched[value_cols], 1L, function(x) all(is.na(x))))
    if (length(no_data) > 0) {
      cli::cli_warn("Removing {length(no_data)} node{?s} missing data for vars requested")
      fetched <- fetched[-no_data, , drop = FALSE]
    }
  }
  fetched
}

#' Give NA fill columns the class used by the graphs that had the variable
#'
#' \code{\link{.fetch_layout_vars}} fills a variable that a graph does not
#' have with a bare logical \code{NA}. \code{rbind} then coerces the whole
#' column to that type, so a factor would come back as character and a
#' \code{POSIXct} as numeric. Filled columns are replaced with \code{NA} of
#' the class found on a graph that did have the variable. Factor columns
#' present on more than one graph get the union of their levels so
#' \code{rbind} does not coerce them to character.
#'
#' @param frames A list of \code{data.frame} objects with identical columns
#'
#' @return \code{frames} with fill columns matched to the other frames
#'
#' @keywords internal
#' @noRd
#'
.match_fill_classes <- function(frames) {
  if (length(frames) < 2) {
    return(frames)
  }
  is_fill <- function(x) is.logical(x) && !is.object(x) && all(is.na(x))
  for (v in names(frames[[1]])) {
    proto_at <- Position(function(df) !is_fill(df[[v]]), frames, nomatch = 0L)
    if (proto_at == 0L) {
      next
    }
    proto <- frames[[proto_at]][[v]]
    for (i in seq_along(frames)) {
      if (i != proto_at && is_fill(frames[[i]][[v]])) {
        frames[[i]][[v]] <- proto[rep(NA_integer_, nrow(frames[[i]]))]
      }
    }
    cols <- lapply(frames, `[[`, v)
    if (all(vapply(cols, is.factor, logical(1)))) {
      ordered <- any(vapply(cols, is.ordered, logical(1)))
      lvls <- unique(unlist(lapply(cols, levels), use.names = FALSE))
      for (i in seq_along(frames)) {
        frames[[i]][[v]] <- factor(
          as.character(frames[[i]][[v]]),
          levels = lvls,
          ordered = ordered
        )
      }
    }
  }
  frames
}

#' Test whether an object is a CellGraph or an unloaded placeholder
#'
#' \code{inherits(NULL, "NULL")} is not portable across R versions, so
#' unloaded graphs are detected with \code{is.null()}.
#'
#' @param x A candidate list element
#'
#' @return \code{TRUE} if \code{x} is \code{NULL} or a \code{CellGraph}
#'
#' @keywords internal
#' @noRd
#'
.is_cellgraph_or_null <- function(x) {
  is.null(x) || inherits(x, "CellGraph")
}

#' Validate a list of CellGraph objects
#'
#' Ensures every element is a \code{CellGraph} or \code{NULL}, and that
#' names are unique and non-missing. An empty list is allowed (no names
#' required). Errors are reported from \code{call}.
#'
#' @param cellgraphs A list
#' @param call Environment to report as the error caller
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
#' @noRd
#'
.validate_cellgraph_list <- function(cellgraphs, call = caller_env()) {
  assert_class(cellgraphs, "list", call = call)
  if (length(cellgraphs) == 0) {
    return(invisible(NULL))
  }
  is_ok <- vapply(cellgraphs, .is_cellgraph_or_null, logical(1))
  if (!all(is_ok)) {
    cli::cli_abort(
      c("x" = "All elements of {.arg cellgraphs} must be {.cls CellGraph} objects or {.cls NULL}."),
      call = call
    )
  }
  nm <- names(cellgraphs)
  if (is.null(nm) || any(is.na(nm) | nm == "")) {
    cli::cli_abort(
      c("x" = "The {.arg cellgraphs} list must have unique, non-missing names."),
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
#' Assay \code{cellgraphs} slots are plain lists. This strips the
#' \code{CellGraphList} class while keeping names and \code{NULL}
#' placeholders.
#'
#' @param x A \code{CellGraphList} or list
#'
#' @return A named list of \code{CellGraph} objects and \code{NULL}s
#'
#' @keywords internal
#' @noRd
#'
.unclass_cellgraph_list <- function(x) {
  if (inherits(x, "CellGraphList")) {
    return(unclass(x))
  }
  x
}
