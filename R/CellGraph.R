#' @include generics.R
#' @importFrom methods setClass setClassUnion setMethod slot slot<- new as slotNames
#' @importClassesFrom Matrix dgCMatrix
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definition
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The CellGraph class
#'
#' The CellGraph class is designed to hold information needed for working with
#' mpx single-cell graphs.
#'
#' @slot cellgraph A \code{tbl_graph} object corresponding to a cell graph
#' @slot counts A \code{matrix}-like object with marker counts
#' @slot layout A \code{list} of \code{tbl_df} objects with coordinates for cell layouts
#'
#' @name CellGraph-class
#' @rdname CellGraph-class
#' @exportClass CellGraph
#' @concept cellgraph
CellGraph <- setClass(
  Class = "CellGraph",
  slots = list(
    cellgraph = "ANY",
    counts = "ANY",
    layout = "ANY"
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a CellGraph object
#'
#' @param cellgraph A \code{tbl_graph} object representing an mpx single-cell graph
#' @param counts A \code{dgCMatrix} with marker counts
#' @param layout A \code{tbl_df} object with cell layout(s)
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @concept cellgraph
#'
#' @return A \code{CellGraph} object
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' edge_list <-
#' ReadMPX_item(
#'   system.file("extdata/five_cells",
#'               "five_cells.pxl",
#'               package = "pixelatorR"),
#'   items = "edgelist"
#' )
#' bipart_graph <-
#'   edge_list %>%
#'   select(upia, upib, marker) %>%
#'   distinct() %>%
#'   as_tbl_graph(directed = FALSE) %>%
#'   mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
#' attr(bipart_graph, "type") <- "bipartite"
#'
#' cg <- CreateCellGraphObject(cellgraph = bipart_graph)
#' cg
#'
#' @export
#'
CreateCellGraphObject <- function (
    cellgraph,
    counts = NULL,
    layout = NULL,
    verbose = FALSE
) {

  # Check input parameters
  stopifnot(
    "'cellgraph' must be a non-empty 'tbl_graph' object" =
      inherits(cellgraph, what = "tbl_graph") &&
      (length(cellgraph) > 0)
  )
  if (!is.null(counts)) {
    stopifnot(
      "'counts' must be a non-empty 'dgCMatrix' object" =
        inherits(counts, what = "dgCMatrix") &&
        (length(counts) > 0)
    )
  }
  if (!is.null(layout)) {
    stopifnot(
      "'layout' must be a non-empty 'tbl_df' object" =
        inherits(layout, what = "tbl_df") &&
        (length(layout) > 0)
    )
  }

  if (!"type" %in% names(attributes(cellgraph))) {
    abort("Graph attribute 'type' is missing.")
  } else {
    if (verbose && check_global_verbosity())
      cli_alert_info("Got a graph of type '{attr(cellgraph, 'type')}'")
  }

  # Add checks for graph types
  if (attr(cellgraph, "type") == "bipartite") {
    stopifnot(
      "Node attribute 'name' is missing" =
        "name" %in% vertex_attr_names(cellgraph),
      "Node attribute 'node_type' is missing" =
        "node_type" %in% vertex_attr_names(cellgraph)
    )
  }
  #TODO: Add check for A-node-projection and linegraph

  # create object
  cellgraph <- new(
    Class = "CellGraph",
    cellgraph = cellgraph,
    counts = counts,
    layout = layout
  )

  return(cellgraph)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get and set CellGraph object data
#'
#' @param object A \code{\link{CellGraph}} object
#' @param slot Information to pull from object (cellgraph, meta_data, layout)
#'
#' @import rlang
#'
#' @rdname CellGraphData
#'
#' @return \code{GetCellGraphData}: A \code{\link{CellGraph}} object slot
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' edge_list <-
#'   ReadMPX_item(
#'     system.file("extdata/five_cells",
#'                 "five_cells.pxl",
#'                 package = "pixelatorR"),
#'     items = "edgelist"
#'   )
#' bipart_graph <-
#'   edge_list %>%
#'   select(upia, upib, marker) %>%
#'   distinct() %>%
#'   as_tbl_graph(directed = FALSE) %>%
#'   mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
#' attr(bipart_graph, "type") <- "bipartite"
#'
#' cg <- CreateCellGraphObject(cellgraph = bipart_graph)
#'
#' # Get slot data
#' CellGraphData(cg, slot = "cellgraph")
#'
#' @export
#'
CellGraphData <- function (
    object,
    slot = "cellgraph"
) {
  if (!inherits(object, what = "CellGraph")) abort(glue("Invalid class {class(object)}"))
  if (!(slot %in% slotNames(x = object))) {
    abort(glue("slot must be one of {paste(slotNames(x = object), collapse = ', ')}"))
  }
  return(slot(object = object, name = slot))
}


#' @param value A new variable to place in \code{slot}
#'
#' @rdname CellGraphData
#'
#' @return \code{CellGraphData<-}: A \code{\link{CellGraph}} with updated data
#'
#' @examples
#' # Set slot data
#' CellGraphData(cg, slot = "cellgraph") <- CellGraphData(cg, slot = "cellgraph")
#'
#'
#' @export
#'
"CellGraphData<-" <- function (
    object,
    slot = "cellgraph",
    value
) {
  if (!inherits(object, what = "CellGraph")) abort(glue("Invalid class {class(object)}"))
  if (!(slot %in% slotNames(x = object))) {
    abort(glue("slot must be one of {paste(slotNames(x = object), collapse = ', ')}"))
  }

  # Get counts and layouts
  cellgraph <- slot(object, name = "cellgraph")
  counts <- slot(object, name = "counts")
  layouts <- slot(object, name = "layout")

  # Validate cellgraph input
  if (slot == "cellgraph") {
    stopifnot(
      "'value' must be a 'tbl_graph' object" =
        inherits(value, what = "tbl_graph")
    )
    if (length(counts) > 0) {
      stopifnot(
        "names in 'value' do not match the row names of 'counts'" =
          all(names(value) == rownames(counts))
      )
    }
    if (length(layouts) > 0) {
      for (layout in names(layouts)) {
        if (length(value) != nrow(layouts[[layout]])) {
          abort(glue("Number of nodes in 'value' does not match the ",
                     "number of rows in the '{layout}' layout table"))
        }
      }
    }
    slot(object, name = "cellgraph") <- cellgraph
  }

  # Validate counts input
  if (slot == "counts") {
    stopifnot(
      "'value' must be a 'dgCMatrix' object" =
        inherits(value, what = "dgCMatrix")
    )
    if (length(layouts) > 0) {
      for (layout in names(layouts)) {
        if (nrow(counts) != nrow(layouts[[layout]])) {
          abort(glue("Number of rows in 'value' does not match the ",
                     "number of rows in the '{layout}' layout table"))
        }
      }
    }
    stopifnot(
      "Number of rows in 'value' do not match the number of nodes in 'cellgraph'" =
        length(cellgraph) == nrow(value)
    )
    slot(object, name = "counts") <- counts
  }

  # Validate layout input
  if (slot == "layout") {
    stopifnot(
      "'value' must be a 'list'" =
        inherits(value, what = "list"),
      "'value' is empty" =
        length(value) > 0,
      "'names' of 'value' cannot be NULL" =
        !is.null(names(value))
    )
    for (layout in names(value)) {
      stopifnot(
        inherits(value[[layout]], what = "tbl_df")
      )
      if (length(cellgraph) != nrow(value[[layout]])) {
        abort(glue("Number of nodes in 'cellgraph' does not match the ",
                   "number of rows in the '{layout}' layout table provided in 'value'"))
      }
      if (length(counts) > 0) {
        if (nrow(counts) != nrow(value[[layout]])) {
          abort(glue("Number of rows in 'counts' does not match ",
                     "the number of rows in '{layout}' table provided in 'value'"))
        }
      }
    }
    slot(object, name = "layout") <- value
  }

  return(object)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Base methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Show method for \code{CellGraph} object
#'
#' @param object A \code{CellGraph} object
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' edge_list <-
#' ReadMPX_item(
#'   system.file("extdata/five_cells",
#'               "five_cells.pxl",
#'               package = "pixelatorR"),
#'   items = "edgelist"
#' )
#' bipart_graph <-
#'   edge_list %>%
#'   select(upia, upib, marker) %>%
#'   distinct() %>%
#'   as_tbl_graph(directed = FALSE) %>%
#'   mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
#' attr(bipart_graph, "type") <- "bipartite"
#'
#' cg <- CreateCellGraphObject(cellgraph = bipart_graph)
#'
#' # Show method
#' cg
#'
setMethod (
  f = "show",
  signature = "CellGraph",
  definition = function(object) {
    graph_type <- attr(slot(object, "cellgraph"), "type")
    if (is.null(slot(object, "counts"))) {
      n_markers <- NULL
    } else {
      n_markers <- ncol(slot(object, "counts"))
    }
    cat(
      "A CellGraph object containing a", graph_type, "graph with",
      slot(object = object, name = "cellgraph") %>% length() %>% col_br_blue(),
      "nodes and",
      slot(object = object, name = "cellgraph") %>% gsize() %>% col_br_blue(),
      "edges "
    )
    if (is.null(n_markers)) {
      cat("\n")
    } else {
      cat("covering", n_markers, "markers\n")
    }
  }
)
