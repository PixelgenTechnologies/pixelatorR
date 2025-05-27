#' @include generics.R
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
#'   ReadMPX_item(
#'     minimal_mpx_pxl_file(),
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
#' cg
#'
#' @export
#'
CreateCellGraphObject <- function(
  cellgraph,
  counts = NULL,
  layout = NULL,
  verbose = FALSE
) {
  # Validate input parameters
  assert_non_empty_object(cellgraph, classes = "tbl_graph")
  assert_non_empty_object(counts, classes = "dgCMatrix", allow_null = TRUE)
  assert_non_empty_object(layout, classes = "tbl_df", allow_null = TRUE)

  if (!"type" %in% names(attributes(cellgraph))) {
    cli::cli_abort(c("x" = "Graph attribute {.str type} is missing."))
  } else {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("Got a graph of type '{attr(cellgraph, 'type')}'")
    }
  }

  # Add checks for graph types
  if (attr(cellgraph, "type") == "bipartite") {
    if (!"name" %in% vertex_attr_names(cellgraph)) {
      cli::cli_abort(
        "x" = "Node attribute {.str name} is missing from the graph"
      )
    }
    if (!"node_type" %in% vertex_attr_names(cellgraph)) {
      cli::cli_abort(
        "x" = "Node attribute {.str node_type} is missing from the graph"
      )
    }
  }
  # TODO: Add check for A-node-projection and linegraph

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
#'     minimal_mpx_pxl_file(),
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
CellGraphData <- function(
  object,
  slot = "cellgraph"
) {
  assert_class(object, "CellGraph")
  assert_single_value(slot, type = "string")
  assert_is_one_of(slot, slotNames(x = object))
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
#' @export
#'
"CellGraphData<-" <- function(
  object,
  slot = "cellgraph",
  value
) {
  assert_class(object, "CellGraph")
  assert_is_one_of(slot, slotNames(x = object))

  # Get counts and layouts
  cellgraph <- slot(object, name = "cellgraph")
  counts <- slot(object, name = "counts")
  layouts <- slot(object, name = "layout")

  # Validate cellgraph input
  if (slot == "cellgraph") {
    assert_class(value, "tbl_graph")
    if (length(counts) > 0) {
      assert_vectors_match(value %N>% pull(name), rownames(counts))
    }
    if (length(layouts) > 0) {
      for (layout in names(layouts)) {
        if (length(value) != nrow(layouts[[layout]])) {
          cli::cli_abort(
            c(
              "x" =
                "Number of nodes in the provided {.cls {class(value)}} does not match ",
              " " = "the number of rows in the 'layout' slot '{layout}' table"
            )
          )
        }
      }
    }
    slot(object, name = "cellgraph") <- cellgraph
  }

  # Validate counts input
  if (slot == "counts") {
    assert_class(value, "dgCMatrix")
    assert_vectors_match(cellgraph %N>% pull(name), rownames(value))
    if (length(layouts) > 0) {
      for (layout in names(layouts)) {
        if (nrow(counts) != nrow(layouts[[layout]])) {
          cli::cli_abort(
            c(
              "x" = "Number of rows ({nrow(counts)}) in the provided {.cls {class(value)}} does not match ",
              " " = "the number of rows ({nrow(layouts[[layout]])}) in the 'layout' slot '{layout}' table"
            )
          )
        }
      }
    }
    slot(object, name = "counts") <- counts
  }

  # Validate layout input
  if (slot == "layout") {
    assert_non_empty_object(value, "list")
    if (is.null(names(value))) {
      cli::cli_abort("The {.cls {class(value)}} must be named")
    }
    for (layout in names(value)) {
      if (!inherits(value[[layout]], what = "tbl_df")) {
        cli::cli_abort(
          c(
            "x" =
              "The '{layout}' layout table must be a {.cls tbl_df}"
          )
        )
      }
      if (length(cellgraph) != nrow(value[[layout]])) {
        cli::cli_abort(
          c(
            "x" = "Number of nodes ({length(cellgraph)}) in the 'cellgraph' slot does not match ",
            " " = "the number of rows ({nrow(value[[layout]])}) in the provided '{layout}' layout table"
          )
        )
      }
      if (length(counts) > 0) {
        if (nrow(counts) != nrow(value[[layout]])) {
          cli::cli_abort(
            c(
              "x" = "Number of rows ({nrow(counts)}) in the 'counts' slot does not match the ",
              " " = "number of rows ({nrow(value[[layout]])}) in the provided '{layout}' layout table"
            )
          )
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

#' CellGraph Methods
#'
#' Methods for \code{\link{CellGraph}} objects for generics defined in other
#' packages
#'
#' @param object A \code{\link{CellGraph}} object
#' @param x A \code{\link{CellGraph}} object
#' @param nodes A character vector of node names
#' @param ... Currently not used
#'
#' @name CellGraph-methods
#' @rdname CellGraph-methods
#'
#' @concept assay
#'
NULL

#' Show method for \code{CellGraph} object
#'
#' @describeIn CellGraph-methods Show a \code{CellGraph} object
#' @method show CellGraph
#' @docType methods
#'
#' @examples
#'
#' library(pixelatorR)
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#' se <- LoadCellGraphs(se, cells = colnames(se)[1], verbose = FALSE)
#' cg <- CellGraphs(se)[[1]]
#'
#' # Show method
#' cg
#'
setMethod(
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
      "A CellGraph object containing a", col_br_blue(graph_type), "graph with",
      slot(object = object, name = "cellgraph") %>% length() %>% col_br_blue(),
      "nodes and",
      slot(object = object, name = "cellgraph") %>% gsize() %>% col_br_blue(),
      "edges"
    )
    if (is.null(n_markers)) {
      cat("\n")
    } else {
      cat("\nNumber of markers: ", col_br_blue(n_markers), "\n")
    }
    if (!is.null(slot(object, "layout"))) {
      cat("Layouts:", col_br_blue(paste(names(slot(object, "layout"))), collapse = ", "), "\n")
    }
  }
)


#' subset method for \code{CellGraph} object
#'
#' @describeIn CellGraph-methods Subset a \code{CellGraph} object
#' @method subset CellGraph
#' @docType methods
#'
#' @examples
#' # Subset
#' cg_small <- subset(cg, nodes = rownames(cg@counts)[1:100])
#' cg_small
#'
#' @return A \code{CellGraph} object containing only the specified nodes.
#'
#' @export
#'
subset.CellGraph <- function(
  x,
  nodes,
  ...
) {
  assert_vector(nodes, type = "character")
  available_nodes <- x@cellgraph %N>% pull(name)
  assert_x_in_y(nodes, available_nodes)
  if (!is.null(x@counts)) {
    x@counts <- x@counts[match(nodes, available_nodes), ]
  }
  if (!is.null(x@layout)) {
    for (lyt in names(x@layout)) {
      x@layout[[lyt]] <- x@layout[[lyt]][match(nodes, available_nodes), ]
    }
  }
  x@cellgraph <- x@cellgraph %N>% filter(name %in% nodes)
  return(x)
}
