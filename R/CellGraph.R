#' @include NodeDimReduc.R
#' @importClassesFrom Matrix dgCMatrix
NULL

# -------------------------------------------------------
# Class definition
# -------------------------------------------------------

#' The CellGraph class
#'
#' The CellGraph class is designed to hold information needed for working with
#' PNA single-cell graphs.
#'
#' Node-level data are always aligned by **node name**, not by row position.
#' When counts, layouts, layers, metadata, or reductions are added, rows are
#' matched to the \code{name} vertex attribute of \code{cellgraph} and stored
#' in that node order. Shuffling the graph or a matrix no longer silently
#' breaks the mapping.
#'
#' Layouts are stored as \code{data.frame} objects with node IDs as row names.
#' A \code{name} column is accepted on input and is converted to row names, so
#' that stored layouts keep only their coordinate columns (typically \code{x},
#' \code{y}, \code{z}). Layouts without node IDs still work if the number of
#' rows matches the graph, which preserves objects created with earlier
#' versions of pixelatorR.
#'
#' @slot cellgraph A \code{tbl_graph} object corresponding to a cell graph
#' @slot counts A \code{matrix}-like object with marker counts (nodes x markers).
#' Row names are node names.
#' @slot layout A named \code{list} of \code{data.frame} objects with coordinates
#' for cell layouts. Row names are node names and the row order matches the
#' graph node order.
#' @slot layers A named \code{list} of additional numeric node matrices
#' (nodes x features), analogous to layers on a Seurat
#' \code{\link[SeuratObject]{Assay5}}. The counts matrix is not stored here;
#' it remains in \code{counts} and is exposed as the \code{"counts"} layer
#' via \code{\link[SeuratObject]{Layers}} / \code{\link[SeuratObject]{LayerData}}.
#' @slot meta.data A \code{data.frame} of node-level metadata (one row per node).
#' Row names are node names. Columns may have mixed types.
#' @slot reductions A named \code{list} of \code{\link{NodeDimReduc}} objects
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
    layout = "ANY",
    layers = "list",
    meta.data = "data.frame",
    reductions = "list"
  ),
  prototype = list(
    cellgraph = NULL,
    counts = NULL,
    layout = NULL,
    layers = list(),
    meta.data = data.frame(),
    reductions = list()
  )
)


# -------------------------------------------------------
# Create methods
# -------------------------------------------------------

#' Create a CellGraph object
#'
#' @param cellgraph A \code{tbl_graph} object representing a PNA single-cell graph
#' @param counts A \code{dgCMatrix} with marker counts. Rows are matched to graph
#' node names (order does not need to match).
#' @param layout A named \code{list} of \code{data.frame} objects with cell
#' layouts. Nodes are identified by row names or by a \code{name} column;
#' otherwise the row order is assumed to follow the graph (legacy behavior).
#' @param layers A named \code{list} of additional numeric node matrices
#' (nodes x features). \code{"counts"} is reserved.
#' @param meta.data A node-level \code{data.frame} or \code{tbl_df}. Either row
#' names or a \code{name} column must identify nodes.
#' @param reductions A named \code{list} of \code{\link{NodeDimReduc}} objects
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
#' # Open a database connection (PXL file)
#' db <- PixelDB$new(minimal_pna_pxl_file())
#'
#' # Select a component ID and load the edgelist
#' sel_comp <- db$cell_meta() %>%
#'   rownames() %>%
#'   head(1)
#' component_edgelist <- db$components_edgelist(
#'   components = sel_comp,
#'   umi_data_type = "suffixed_string"
#' ) %>%
#'   select(umi1, umi2)
#'
#' # Define node types for the bipartite graph
#' umi_node_type <- bind_rows(
#'   component_edgelist %>% select(name = umi1) %>% mutate(node_type = "umi1"),
#'   component_edgelist %>% select(name = umi2) %>% mutate(node_type = "umi2")
#' ) %>%
#'   distinct()
#'
#' # Create a bipartite graph from the edgelist and add node types
#' component_graph <- as_tbl_graph(component_edgelist, directed = FALSE) %N>%
#'   left_join(umi_node_type, by = "name")
#'
#' # Set the graph type attribute to "bipartite"
#' attr(component_graph, "type") <- "bipartite"
#'
#' # Create a CellGraph object with just the graph
#' cg <- CreateCellGraphObject(cellgraph = component_graph)
#' cg
#'
#' # Load cell count matrix
#' counts <- db$components_marker_counts(
#'   components = sel_comp, as_sparse = TRUE
#' )[[1]]
#'
#' # Counts are matched by row name; pre-ordering is optional
#' cg <- CreateCellGraphObject(cellgraph = component_graph, counts = counts)
#' cg
#'
#' # Create a CellGraph object with counts and layout
#' layout <- db$components_layout(
#'   components = sel_comp
#' )[[1]]
#'
#' # Layouts with a name column or node row names are matched automatically
#' cg <- CreateCellGraphObject(
#'   cellgraph = component_graph,
#'   counts = counts,
#'   layout = list(wpmds_3d = layout)
#' )
#' cg
#'
#' @export
#'
CreateCellGraphObject <- function(
  cellgraph,
  counts = NULL,
  layout = NULL,
  layers = NULL,
  meta.data = NULL,
  reductions = NULL,
  verbose = FALSE
) {
  # Validate input parameters
  assert_non_empty_object(cellgraph, classes = "tbl_graph")
  assert_non_empty_object(counts, classes = "dgCMatrix", allow_null = TRUE)
  assert_non_empty_object(layout, classes = "list", allow_null = TRUE)
  assert_class(layers, classes = "list", allow_null = TRUE)
  assert_class(reductions, classes = "list", allow_null = TRUE)

  .validate_cellgraph(cellgraph, verbose = verbose)

  node_names <- .cg_node_names(cellgraph)

  new(
    Class = "CellGraph",
    cellgraph = cellgraph,
    counts = .align_counts(counts, node_names),
    layout = .align_layout_list(layout, node_names),
    layers = .align_layers(layers, node_names),
    meta.data = .align_meta_data(meta.data, node_names),
    reductions = .align_reductions(reductions, node_names)
  )
}


# -------------------------------------------------------
# Get methods
# -------------------------------------------------------

#' Get and set CellGraph object data
#'
#' @param object A \code{\link{CellGraph}} object
#' @param slot Information to pull from object (\code{cellgraph}, \code{counts},
#' \code{layout}, \code{layers}, \code{meta.data}, \code{reductions}).
#' \code{meta_data} is accepted as an alias for \code{meta.data}.
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
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#' se <- LoadCellGraphs(se, cells = colnames(se)[1], verbose = FALSE)
#' cg <- CellGraphs(se)[[1]]
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
  object <- .upgrade_cellgraph(object)
  assert_single_value(slot, type = "string")
  slot <- .normalize_cellgraph_slot_name(slot)
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
  object <- .upgrade_cellgraph(object)
  slot <- .normalize_cellgraph_slot_name(slot)
  assert_is_one_of(slot, slotNames(x = object))

  node_names <- .cg_node_names(slot(object, name = "cellgraph"))

  if (slot == "cellgraph") {
    assert_class(value, "tbl_graph")
    .validate_cellgraph(value, verbose = FALSE)
    object <- .ensure_node_ids_on_slots(object)
    object <- .remap_cellgraph_nodes(object, .cg_node_names(value))
    slot(object, name = "cellgraph") <- value
    return(object)
  }

  if (slot == "counts") {
    slot(object, name = "counts") <- .align_counts(value, node_names)
    return(object)
  }

  if (slot == "layout") {
    slot(object, name = "layout") <- .align_layout_list(value, node_names)
    return(object)
  }

  if (slot == "layers") {
    slot(object, name = "layers") <- .align_layers(value, node_names)
    return(object)
  }

  if (slot == "meta.data") {
    slot(object, name = "meta.data") <- .align_meta_data(value, node_names)
    return(object)
  }

  if (slot == "reductions") {
    slot(object, name = "reductions") <- .align_reductions(value, node_names)
    return(object)
  }

  return(object)
}


# -------------------------------------------------------
# Seurat-style accessors
# -------------------------------------------------------

#' @param search Optional layer name or pattern passed to \code{\link[SeuratObject]{Layers}}
#'
#' @rdname CellGraph-methods
#' @method Layers CellGraph
#' @export
#'
Layers.CellGraph <- function(object, search = NULL, ...) {
  object <- .upgrade_cellgraph(object)
  lyrs <- names(slot(object, "layers"))
  if (!is.null(slot(object, "counts"))) {
    lyrs <- c("counts", lyrs)
  }
  if (!is.null(search)) {
    exact <- intersect(lyrs, search)
    if (length(exact) > 0) {
      return(exact)
    }
    lyrs <- lyrs[grepl(paste(search, collapse = "|"), lyrs)]
    if (length(lyrs) == 0) {
      return(NULL)
    }
  }
  lyrs
}

#' @param layer Name of a node matrix layer. Use \code{"counts"} for the
#' counts slot.
#'
#' @rdname CellGraph-methods
#' @method LayerData CellGraph
#' @export
#'
LayerData.CellGraph <- function(object, layer = "counts", ...) {
  object <- .upgrade_cellgraph(object)
  assert_single_value(layer, type = "string")
  if (identical(layer, "counts")) {
    return(slot(object, "counts"))
  }
  layers <- slot(object, "layers")
  if (!layer %in% names(layers)) {
    cli::cli_abort(
      c(
        "x" = "Unknown layer {.val {layer}}.",
        "i" = "Available layers: {.val {Layers(object)}}"
      )
    )
  }
  layers[[layer]]
}

#' @rdname CellGraph-methods
#' @method LayerData<- CellGraph
#' @export
#'
"LayerData<-.CellGraph" <- function(object, layer = "counts", ..., value) {
  object <- .upgrade_cellgraph(object)
  assert_single_value(layer, type = "string")
  node_names <- .cg_node_names(slot(object, "cellgraph"))
  if (identical(layer, "counts")) {
    slot(object, "counts") <- .align_counts(value, node_names)
    return(object)
  }
  if (is.null(value)) {
    slot(object, "layers")[[layer]] <- NULL
    return(object)
  }
  slot(object, "layers")[[layer]] <- .align_node_matrix(value, node_names, arg = layer)
  object
}

#' @param reduction Name of a stored \code{\link{NodeDimReduc}}. Defaults to
#' the first reduction when \code{NULL}.
#'
#' @rdname CellGraph-methods
#' @method Embeddings CellGraph
#' @export
#'
Embeddings.CellGraph <- function(object, reduction = NULL, ...) {
  Embeddings(.get_cellgraph_reduction(object, reduction))
}

#' @rdname CellGraph-methods
#' @method Loadings CellGraph
#' @export
#'
Loadings.CellGraph <- function(object, reduction = NULL, ...) {
  Loadings(.get_cellgraph_reduction(object, reduction), ...)
}

#' @rdname CellGraph-methods
#' @method Stdev CellGraph
#' @export
#'
Stdev.CellGraph <- function(object, reduction = NULL, ...) {
  Stdev(.get_cellgraph_reduction(object, reduction))
}

#' @rdname CellGraph-methods
#' @method Cells CellGraph
#' @export
#'
Cells.CellGraph <- function(x, ...) {
  .cg_node_names(slot(.upgrade_cellgraph(x), "cellgraph"))
}

#' @param metadata A vector, matrix, or \code{data.frame} of node metadata
#' @param col.name Name of the metadata column when \code{metadata} is a vector
#'
#' @rdname CellGraph-methods
#' @method AddMetaData CellGraph
#' @export
#'
AddMetaData.CellGraph <- function(object, metadata, col.name = NULL, ...) {
  object <- .upgrade_cellgraph(object)
  node_names <- .cg_node_names(slot(object, "cellgraph"))

  if (is.null(metadata)) {
    return(object)
  }

  if (is.atomic(metadata) && is.null(dim(metadata))) {
    if (is.null(col.name)) {
      cli::cli_abort("{.arg col.name} must be provided when {.arg metadata} is a vector.")
    }
    meta_names <- names(metadata)
    if (is.null(meta_names)) {
      if (length(metadata) != length(node_names)) {
        cli::cli_abort(
          c(
            "x" = "Length of {.arg metadata} ({length(metadata)}) must match",
            " " = "the number of nodes ({length(node_names)}) when it is unnamed."
          )
        )
      }
      meta_names <- node_names
    }
    metadata <- data.frame(
      x = unname(metadata),
      row.names = meta_names,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    colnames(metadata) <- col.name
  }

  new_meta <- .align_meta_data(metadata, node_names)
  old_meta <- slot(object, "meta.data")
  if (ncol(old_meta) == 0) {
    slot(object, "meta.data") <- new_meta
    return(object)
  }
  overlap <- intersect(colnames(old_meta), colnames(new_meta))
  if (length(overlap) > 0) {
    old_meta <- old_meta[, setdiff(colnames(old_meta), overlap), drop = FALSE]
  }
  slot(object, "meta.data") <- cbind(old_meta, new_meta)
  slot(object, "meta.data") <- .align_meta_data(slot(object, "meta.data"), node_names)
  object
}


# -------------------------------------------------------
# Base methods
# -------------------------------------------------------

#' CellGraph Methods
#'
#' Methods for \code{\link{CellGraph}} objects for generics defined in other
#' packages
#'
#' @param object A \code{\link{CellGraph}} object
#' @param x A \code{\link{CellGraph}} object
#' @param i Name of a stored reduction
#' @param j Unused
#' @param nodes A character vector of node names
#' @param value Replacement value
#' @param drop Unused
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
    object <- .upgrade_cellgraph(object)
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
    extra_layers <- names(slot(object, "layers"))
    if (length(extra_layers) > 0) {
      cat("Layers:", col_br_blue(paste(extra_layers, collapse = ", ")), "\n")
    }
    meta_cols <- colnames(slot(object, "meta.data"))
    if (length(meta_cols) > 0) {
      cat("Node metadata:", col_br_blue(paste(meta_cols, collapse = ", ")), "\n")
    }
    dr_names <- names(slot(object, "reductions"))
    if (length(dr_names) > 0) {
      cat("Reductions:", col_br_blue(paste(dr_names, collapse = ", ")), "\n")
    }
  }
)

#' Extract a dimensionality reduction from a CellGraph
#'
#' @describeIn CellGraph-methods Extract a \code{NodeDimReduc} by name
#'
#' @export
#'
setMethod(
  f = "[[",
  signature = c("x" = "CellGraph", "i" = "character", "j" = "missing"),
  definition = function(x, i, j, ..., drop = TRUE) {
    x <- .upgrade_cellgraph(x)
    reductions <- slot(x, "reductions")
    if (!i %in% names(reductions)) {
      cli::cli_abort(
        c(
          "x" = "Unknown reduction {.val {i}}.",
          "i" = "Available reductions: {.val {names(reductions)}}"
        )
      )
    }
    reductions[[i]]
  }
)

#' Add or replace a dimensionality reduction in a CellGraph
#'
#' @describeIn CellGraph-methods Add or replace a \code{NodeDimReduc}
#'
#' @export
#'
setMethod(
  f = "[[<-",
  signature = c("x" = "CellGraph", "i" = "character", "j" = "missing", "value" = "ANY"),
  definition = function(x, i, j, ..., value) {
    x <- .upgrade_cellgraph(x)
    node_names <- .cg_node_names(slot(x, "cellgraph"))
    if (is.null(value)) {
      slot(x, "reductions")[[i]] <- NULL
      return(x)
    }
    slot(x, "reductions")[[i]] <- .align_node_dimreduc(value, node_names, arg = i)
    x
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
  x <- .upgrade_cellgraph(x)
  assert_vector(nodes, type = "character")
  available_nodes <- .cg_node_names(slot(x, "cellgraph"))
  assert_x_in_y(nodes, available_nodes)

  x <- .ensure_node_ids_on_slots(x)
  graph_type <- attr(x@cellgraph, "type")
  x@cellgraph <- x@cellgraph %N>% filter(name %in% nodes)
  attr(x@cellgraph, "type") <- graph_type
  .remap_cellgraph_nodes(x, .cg_node_names(x@cellgraph))
}


# -------------------------------------------------------
# Internal helpers
# -------------------------------------------------------

#' Normalize CellGraphData slot names
#'
#' @noRd
#'
.normalize_cellgraph_slot_name <- function(slot) {
  if (identical(slot, "meta_data")) {
    return("meta.data")
  }
  slot
}

#' Validate a tbl_graph for use in a CellGraph
#'
#' @noRd
#'
.validate_cellgraph <- function(cellgraph, verbose = FALSE, call = caller_env()) {
  if (!"type" %in% names(attributes(cellgraph))) {
    cli::cli_abort(c("x" = "Graph attribute {.str type} is missing."), call = call)
  } else if (verbose && check_global_verbosity()) {
    cli::cli_alert_info("Got a graph of type '{attr(cellgraph, 'type')}'")
  }

  if (attr(cellgraph, "type") == "bipartite") {
    if (!"name" %in% vertex_attr_names(cellgraph)) {
      cli::cli_abort("x" = "Node attribute {.str name} is missing from the graph", call = call)
    }
    if (!"node_type" %in% vertex_attr_names(cellgraph)) {
      cli::cli_abort("x" = "Node attribute {.str node_type} is missing from the graph", call = call)
    }
  }
  # TODO: Add check for A-node-projection and linegraph

  node_names <- .cg_node_names(cellgraph)
  if (anyDuplicated(node_names)) {
    cli::cli_abort(c("x" = "Node names in {.arg cellgraph} must be unique."), call = call)
  }
}

#' Node names for a tbl_graph
#'
#' @noRd
#'
.cg_node_names <- function(cellgraph) {
  if ("name" %in% vertex_attr_names(cellgraph)) {
    return(as.character(cellgraph %N>% pull(name)))
  }
  as.character(seq_len(length(cellgraph)))
}

#' Upgrade CellGraph objects created with fewer slots
#'
#' @noRd
#'
.upgrade_cellgraph <- function(object) {
  if (!is(object, "CellGraph")) {
    return(object)
  }
  if (!methods::.hasSlot(object, "layers")) {
    object <- new(
      Class = "CellGraph",
      cellgraph = slot(object, "cellgraph"),
      counts = slot(object, "counts"),
      layout = slot(object, "layout")
    )
  }
  if (is.null(slot(object, "layers"))) {
    slot(object, "layers") <- list()
  }
  if (is.null(slot(object, "reductions"))) {
    slot(object, "reductions") <- list()
  }
  if (is.null(slot(object, "meta.data"))) {
    slot(object, "meta.data") <- data.frame(row.names = .cg_node_names(slot(object, "cellgraph")))
  }
  object
}

#' Match row names of a matrix to node names
#'
#' @noRd
#'
.align_matrix_rows <- function(mat, node_names, arg = "matrix", call = caller_env()) {
  if (is.null(mat)) {
    return(NULL)
  }
  if (is.null(rownames(mat))) {
    if (nrow(mat) != length(node_names)) {
      cli::cli_abort(
        c(
          "x" = "{.arg {arg}} has no row names and {nrow(mat)} row{?s},",
          " " = "but the graph has {length(node_names)} node{?s}."
        ),
        call = call
      )
    }
    rownames(mat) <- node_names
    return(mat)
  }
  rownames(mat) <- as.character(rownames(mat))
  if (anyDuplicated(rownames(mat))) {
    cli::cli_abort(c("x" = "Row names in {.arg {arg}} must be unique."), call = call)
  }
  missing_nodes <- setdiff(node_names, rownames(mat))
  if (length(missing_nodes) > 0) {
    cli::cli_abort(
      c(
        "x" = "{.arg {arg}} is missing {length(missing_nodes)} node{?s} present in the graph.",
        "i" = "Example: {.val {head(missing_nodes, 3)}}"
      ),
      call = call
    )
  }
  mat[node_names, , drop = FALSE]
}

#' Align the counts matrix
#'
#' @noRd
#'
.align_counts <- function(counts, node_names, call = caller_env()) {
  if (is.null(counts)) {
    return(NULL)
  }
  assert_class(counts, "dgCMatrix", call = call)
  .align_matrix_rows(counts, node_names, arg = "counts", call = call)
}

#' Align a numeric node matrix used as a layer
#'
#' @noRd
#'
.align_node_matrix <- function(mat, node_names, arg = "layer", call = caller_env()) {
  if (inherits(mat, "data.frame")) {
    mat <- as.matrix(mat)
  }
  if (!(is.matrix(mat) || inherits(mat, "Matrix"))) {
    cli::cli_abort(
      c(
        "x" = "{.arg {arg}} must be a numeric matrix.",
        "i" = "Got a {.cls {class(mat)}}."
      ),
      call = call
    )
  }
  if (nrow(mat) > 0 && ncol(mat) > 0 && !is.numeric(mat[1, 1])) {
    cli::cli_abort(c("x" = "{.arg {arg}} must be numeric."), call = call)
  }
  .align_matrix_rows(mat, node_names, arg = arg, call = call)
}

#' Align a named list of layers
#'
#' @noRd
#'
.align_layers <- function(layers, node_names, call = caller_env()) {
  if (is.null(layers) || length(layers) == 0) {
    return(list())
  }
  if (is.null(names(layers)) || any(names(layers) == "")) {
    cli::cli_abort("The {.arg layers} list must be named.", call = call)
  }
  if ("counts" %in% names(layers)) {
    cli::cli_abort(
      c(
        "x" = "{.str counts} is a reserved layer name.",
        "i" = "Store the count matrix in the {.arg counts} argument / slot."
      ),
      call = call
    )
  }
  aligned <- lapply(names(layers), function(nm) {
    .align_node_matrix(layers[[nm]], node_names, arg = nm, call = call)
  })
  names(aligned) <- names(layers)
  aligned
}

#' Node identifiers stored as row names
#'
#' Returns \code{NULL} for automatic row names, which \code{data.frame} uses
#' when no identifiers have been set (tibbles never carry row names, so they
#' also end up with automatic row names once coerced).
#'
#' @noRd
#'
.explicit_rownames <- function(x) {
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (.row_names_info(x) < 0L) {
    return(NULL)
  }
  as.character(attr(x, "row.names"))
}

#' Align a layout table to graph node order
#'
#' @noRd
#'
.align_layout <- function(layout, node_names, layout_name = "layout", call = caller_env()) {
  if (!inherits(layout, "data.frame")) {
    cli::cli_abort(
      c("x" = "The '{layout_name}' layout table must be a {.cls data.frame}"),
      call = call
    )
  }

  layout <- as.data.frame(layout, stringsAsFactors = FALSE, check.names = FALSE)
  if ("name" %in% colnames(layout)) {
    layout_names <- as.character(layout$name)
    layout <- layout[, setdiff(colnames(layout), "name"), drop = FALSE]
  } else {
    layout_names <- .explicit_rownames(layout)
  }

  if (is.null(layout_names)) {
    if (nrow(layout) != length(node_names)) {
      cli::cli_abort(
        c(
          "x" = "Number of nodes ({length(node_names)}) in the 'cellgraph' slot does not match ",
          " " = "the number of rows ({nrow(layout)}) in the '{layout_name}' layout table"
        ),
        call = call
      )
    }
    layout_names <- node_names
  }

  if (anyDuplicated(layout_names)) {
    cli::cli_abort(
      c("x" = "Node names in the '{layout_name}' layout must be unique."),
      call = call
    )
  }
  missing_nodes <- setdiff(node_names, layout_names)
  if (length(missing_nodes) > 0) {
    cli::cli_abort(
      c(
        "x" = "The '{layout_name}' layout is missing {length(missing_nodes)} node{?s}.",
        "i" = "Example: {.val {head(missing_nodes, 3)}}"
      ),
      call = call
    )
  }
  layout <- layout[match(node_names, layout_names), , drop = FALSE]
  rownames(layout) <- node_names
  layout
}

#' Align a named list of layouts
#'
#' @noRd
#'
.align_layout_list <- function(layout, node_names, call = caller_env()) {
  if (is.null(layout)) {
    return(NULL)
  }
  assert_non_empty_object(layout, "list", call = call)
  if (is.null(names(layout)) || any(names(layout) == "")) {
    cli::cli_abort("The {.arg layout} list must be named.", call = call)
  }
  aligned <- lapply(names(layout), function(nm) {
    .align_layout(layout[[nm]], node_names, layout_name = nm, call = call)
  })
  names(aligned) <- names(layout)
  aligned
}

#' Align node metadata
#'
#' @noRd
#'
.align_meta_data <- function(meta, node_names, call = caller_env()) {
  if (is.null(meta) || (nrow(as.data.frame(meta)) == 0 && ncol(as.data.frame(meta)) == 0)) {
    return(data.frame(row.names = node_names))
  }
  meta <- as.data.frame(meta, stringsAsFactors = FALSE, check.names = FALSE)
  meta_names <- .explicit_rownames(meta)
  if ("name" %in% colnames(meta) && is.null(meta_names)) {
    meta_names <- as.character(meta$name)
    meta$name <- NULL
  }
  if (is.null(meta_names)) {
    if (nrow(meta) != length(node_names)) {
      cli::cli_abort(
        c(
          "x" = "{.arg meta.data} has no node identifiers and {nrow(meta)} row{?s},",
          " " = "but the graph has {length(node_names)} node{?s}."
        ),
        call = call
      )
    }
    meta_names <- node_names
  }
  if (anyDuplicated(meta_names)) {
    cli::cli_abort(c("x" = "Node names in {.arg meta.data} must be unique."), call = call)
  }
  rownames(meta) <- meta_names
  missing_nodes <- setdiff(node_names, meta_names)
  if (length(missing_nodes) > 0) {
    cli::cli_abort(
      c(
        "x" = "{.arg meta.data} is missing {length(missing_nodes)} node{?s} present in the graph.",
        "i" = "Example: {.val {head(missing_nodes, 3)}}"
      ),
      call = call
    )
  }
  meta[node_names, , drop = FALSE]
}

#' Align a named list of NodeDimReduc objects
#'
#' @noRd
#'
.align_reductions <- function(reductions, node_names, call = caller_env()) {
  if (is.null(reductions) || length(reductions) == 0) {
    return(list())
  }
  if (is.null(names(reductions)) || any(names(reductions) == "")) {
    cli::cli_abort("The {.arg reductions} list must be named.", call = call)
  }
  aligned <- lapply(names(reductions), function(nm) {
    .align_node_dimreduc(reductions[[nm]], node_names, arg = nm, call = call)
  })
  names(aligned) <- names(reductions)
  aligned
}

#' Fetch a named reduction from a CellGraph
#'
#' @noRd
#'
.get_cellgraph_reduction <- function(object, reduction = NULL, call = caller_env()) {
  object <- .upgrade_cellgraph(object)
  reductions <- slot(object, "reductions")
  if (length(reductions) == 0) {
    cli::cli_abort(c("x" = "This {.cls CellGraph} has no reductions."), call = call)
  }
  if (is.null(reduction)) {
    reduction <- names(reductions)[1]
  }
  assert_single_value(reduction, type = "string", call = call)
  if (!reduction %in% names(reductions)) {
    cli::cli_abort(
      c(
        "x" = "Unknown reduction {.val {reduction}}.",
        "i" = "Available reductions: {.val {names(reductions)}}"
      ),
      call = call
    )
  }
  reductions[[reduction]]
}

#' Attach node identities to slots that still rely on positional alignment
#'
#' @noRd
#'
.ensure_node_ids_on_slots <- function(object) {
  node_names <- .cg_node_names(slot(object, "cellgraph"))

  counts <- slot(object, "counts")
  if (!is.null(counts) && is.null(rownames(counts)) && nrow(counts) == length(node_names)) {
    rownames(counts) <- node_names
    slot(object, "counts") <- counts
  }

  layouts <- slot(object, "layout")
  if (!is.null(layouts) && length(layouts) > 0) {
    slot(object, "layout") <- lapply(layouts, function(ly) {
      has_names <- !is.null(.explicit_rownames(ly)) || "name" %in% colnames(ly)
      if (!has_names && nrow(ly) == length(node_names)) {
        ly <- as.data.frame(ly, stringsAsFactors = FALSE, check.names = FALSE)
        rownames(ly) <- node_names
      }
      ly
    })
  }

  meta <- slot(object, "meta.data")
  if (ncol(meta) > 0 && is.null(.explicit_rownames(meta)) && nrow(meta) == length(node_names)) {
    rownames(meta) <- node_names
    slot(object, "meta.data") <- meta
  }

  layers <- slot(object, "layers")
  if (length(layers) > 0) {
    slot(object, "layers") <- lapply(layers, function(mat) {
      if (is.null(rownames(mat)) && nrow(mat) == length(node_names)) {
        rownames(mat) <- node_names
      }
      mat
    })
  }

  object
}

#' Reorder or subset all node-level slots to \code{node_names}
#'
#' @noRd
#'
.remap_cellgraph_nodes <- function(object, node_names, call = caller_env()) {
  counts <- slot(object, "counts")
  if (!is.null(counts)) {
    slot(object, "counts") <- .align_counts(counts, node_names, call = call)
  }

  layouts <- slot(object, "layout")
  if (!is.null(layouts) && length(layouts) > 0) {
    slot(object, "layout") <- .align_layout_list(layouts, node_names, call = call)
  }

  layers <- slot(object, "layers")
  if (length(layers) > 0) {
    slot(object, "layers") <- .align_layers(layers, node_names, call = call)
  }

  meta <- slot(object, "meta.data")
  if (ncol(meta) == 0 && (nrow(meta) == 0 || !all(node_names %in% rownames(meta)))) {
    slot(object, "meta.data") <- data.frame(row.names = node_names)
  } else {
    slot(object, "meta.data") <- .align_meta_data(meta, node_names, call = call)
  }

  reductions <- slot(object, "reductions")
  if (length(reductions) > 0) {
    slot(object, "reductions") <- .align_reductions(reductions, node_names, call = call)
  }

  object
}
