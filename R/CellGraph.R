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
#' A \code{CellGraph} contains counts and a graph, and optionally layouts, layers,
#' metadata, or reductions.
#'
#' Node-level variable names must be unique across the graph node table,
#' \code{meta.data}, and reduction embeddings, and must not overlap count or
#' layer features. Count and layer matrices are the exception: they may share
#' feature names because callers select a layer explicitly.
#'
#' Objects serialized before these extra slots existed are not upgraded.
#' Using one aborts with a message that names the missing slots and points
#' to a \code{pixelatorR} version that still reads the old class.
#'
#' @slot cellgraph A \code{tbl_graph} object corresponding to a cell graph
#' @slot nodes Character vector of node IDs in graph order. This is the map
#' used to align counts, layouts, layers, metadata, and reductions. Those
#' tables are stored in this order without copying the IDs as row names.
#' @slot counts A \code{matrix}-like object with marker counts (nodes x markers).
#' Rows follow \code{nodes}. The counts matrix can be extracted as the
#' \code{"counts"} layer via \code{\link[SeuratObject]{Layers}} /
#' \code{\link[SeuratObject]{LayerData}}.
#' @slot layout A named \code{list} of \code{data.frame} objects with coordinates
#' for cell layouts. Rows follow \code{nodes}. A \code{name} column or explicit
#' row names are accepted on input and used only to reorder; MPX bipartite
#' layouts may omit \code{-A}/\code{-B} suffixes. Stored layouts keep
#' coordinate columns (typically \code{x}, \code{y}, \code{z}). Layouts without
#' node IDs still work if the number of rows matches the graph.
#' @slot layers A named \code{list} of additional numeric node matrices
#' (nodes x features), analogous to layers on a Seurat
#' \code{\link[SeuratObject]{Assay5}}. A layer can be extracted
#' via \code{\link[SeuratObject]{Layers}} / \code{\link[SeuratObject]{LayerData}}.
#' @slot meta.data A \code{data.frame} of node-level metadata (one row per node).
#' Rows follow \code{nodes}. Columns may have mixed types.
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
    nodes = "character",
    counts = "ANY",
    layout = "ANY",
    layers = "list",
    meta.data = "data.frame",
    reductions = "list"
  )
)

#' Initialize a CellGraph object
#'
#' Supplies default values for the extra slots so
#' \code{methods::new("CellGraph")} and constructors can omit
#' \code{nodes}, \code{layers}, \code{meta.data}, and \code{reductions}.
#' When a graph is provided, \code{nodes} is filled from the graph and an
#' empty \code{meta.data} table has one row per node.
#'
#' @param .Object A \code{CellGraph} instance being constructed
#' @param cellgraph A \code{tbl_graph}, or \code{NULL}
#' @param nodes Character vector of node IDs
#' @param counts A count matrix, or \code{NULL}
#' @param layout A named list of layout tables, or \code{NULL}
#' @param layers A named list of extra node matrices
#' @param meta.data A node-level \code{data.frame}
#' @param reductions A named list of \code{NodeDimReduc} objects
#' @param ... Passed to the next \code{initialize} method
#'
#' @return A \code{CellGraph} object
#'
#' @keywords internal
#' @noRd
#'
setMethod(
  f = "initialize",
  signature = "CellGraph",
  definition = function(
    .Object,
    cellgraph = NULL,
    nodes = character(),
    counts = NULL,
    layout = NULL,
    layers = list(),
    meta.data = data.frame(),
    reductions = list(),
    ...
  ) {
    if (is.null(layers)) {
      layers <- list()
    }
    if (is.null(reductions)) {
      reductions <- list()
    }
    if (is.null(meta.data)) {
      meta.data <- data.frame()
    }
    if (is.null(nodes)) {
      nodes <- character()
    }
    if (!is.null(cellgraph) && length(nodes) == 0) {
      nodes <- .cg_node_names(cellgraph)
    }
    if (!is.null(cellgraph) && ncol(as.data.frame(meta.data)) == 0) {
      meta.data <- .empty_node_meta(length(nodes))
    }
    callNextMethod(
      .Object,
      cellgraph = cellgraph,
      nodes = nodes,
      counts = counts,
      layout = layout,
      layers = layers,
      meta.data = meta.data,
      reductions = reductions,
      ...
    )
  }
)


# -------------------------------------------------------
# Create methods
# -------------------------------------------------------

#' Create a CellGraph object
#'
#' @details
#' Node-level variable names must not clash between the graph node table,
#' \code{meta.data}, reduction embeddings, and matrix features. Count and
#' layer matrices may share feature names because methods such as
#' \code{\link[SeuratObject]{FetchData}} select a specific layer.
#'
#' @param cellgraph A \code{tbl_graph} object representing a PNA single-cell graph
#' @param counts A \code{dgCMatrix} with marker counts. Rows are matched to graph
#' node names (order does not need to match).
#' @param layout A named \code{list} of \code{data.frame} objects with cell
#' layouts. Nodes are identified by row names or by a \code{name} column;
#' otherwise the row order is assumed to follow the graph. MPX bipartite
#' layouts may use unsuffixed names while graph nodes keep \code{-A}/\code{-B};
#' those names are matched after stripping the suffix, as in
#' \code{\link{LoadCellGraphs}}. Stored layouts keep graph node order and do
#' not copy node IDs as row names.
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
#' # Create a CellGraph object with graph and counts
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
  cellgraph <- .cg_graph_with_node_names(cellgraph)

  node_names <- .cg_node_names(cellgraph)

  object <- new(
    Class = "CellGraph",
    cellgraph = cellgraph,
    nodes = node_names,
    counts = .align_counts(counts, node_names),
    layout = .align_layout_list(layout, node_names),
    layers = .align_layers(layers, node_names),
    meta.data = .align_meta_data(meta.data, node_names),
    reductions = .align_reductions(reductions, node_names)
  )
  .validate_cellgraph_data_names(object)
  object
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
  assert_single_value(slot, type = "string")
  .assert_current_cellgraph(object)
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
  .assert_current_cellgraph(object)
  slot <- .normalize_cellgraph_slot_name(slot)
  assert_is_one_of(slot, slotNames(x = object))
  node_names <- .cg_node_map(object)

  if (slot == "nodes") {
    cli::cli_abort(
      c(
        "x" = "The {.field nodes} map cannot be replaced directly.",
        "i" = "Replace {.field cellgraph} or use {.fn subset}."
      )
    )
  }

  if (slot == "cellgraph") {
    assert_class(value, "tbl_graph")
    .validate_cellgraph(value, verbose = FALSE)
    value <- .cg_graph_with_node_names(value)
    object <- .ensure_node_ids_on_slots(object)
    object <- .remap_cellgraph_nodes(object, .cg_node_names(value))
    slot(object, name = "cellgraph") <- value
    .validate_cellgraph_data_names(object)
    return(object)
  }

  if (slot == "counts") {
    slot(object, name = "counts") <- .align_counts(value, node_names)
    .validate_cellgraph_data_names(object)
    return(object)
  }

  if (slot == "layout") {
    slot(object, name = "layout") <- .align_layout_list(value, node_names)
    return(object)
  }

  if (slot == "layers") {
    slot(object, name = "layers") <- .align_layers(value, node_names)
    .validate_cellgraph_data_names(object)
    return(object)
  }

  if (slot == "meta.data") {
    slot(object, name = "meta.data") <- .align_meta_data(value, node_names)
    .validate_cellgraph_data_names(object)
    return(object)
  }

  if (slot == "reductions") {
    slot(object, name = "reductions") <- .align_reductions(value, node_names)
    .validate_cellgraph_data_names(object)
    return(object)
  }

  return(object)
}


# -------------------------------------------------------
# Seurat-style accessors
# -------------------------------------------------------

#' @param search Optional layer name or pattern passed to \code{\link[SeuratObject]{Layers}}
#' @param layer Name of a node matrix layer. Use \code{"counts"} for the
#' counts slot. For \code{FetchData}, \code{NULL} (default) selects
#' \code{"counts"} when present, otherwise the first extra layer.
#'
#' @rdname CellGraph-methods
#' @method Layers CellGraph
#' @export
#'
Layers.CellGraph <- function(object, search = NULL, ...) {
  .assert_current_cellgraph(object)
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

#' @rdname CellGraph-methods
#' @method LayerData CellGraph
#' @export
#'
LayerData.CellGraph <- function(object, layer = "counts", ...) {
  assert_single_value(layer, type = "string")
  .assert_current_cellgraph(object)
  if (identical(layer, "counts")) {
    return(.with_row_ids(slot(object, "counts"), .cg_node_map(object)))
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
  .with_row_ids(layers[[layer]], .cg_node_map(object))
}

#' @rdname CellGraph-methods
#' @method LayerData<- CellGraph
#' @export
#'
"LayerData<-.CellGraph" <- function(object, layer = "counts", ..., value) {
  assert_single_value(layer, type = "string")
  .assert_current_cellgraph(object)
  node_names <- .cg_node_map(object)
  if (identical(layer, "counts")) {
    slot(object, "counts") <- .align_counts(value, node_names)
    .validate_cellgraph_data_names(object)
    return(object)
  }
  if (is.null(value)) {
    slot(object, "layers")[[layer]] <- NULL
    return(object)
  }
  slot(object, "layers")[[layer]] <- .align_node_matrix(value, node_names, arg = layer)
  .validate_cellgraph_data_names(object)
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
  .with_row_ids(
    Embeddings(.get_cellgraph_reduction(object, reduction)),
    .cg_node_map(object)
  )
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
  .assert_current_cellgraph(x)
  .cg_node_map(x)
}

#' @param metadata A vector, matrix, or \code{data.frame} of node metadata.
#' Nodes are matched by name, so metadata may cover a subset of the graph;
#' the remaining nodes get \code{NA}. Names that are not graph nodes are
#' dropped.
#' @param col.name Name of the metadata column when \code{metadata} is a vector
#'
#' @rdname CellGraph-methods
#' @method AddMetaData CellGraph
#' @export
#'
AddMetaData.CellGraph <- function(object, metadata, col.name = NULL, ...) {
  .assert_current_cellgraph(object)
  node_names <- .cg_node_map(object)

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
    if (anyDuplicated(meta_names)) {
      cli::cli_abort(c("x" = "Names in {.arg metadata} must be unique."))
    }
    metadata <- data.frame(
      x = unname(metadata),
      row.names = meta_names,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    colnames(metadata) <- col.name
  }

  new_meta <- .fill_meta_data(metadata, node_names)
  old_meta <- slot(object, "meta.data")
  if (ncol(old_meta) == 0) {
    slot(object, "meta.data") <- new_meta
    .validate_cellgraph_data_names(object)
    return(object)
  }
  overlap <- intersect(colnames(old_meta), colnames(new_meta))
  if (length(overlap) > 0) {
    old_meta <- old_meta[, setdiff(colnames(old_meta), overlap), drop = FALSE]
  }
  slot(object, "meta.data") <- .cbind_keep_names(old_meta, new_meta)
  slot(object, "meta.data") <- .align_meta_data(slot(object, "meta.data"), node_names)
  .validate_cellgraph_data_names(object)
  object
}

#' @param vars Variables to fetch: marker names, node metadata columns,
#' graph vertex attributes, or reduction embedding columns (for example
#' \code{"PC_1"}).
#' @param cells For \code{FetchData.CellGraph}, nodes to collect (default is
#' all nodes). Numeric indices are allowed, matching
#' \code{\link[SeuratObject]{FetchData}}. For \code{FetchData.CellGraphList},
#' component IDs (default is all loaded graphs). Unloaded graphs raise an
#' error when they are included in \code{cells}.
#' @param clean If \code{TRUE}, remove nodes that are missing data for every
#' requested variable. \code{FetchData.CellGraph} defaults to \code{TRUE}.
#' \code{FetchData.CellGraphList} defaults to \code{FALSE} so graphs that
#' lack the requested variables still appear with \code{NA} values.
#'
#' @details
#' Variable names must be unique across the graph node table,
#' \code{meta.data}, reduction embeddings, and matrix features. Count and
#' layer matrices may share feature names because \code{layer} selects the
#' matrix to search.
#'
#' @describeIn CellGraph-methods Pull node-level data from a \code{CellGraph}
#' @method FetchData CellGraph
#' @export
#'
FetchData.CellGraph <- function(
  object,
  vars,
  cells = NULL,
  layer = NULL,
  clean = TRUE,
  ...
) {
  .validate_cellgraph_data_names(object)
  node_names <- Cells(object)

  if (isTRUE(clean)) {
    clean <- "all"
  } else if (isFALSE(clean)) {
    clean <- "none"
  }
  clean <- rlang::arg_match0(clean, values = c("all", "none"))

  cells <- cells %||% node_names
  if (is.numeric(cells)) {
    cells <- node_names[cells]
  }
  assert_vector(cells, type = "character", n = 1)
  cells <- unique(as.character(cells))
  missing_cells <- setdiff(cells, node_names)
  cells <- intersect(cells, node_names)
  if (length(cells) == 0) {
    cli::cli_abort(c("x" = "None of the requested nodes were found in this {.cls CellGraph}."))
  }
  if (length(missing_cells) > 0) {
    cli::cli_warn("Removing {length(missing_cells)} node{?s} not present in this {.cls CellGraph}.")
  }

  if (is.null(vars) || length(vars) == 0) {
    return(data.frame(row.names = cells))
  }
  assert_vector(vars, type = "character", n = 1)
  vars <- as.character(vars)

  data_fetched <- data.frame(row.names = cells)

  # Pull vars from node metadata
  meta <- slot(object, "meta.data")
  meta_vars <- intersect(vars, colnames(meta))
  if (length(meta_vars) > 0) {
    data_fetched <- .add_fetched_cols(
      data_fetched,
      .with_row_ids(meta[.row_index(cells, node_names), meta_vars, drop = FALSE], cells)
    )
  }

  # Pull remaining vars from graph vertex attributes
  graph_meta <- .cg_vertex_attr_df(slot(object, "cellgraph"))
  graph_vars <- setdiff(intersect(vars, colnames(graph_meta)), names(data_fetched))
  if (length(graph_vars) > 0) {
    data_fetched <- .add_fetched_cols(data_fetched, graph_meta[cells, graph_vars, drop = FALSE])
  }

  # Pull keyed embedding columns from reductions
  remaining <- setdiff(vars, names(data_fetched))
  if (length(remaining) > 0) {
    reductions <- slot(object, "reductions")
    for (nm in names(reductions)) {
      remaining <- setdiff(vars, names(data_fetched))
      if (length(remaining) == 0) {
        break
      }
      data_fetched <- .add_fetched_cols(
        data_fetched,
        .fetch_nodedimreduc_vars(reductions[[nm]], remaining, cells, node_names)
      )
    }
  }

  # Pull remaining vars from a node layer (markers / extra layers)
  remaining <- setdiff(vars, names(data_fetched))
  available_layers <- Layers(object)
  if (length(remaining) > 0 && length(available_layers) > 0) {
    if (is.null(layer)) {
      layer <- if ("counts" %in% available_layers) "counts" else available_layers[[1]]
    }
    assert_single_value(layer, type = "string")
    if (!layer %in% available_layers) {
      cli::cli_abort(
        c(
          "x" = "Unknown layer {.val {layer}}.",
          "i" = "Available layers: {.val {available_layers}}"
        )
      )
    }
    data_fetched <- .add_fetched_cols(
      data_fetched,
      .fetch_layer_vars(object, layer, remaining, cells)
    )
    remaining <- setdiff(vars, names(data_fetched))
    other_layers <- setdiff(available_layers, layer)
    if (length(remaining) > 0 && length(other_layers) > 0) {
      data_fetched <- .add_fetched_cols(
        data_fetched,
        .fetch_vars_from_other_layers(object, remaining, cells, other_layers)
      )
    }
  } else if (!is.null(layer) && length(available_layers) == 0) {
    cli::cli_abort(
      c(
        "x" = "Unknown layer {.val {layer}}.",
        "i" = "This {.cls CellGraph} has no layers."
      )
    )
  }

  vars_missing <- setdiff(vars, names(data_fetched))
  m2 <- if (length(vars_missing) > 10) {
    paste0(" (10 out of ", length(vars_missing), " shown)")
  } else {
    ""
  }
  if (length(vars_missing) == length(vars)) {
    cli::cli_abort(
      c("x" = "None of the requested variables were found{m2}: {.val {head(vars_missing, 10)}}")
    )
  } else if (length(vars_missing) > 0) {
    cli::cli_warn("The following requested variables were not found{m2}: {.val {head(vars_missing, 10)}}")
  }

  found <- intersect(vars, names(data_fetched))
  data_fetched <- data_fetched[, found, drop = FALSE]

  if (identical(clean, "all")) {
    no_data <- which(apply(data_fetched, 1L, function(x) all(is.na(x))))
    if (length(no_data) > 0) {
      cli::cli_warn("Removing {length(no_data)} node{?s} missing data for vars requested")
      data_fetched <- data_fetched[-no_data, , drop = FALSE]
    }
  }
  data_fetched
}


# -------------------------------------------------------
# Base methods
# -------------------------------------------------------

#' CellGraph Methods
#'
#' Methods for \code{\link{CellGraph}} objects for generics defined in other
#' packages
#'
#' @param object A \code{\link{CellGraph}} or \code{\link{CellGraphList}} object
#' @param x A \code{\link{CellGraph}} object
#' @param i Name of a stored reduction
#' @param j,drop Required by the S4 \code{[[} generic and ignored
#' @param nodes A character vector of node names
#' @param value Replacement value
#' @param ... Currently not used
#'
#' @return \code{FetchData.CellGraph}: a \code{data.frame} with nodes as rows
#' and requested variables as columns. \code{FetchData.CellGraphList}: a
#' \code{data.frame} with a \code{component} column identifying the source graph
#' and the requested variables. Row names are \code{component:node}.
#' \code{subset}: a \code{CellGraph} object containing only the specified nodes.
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
#' # Fetch marker counts, node attributes, or embeddings
#' head(SeuratObject::FetchData(cg, vars = colnames(cg@counts)[1]))
#'
#' # FetchData.CellGraphList combines loaded graphs without requiring a layout
#' cgl <- CellGraphs(se)
#' head(SeuratObject::FetchData(cgl, vars = colnames(cg@counts)[1]))
#'
setMethod(
  f = "show",
  signature = "CellGraph",
  definition = function(object) {
    .assert_current_cellgraph(object)
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
    .assert_current_cellgraph(x)
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
    .assert_current_cellgraph(x)
    node_names <- .cg_node_map(x)
    if (is.null(value)) {
      slot(x, "reductions")[[i]] <- NULL
      return(x)
    }
    slot(x, "reductions")[[i]] <- .align_node_dimreduc(value, node_names, arg = i)
    .validate_cellgraph_data_names(x)
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
#' cg_small <- subset(cg, nodes = Cells(cg)[1:100])
#' cg_small
#'
#' @export
#'
subset.CellGraph <- function(
  x,
  nodes,
  ...
) {
  assert_vector(nodes, type = "character", n = 1)
  .assert_current_cellgraph(x)
  available_nodes <- .cg_node_map(x)
  assert_x_in_y(nodes, available_nodes)

  x <- .ensure_node_ids_on_slots(x)
  graph_type <- attr(x@cellgraph, "type")
  x@cellgraph <- .cg_graph_with_node_names(x@cellgraph) %N>% filter(name %in% nodes)
  attr(x@cellgraph, "type") <- graph_type
  .remap_cellgraph_nodes(x, .cg_node_names(x@cellgraph))
}


# -------------------------------------------------------
# Internal helpers
# -------------------------------------------------------

#' Normalize CellGraphData slot names
#'
#' Maps the \code{meta_data} alias to the \code{meta.data} slot name so
#' getters and setters accept either spelling.
#'
#' @param slot Character slot name from \code{CellGraphData}
#'
#' @return The canonical slot name
#'
#' @keywords internal
#' @noRd
#'
.normalize_cellgraph_slot_name <- function(slot) {
  if (identical(slot, "meta_data")) {
    return("meta.data")
  }
  slot
}

#' Canonical node IDs for a CellGraph
#'
#' Prefers the \code{nodes} slot. Falls back to graph vertex names when
#' the slot is empty.
#'
#' @param object A \code{CellGraph}
#'
#' @return Character vector of node IDs
#'
#' @keywords internal
#' @noRd
#'
.cg_node_map <- function(object) {
  nodes <- slot(object, "nodes")
  if (length(nodes) > 0) {
    return(nodes)
  }
  cellgraph <- slot(object, "cellgraph")
  if (is.null(cellgraph)) {
    return(character())
  }
  .cg_node_names(cellgraph)
}

#' Empty node metadata with one row per node
#'
#' Uses compact automatic row names so node IDs are not stored twice.
#'
#' @param n Number of nodes
#'
#' @return A zero-column \code{data.frame} with \code{n} rows
#'
#' @keywords internal
#' @noRd
#'
.empty_node_meta <- function(n) {
  as.data.frame(matrix(nrow = n, ncol = 0))
}

#' Drop stored row names from a node-level table
#'
#' @param x A matrix or data frame, or \code{NULL}
#'
#' @return \code{x} with \code{rownames} unset
#'
#' @keywords internal
#' @noRd
#'
.drop_row_ids <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  rownames(x) <- NULL
  x
}

#' Attach node IDs as row names for user-facing extracts
#'
#' @param x A matrix or data frame, or \code{NULL}
#' @param node_names Character vector of node IDs
#'
#' @return \code{x} with \code{rownames} set to \code{node_names}
#'
#' @keywords internal
#' @noRd
#'
.with_row_ids <- function(x, node_names) {
  if (is.null(x)) {
    return(NULL)
  }
  rownames(x) <- node_names
  x
}

#' Match requested node IDs to a name order
#'
#' @param cells Requested node IDs
#' @param node_names Node IDs in stored row order
#'
#' @return Integer indices into \code{node_names}
#'
#' @keywords internal
#' @noRd
#'
.row_index <- function(cells, node_names) {
  match(cells, node_names)
}

#' Validate a tbl_graph for use in a CellGraph
#'
#' Checks that the graph has a \code{type} attribute, unique node names,
#' and (for bipartite graphs) \code{name} and \code{node_type} vertex
#' attributes. Errors are reported from \code{call} so they point at the
#' user-facing constructor or setter, not this helper.
#'
#' @param cellgraph A \code{tbl_graph}
#' @param verbose Print the detected graph type
#' @param call Environment to report as the error caller
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
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
      cli::cli_abort(c("x" = "Node attribute {.str name} is missing from the graph"), call = call)
    }
    if (!"node_type" %in% vertex_attr_names(cellgraph)) {
      cli::cli_abort(c("x" = "Node attribute {.str node_type} is missing from the graph"), call = call)
    }
  }
  # TODO: Add check for A-node-projection and linegraph

  node_names <- .cg_node_names(cellgraph)
  if (anyDuplicated(node_names)) {
    cli::cli_abort(c("x" = "Node names in {.arg cellgraph} must be unique."), call = call)
  }
}

#' Validate node-level variable names across CellGraph data sources
#'
#' Prevents ambiguous lookups in \code{FetchData.CellGraph}. Names in the
#' graph node table, \code{meta.data}, and each reduction's embeddings must
#' be unique across all sources and must not overlap matrix feature names.
#' Count and layer matrices are grouped into one source because overlap
#' between those matrices is explicitly supported: callers disambiguate them
#' with the \code{layer} argument.
#'
#' Layout coordinate names are excluded because layouts are selected by
#' \code{layout_method} and are not searched by \code{FetchData.CellGraph}.
#'
#' @param object A current \code{CellGraph} object
#' @param call Environment to report as the error caller
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
#' @noRd
#'
.validate_cellgraph_data_names <- function(object, call = caller_env()) {
  .assert_current_cellgraph(object, call = call)
  graph <- slot(object, "cellgraph")
  graph_names <- if (is.null(graph)) {
    character()
  } else {
    igraph::vertex_attr_names(graph)
  }

  meta <- slot(object, "meta.data")
  meta_names <- if (is.null(meta)) character() else colnames(meta) %||% character()

  reductions <- slot(object, "reductions")
  reduction_sources <- list()
  for (reduction_name in names(reductions)) {
    reduction_sources[paste0("reduction '", reduction_name, "'")] <- list(
      colnames(Embeddings(reductions[[reduction_name]])) %||% character()
    )
  }

  counts <- slot(object, "counts")
  matrix_sources <- list()
  if (!is.null(counts)) {
    matrix_sources["counts"] <- list(colnames(counts) %||% character())
  }
  layers <- slot(object, "layers")
  for (layer_name in names(layers)) {
    matrix_sources[paste0("layer '", layer_name, "'")] <- list(
      colnames(layers[[layer_name]]) %||% character()
    )
  }
  for (source in names(matrix_sources)) {
    duplicated_names <- unique(
      matrix_sources[[source]][duplicated(matrix_sources[[source]])]
    )
    if (length(duplicated_names) > 0) {
      cli::cli_abort(
        c(
          "x" = "Feature names in {source} must be unique.",
          "i" = "Duplicated name{?s}: {.val {duplicated_names}}"
        ),
        call = call
      )
    }
  }
  matrix_names <- unique(unlist(matrix_sources, use.names = FALSE))

  sources <- c(
    list(
      "cellgraph node table" = graph_names,
      "meta.data" = meta_names
    ),
    reduction_sources,
    list("counts/layers" = unique(matrix_names))
  )
  sources <- sources[vapply(sources, length, integer(1)) > 0]

  exclusive_sources <- setdiff(names(sources), "counts/layers")
  for (source in exclusive_sources) {
    duplicated_names <- unique(sources[[source]][duplicated(sources[[source]])])
    if (length(duplicated_names) > 0) {
      cli::cli_abort(
        c(
          "x" = "Node-level variable names in {source} must be unique.",
          "i" = "Duplicated name{?s}: {.val {duplicated_names}}"
        ),
        call = call
      )
    }
  }

  if (length(sources) > 1) {
    source_pairs <- utils::combn(names(sources), 2, simplify = FALSE)
    for (pair in source_pairs) {
      overlap <- intersect(sources[[pair[[1]]]], sources[[pair[[2]]]])
      if (length(overlap) > 0) {
        source_x <- pair[[1]]
        source_y <- pair[[2]]
        cli::cli_abort(
          c(
            "x" = paste0(
              "{cli::qty(length(overlap))}Node-level variable name{?s} ",
              "{.val {overlap}} {?is/are} present in both {source_x} and {source_y}."
            ),
            "i" = paste0(
              "Names must be unique across the cellgraph node table, ",
              "meta.data, reductions, and matrix features."
            ),
            "i" = "Only counts and layer matrices may share feature names."
          ),
          call = call
        )
      }
    }
  }

  invisible(NULL)
}

#' Reject CellGraph objects created before the current slot layout
#'
#' Objects serialized by older versions only have \code{cellgraph},
#' \code{counts}, and \code{layout}. Reading \code{layers},
#' \code{meta.data}, or \code{reductions} on those instances fails with
#' \code{no slot of name ...}, which says nothing about the cause, so
#' check for the slots up front and report what to do instead.
#'
#' @param object A \code{CellGraph}, or another object (checked and ignored)
#' @param call Environment to report as the error caller
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
#' @noRd
#'
.assert_current_cellgraph <- function(object, call = caller_env()) {
  if (!is(object, "CellGraph")) {
    return(invisible(NULL))
  }
  added_slots <- c("layers", "meta.data", "reductions", "nodes")
  missing_slots <- added_slots[!vapply(added_slots, function(nm) {
    .hasSlot(object, nm)
  }, logical(1))]
  if (length(missing_slots) == 0) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "x" = "This {.cls CellGraph} has no {.field {missing_slots}} slot{?s}.",
      "i" = "It was saved by {.pkg pixelatorR} 0.21.0 or earlier, before
             {.cls CellGraph} gained these slots. Such objects are not upgraded.",
      "i" = "Load the cell graphs again from the PXL file with {.fn LoadCellGraphs},",
      " " = "or read the object with a version that still supports the old class:",
      " " = "{.code remotes::install_github(\"PixelgenTechnologies/pixelatorR@v0.20.1\")}"
    ),
    call = call
  )
}

#' Node names for a tbl_graph
#'
#' Reads the \code{name} vertex attribute when present; otherwise uses
#' \code{"1"}, \code{"2"}, ... in node order.
#'
#' @param cellgraph A \code{tbl_graph}
#'
#' @return A character vector of node names, one per vertex
#'
#' @keywords internal
#' @noRd
#'
.cg_node_names <- function(cellgraph) {
  if ("name" %in% vertex_attr_names(cellgraph)) {
    return(as.character(cellgraph %N>% pull(name)))
  }
  as.character(seq_along(cellgraph))
}

#' Ensure a tbl_graph has a name vertex attribute
#'
#' Graphs without \code{name} are identified by \code{"1"}, \code{"2"}, ...
#' in node order. That identity is written onto the graph so later subsetting
#' can filter on \code{name} and keep the original IDs after nodes are dropped.
#'
#' @param cellgraph A \code{tbl_graph}
#'
#' @return \code{cellgraph} with a \code{name} vertex attribute
#'
#' @keywords internal
#' @noRd
#'
.cg_graph_with_node_names <- function(cellgraph) {
  if ("name" %in% vertex_attr_names(cellgraph)) {
    return(cellgraph)
  }
  graph_type <- attr(cellgraph, "type")
  cellgraph <- cellgraph %N>% mutate(name = .cg_node_names(cellgraph))
  attr(cellgraph, "type") <- graph_type
  cellgraph
}

#' Match row names of a matrix to node names
#'
#' Reorders rows to \code{node_names}. If the matrix has no row names and
#' \code{nrow} matches the graph, names are assigned in current row order.
#' Duplicate or missing node names abort with an error attributed to
#' \code{call}.
#'
#' @param mat A matrix-like object, or \code{NULL}
#' @param node_names Character vector of graph node names (target row order)
#' @param arg Name of the argument to cite in error messages
#' @param call Environment to report as the error caller
#'
#' @return \code{mat} with rows in \code{node_names} order, or \code{NULL}
#'
#' @keywords internal
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
    return(.drop_row_ids(mat))
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
  .drop_row_ids(mat[node_names, , drop = FALSE])
}

#' Align the counts matrix
#'
#' Requires a \code{dgCMatrix} (when not \code{NULL}) and matches rows to
#' graph node names.
#'
#' @param counts A \code{dgCMatrix}, or \code{NULL}
#' @param node_names Character vector of graph node names
#' @param call Environment to report as the error caller
#'
#' @return \code{counts} aligned to \code{node_names}, or \code{NULL}
#'
#' @keywords internal
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
#' Coerces data frames to matrices, requires numeric values, and matches
#' rows to graph node names. Used for extra \code{layers} on a
#' \code{CellGraph}.
#'
#' @param mat A matrix, \code{Matrix}, or data frame
#' @param node_names Character vector of graph node names
#' @param arg Name of the argument to cite in error messages
#' @param call Environment to report as the error caller
#'
#' @return \code{mat} as a matrix aligned to \code{node_names}
#'
#' @keywords internal
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
#' Each element is aligned with \code{\link{.align_node_matrix}}. The name
#' \code{"counts"} is reserved for the \code{counts} slot.
#'
#' @param layers A named list of numeric node matrices, or \code{NULL}
#' @param node_names Character vector of graph node names
#' @param call Environment to report as the error caller
#'
#' @return A named list of aligned matrices, or an empty list
#'
#' @keywords internal
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
#' also end up with automatic row names once coerced). Explicit row names
#' are returned as a character vector.
#'
#' @param x A data frame, tibble, or object coercible to a data frame
#'
#' @return Character row names, or \code{NULL} if they are automatic
#'
#' @keywords internal
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

#' Strip the \code{-A}/\code{-B} suffix used on MPX bipartite graph nodes.
#'
#' Layout tables in PXL files identify nodes without that suffix. Same
#' pattern as \code{\link{LoadCellGraphs}} / \code{\link{WriteMPX_pxl_file}}.
#'
#' @param x Character node names
#'
#' @return \code{x} with the first \code{-A} or \code{-B} suffix removed
#'
#' @keywords internal
#' @noRd
#'
.strip_bipartite_node_suffix <- function(x) {
  stringr::str_replace(x, "-[A|B]", "")
}

#' Align a layout table to graph node order
#'
#' Accepts a data frame or matrix. Nodes are identified by row names or a
#' \code{name} column (which is then dropped so only coordinates remain).
#' If neither is present and \code{nrow} matches the graph, rows are assumed
#' to follow node order. When identifiers are present but do not match graph
#' node IDs, \code{-A}/\code{-B} suffixes are stripped from the graph names
#' so MPX bipartite layouts can share one row between the two partitions.
#' The result is a base \code{data.frame} in \code{node_names} order without
#' stored node IDs as row names.
#'
#' @param layout A data frame or matrix of coordinates
#' @param node_names Character vector of graph node names
#' @param layout_name Name of this layout, used in error messages
#' @param call Environment to report as the error caller
#'
#' @return A \code{data.frame} of coordinates in \code{node_names} order
#'
#' @keywords internal
#' @noRd
#'
.align_layout <- function(layout, node_names, layout_name = "layout", call = caller_env()) {
  # Matrix layouts are accepted so that coordinates produced by layout
  # functions can be stored directly; row names are kept as node identifiers
  if (is.matrix(layout) || inherits(layout, "Matrix")) {
    layout <- as.data.frame(as.matrix(layout))
  }
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
  layout_keys <- node_names
  missing_nodes <- setdiff(node_names, layout_names)
  if (length(missing_nodes) > 0) {
    stripped <- .strip_bipartite_node_suffix(node_names)
    if (all(stripped %in% layout_names)) {
      layout_keys <- stripped
    } else {
      cli::cli_abort(
        c(
          "x" = "The '{layout_name}' layout is missing {length(missing_nodes)} node{?s}.",
          "i" = "Example: {.val {head(missing_nodes, 3)}}"
        ),
        call = call
      )
    }
  }
  layout <- layout[match(layout_keys, layout_names), , drop = FALSE]
  .drop_row_ids(layout)
}

#' Align a named list of layouts
#'
#' Runs \code{\link{.align_layout}} on each named element. \code{NULL} is
#' returned unchanged (no layouts stored).
#'
#' @param layout A named list of layout tables, or \code{NULL}
#' @param node_names Character vector of graph node names
#' @param call Environment to report as the error caller
#'
#' @return A named list of aligned layout data frames, or \code{NULL}
#'
#' @keywords internal
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
#' Matches rows to graph nodes via explicit row names or a \code{name}
#' column. An empty table becomes a zero-column \code{data.frame} with
#' \code{node_names} as row names.
#'
#' @param meta A data frame or tibble, or \code{NULL}
#' @param node_names Character vector of graph node names
#' @param call Environment to report as the error caller
#'
#' @return A \code{data.frame} with row names \code{node_names}
#'
#' @keywords internal
#' @noRd
#'
.align_meta_data <- function(meta, node_names, call = caller_env()) {
  if (is.null(meta)) {
    return(.empty_node_meta(length(node_names)))
  }
  meta <- as.data.frame(meta, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(meta) == 0) {
    return(.empty_node_meta(length(node_names)))
  }
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
  .drop_row_ids(meta[node_names, , drop = FALSE])
}

#' Align node metadata, filling nodes that are not covered
#'
#' Like \code{.align_meta_data()}, but partial tables are allowed:
#' graph nodes without a row get \code{NA} and rows that do not name a
#' graph node are dropped. This follows
#' \code{\link[SeuratObject]{AddMetaData}}, which annotates a subset of
#' cells without touching the rest.
#'
#' @param meta A data frame or tibble
#' @param node_names Character vector of graph node names
#' @param call Environment to report as the error caller
#'
#' @return A \code{data.frame} with row names \code{node_names}
#'
#' @keywords internal
#' @noRd
#'
.fill_meta_data <- function(meta, node_names, call = caller_env()) {
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
          "x" = "{.arg metadata} has no node identifiers and {nrow(meta)} row{?s},",
          " " = "but the graph has {length(node_names)} node{?s}."
        ),
        call = call
      )
    }
    meta_names <- node_names
  }
  if (anyDuplicated(meta_names)) {
    cli::cli_abort(c("x" = "Node names in {.arg metadata} must be unique."), call = call)
  }
  if (length(intersect(node_names, meta_names)) == 0) {
    cli::cli_abort(
      c("x" = "No node in {.arg metadata} is present in this {.cls CellGraph}."),
      call = call
    )
  }
  filled <- meta[match(node_names, meta_names), , drop = FALSE]
  .drop_row_ids(filled)
}

#' Align a named list of NodeDimReduc objects
#'
#' Each reduction's embeddings are matched to graph node names.
#'
#' @param reductions A named list of \code{NodeDimReduc} objects, or \code{NULL}
#' @param node_names Character vector of graph node names
#' @param call Environment to report as the error caller
#'
#' @return A named list of aligned \code{NodeDimReduc} objects, or an empty list
#'
#' @keywords internal
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

#' Bind data.frames while keeping non-syntactic column names
#'
#' \code{cbind()} rewrites names such as \code{HLA-DR} to \code{HLA.DR}.
#' This helper cbinds, then restores the original column names.
#'
#' @param ... Data frames to bind column-wise
#'
#' @return A data frame, a single input if only one has columns, or \code{NULL}
#'
#' @keywords internal
#' @noRd
#'
.cbind_keep_names <- function(...) {
  dfs <- list(...)
  dfs <- dfs[vapply(dfs, function(d) {
    is.data.frame(d) && ncol(d) > 0
  }, logical(1))]
  if (length(dfs) == 0) {
    return(NULL)
  }
  if (length(dfs) == 1) {
    return(dfs[[1]])
  }
  kept_names <- unlist(lapply(dfs, names), use.names = FALSE)
  out <- do.call(cbind, dfs)
  names(out) <- kept_names
  out
}

#' Bind fetched columns onto a node-level data.frame
#'
#' Pads missing rows with \code{NA}, aligns by row name, and cbinds while
#' keeping non-syntactic names. Used by \code{FetchData.CellGraph}.
#'
#' @param data_fetched Accumulated result with nodes as row names
#' @param new_df Columns to add, with nodes as row names
#'
#' @return \code{data_fetched} with columns from \code{new_df} appended
#'
#' @keywords internal
#' @noRd
#'
.add_fetched_cols <- function(data_fetched, new_df) {
  if (is.null(new_df) || ncol(new_df) == 0) {
    return(data_fetched)
  }
  missing_rows <- setdiff(rownames(data_fetched), rownames(new_df))
  if (length(missing_rows) > 0) {
    pad <- new_df[rep(NA_integer_, length(missing_rows)), , drop = FALSE]
    rownames(pad) <- missing_rows
    new_df <- rbind(new_df, pad)
  }
  new_df <- new_df[rownames(data_fetched), , drop = FALSE]
  if (ncol(data_fetched) == 0) {
    return(new_df)
  }
  .cbind_keep_names(data_fetched, new_df)
}

#' Fetch marker columns from a CellGraph layer
#'
#' Pulls requested feature names from \code{LayerData()}. Columns that
#' match \code{vars} are returned without passing marker names through
#' \code{check.names}. Cross-source name collisions are rejected before
#' this helper is called.
#'
#' @param object A \code{CellGraph}
#' @param layer Layer name (for example \code{"counts"})
#' @param vars Character vector of requested variable names
#' @param cells Node names to keep as rows
#'
#' @return A data frame of selected columns, or \code{NULL} if none match
#'
#' @keywords internal
#' @noRd
#'
.fetch_layer_vars <- function(object, layer, vars, cells) {
  mat <- LayerData(object, layer = layer)
  if (is.null(mat) || ncol(mat) == 0) {
    return(NULL)
  }
  feature_vars <- intersect(vars, colnames(mat))
  if (length(feature_vars) == 0) {
    return(NULL)
  }
  as.data.frame(
    as.matrix(mat[.row_index(cells, rownames(mat)), feature_vars, drop = FALSE]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Search remaining vars in layers other than the default
#'
#' Features found in exactly one non-default layer are included (with a
#' warning). Features present in more than one extra layer are skipped.
#'
#' @param object A \code{CellGraph}
#' @param vars Character vector of names not found in the default layer
#' @param cells Node names to keep as rows
#' @param other_layers Character vector of layer names to search
#'
#' @return A data frame of recovered columns, or \code{NULL}
#'
#' @keywords internal
#' @noRd
#'
.fetch_vars_from_other_layers <- function(object, vars, cells, other_layers) {
  vars_alt <- vector("list", length(vars))
  names(vars_alt) <- vars
  for (lyr in other_layers) {
    mat <- LayerData(object, layer = lyr)
    if (is.null(mat) || ncol(mat) == 0) {
      next
    }
    for (var in intersect(vars, colnames(mat))) {
      vars_alt[[var]] <- c(vars_alt[[var]], lyr)
    }
  }
  n_hits <- vapply(vars_alt, length, integer(1))
  vars_many <- names(vars_alt)[n_hits > 1]
  if (length(vars_many) > 0) {
    cli::cli_warn(
      "Found the following features in more than one layer besides the default;
      they will not be included: {.val {vars_many}}"
    )
  }
  vars_one <- vars_alt[n_hits == 1]
  if (length(vars_one) == 0) {
    return(NULL)
  }
  pieces <- lapply(names(vars_one), function(var) {
    lyr <- vars_one[[var]]
    cli::cli_warn("Could not find {.val {var}} in the default layer, found in {.val {lyr}} instead")
    .fetch_layer_vars(object, lyr, var, cells)
  })
  Reduce(function(left, right) {
    if (is.null(left)) {
      return(right)
    }
    if (is.null(right)) {
      return(left)
    }
    .add_fetched_cols(left, right)
  }, pieces)
}

#' Vertex attributes of a CellGraph as a node-level data.frame
#'
#' Collects atomic, unnamed vertex attributes (such as \code{name} and
#' \code{node_type}) into a data frame keyed by node name for
#' \code{FetchData.CellGraph}.
#'
#' @param cellgraph A \code{tbl_graph}
#'
#' @return A data frame with one row per node
#'
#' @keywords internal
#' @noRd
#'
.cg_vertex_attr_df <- function(cellgraph) {
  node_names <- .cg_node_names(cellgraph)
  attr_names <- igraph::vertex_attr_names(cellgraph)
  if (length(attr_names) == 0) {
    return(data.frame(row.names = node_names))
  }
  attrs <- lapply(attr_names, function(nm) {
    igraph::vertex_attr(cellgraph, name = nm)
  })
  names(attrs) <- attr_names
  keep <- vapply(attrs, function(x) {
    is.atomic(x) && is.null(dim(x)) && length(x) == length(node_names)
  }, logical(1))
  if (!any(keep)) {
    return(data.frame(row.names = node_names))
  }
  df <- as.data.frame(attrs[keep], stringsAsFactors = FALSE, check.names = FALSE)
  rownames(df) <- node_names
  df
}

#' Fetch embedding columns from a NodeDimReduc
#'
#' Selects columns whose names match \code{vars}, including keyed names
#' such as \code{PC_1}.
#'
#' @param object A \code{NodeDimReduc}
#' @param vars Character vector of requested variable names
#' @param cells Node names to keep as rows
#' @param node_names Node IDs in embedding row order
#'
#' @return A data frame of embedding columns, or \code{NULL} if none match
#'
#' @keywords internal
#' @noRd
#'
.fetch_nodedimreduc_vars <- function(object, vars, cells, node_names) {
  key <- Key(object)
  emb <- Embeddings(object)
  if (is.null(emb) || ncol(emb) == 0) {
    return(NULL)
  }
  keyed <- character()
  if (length(key) && nzchar(key)) {
    keyed <- grep(paste0("^", key), vars, value = TRUE)
  }
  keyed <- unique(c(keyed, intersect(vars, colnames(emb))))
  if (length(keyed) == 0) {
    return(NULL)
  }
  missing <- setdiff(keyed, colnames(emb))
  keyed <- setdiff(keyed, missing)
  if (length(keyed) == 0) {
    return(NULL)
  }
  if (!is.null(rownames(emb))) {
    node_names <- as.character(rownames(emb))
  }
  if (length(node_names) != nrow(emb)) {
    return(NULL)
  }
  cells_keep <- intersect(cells, node_names)
  if (length(cells_keep) == 0) {
    return(NULL)
  }
  as.data.frame(
    .with_row_ids(emb[.row_index(cells_keep, node_names), keyed, drop = FALSE], cells_keep),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Fetch a named reduction from a CellGraph
#'
#' When \code{reduction} is \code{NULL}, the first stored reduction is used.
#'
#' @param object A \code{CellGraph}
#' @param reduction Name of a stored \code{NodeDimReduc}, or \code{NULL}
#' @param call Environment to report as the error caller
#'
#' @return A \code{NodeDimReduc} object
#'
#' @keywords internal
#' @noRd
#'
.get_cellgraph_reduction <- function(object, reduction = NULL, call = caller_env()) {
  .assert_current_cellgraph(object, call = call)
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
#' Used before subsetting. Layouts without row names (or with a \code{name}
#' column) are converted to node-ID row names. Count, layer, and metadata
#' matrices that lack row names but match the graph length get node names
#' assigned in current order.
#'
#' @param object A \code{CellGraph}
#'
#' @return \code{object} with node IDs on slots that were positional
#'
#' @keywords internal
#' @noRd
#'
.ensure_node_ids_on_slots <- function(object) {
  graph_type <- attr(slot(object, "cellgraph"), "type")
  slot(object, "cellgraph") <- .cg_graph_with_node_names(slot(object, "cellgraph"))
  attr(slot(object, "cellgraph"), "type") <- graph_type
  node_names <- .cg_node_names(slot(object, "cellgraph"))
  if (length(slot(object, "nodes")) == 0) {
    slot(object, "nodes") <- node_names
  }

  counts <- slot(object, "counts")
  if (!is.null(counts) && !is.null(rownames(counts))) {
    slot(object, "counts") <- .align_counts(counts, node_names)
  }

  layouts <- slot(object, "layout")
  if (!is.null(layouts) && length(layouts) > 0) {
    layout_names <- names(layouts)
    slot(object, "layout") <- lapply(seq_along(layouts), function(i) {
      ly <- layouts[[i]]
      nm <- if (!is.null(layout_names) && nzchar(layout_names[[i]])) {
        layout_names[[i]]
      } else {
        "layout"
      }
      if (!is.null(.explicit_rownames(ly)) || "name" %in% colnames(ly)) {
        ly <- .align_layout(ly, node_names, layout_name = nm)
      }
      ly
    })
    names(slot(object, "layout")) <- layout_names
  }

  meta <- slot(object, "meta.data")
  if (ncol(meta) == 0) {
    slot(object, "meta.data") <- .empty_node_meta(length(node_names))
  } else if (!is.null(.explicit_rownames(meta))) {
    slot(object, "meta.data") <- .align_meta_data(meta, node_names)
  }

  layers <- slot(object, "layers")
  if (length(layers) > 0) {
    slot(object, "layers") <- lapply(names(layers), function(nm) {
      mat <- layers[[nm]]
      if (!is.null(rownames(mat))) {
        mat <- .align_node_matrix(mat, node_names, arg = nm)
      }
      mat
    })
    names(slot(object, "layers")) <- names(layers)
  }

  object
}

#' Reorder or subset all node-level slots to \code{node_names}
#'
#' After the graph is filtered or replaced, counts, layouts, layers,
#' metadata, and reductions are aligned to the remaining node names so
#' they stay in sync with the graph. Stored tables are subset by the
#' central \code{nodes} map, not by copied row names.
#'
#' @param object A \code{CellGraph}
#' @param node_names Character vector of node names to keep, in graph order
#' @param call Environment to report as the error caller
#'
#' @return \code{object} with every node-level slot aligned to \code{node_names}
#'
#' @keywords internal
#' @noRd
#'
.remap_cellgraph_nodes <- function(object, node_names, call = caller_env()) {
  old_nodes <- .cg_node_map(object)
  idx <- match(node_names, old_nodes)
  if (anyNA(idx)) {
    cli::cli_abort(
      c("x" = "Cannot remap nodes that are missing from the {.field nodes} map."),
      call = call
    )
  }

  counts <- slot(object, "counts")
  if (!is.null(counts)) {
    if (!is.null(rownames(counts))) {
      slot(object, "counts") <- .align_counts(counts, node_names, call = call)
    } else {
      slot(object, "counts") <- .drop_row_ids(counts[idx, , drop = FALSE])
    }
  }

  layouts <- slot(object, "layout")
  if (!is.null(layouts) && length(layouts) > 0) {
    if (any(vapply(layouts, function(ly) {
      !is.null(.explicit_rownames(ly)) || "name" %in% colnames(ly)
    }, logical(1)))) {
      slot(object, "layout") <- .align_layout_list(layouts, node_names, call = call)
    } else {
      slot(object, "layout") <- lapply(layouts, function(ly) {
        .drop_row_ids(ly[idx, , drop = FALSE])
      })
    }
  }

  layers <- slot(object, "layers")
  if (length(layers) > 0) {
    slot(object, "layers") <- lapply(names(layers), function(nm) {
      mat <- layers[[nm]]
      if (!is.null(rownames(mat))) {
        .align_node_matrix(mat, node_names, arg = nm, call = call)
      } else {
        .drop_row_ids(mat[idx, , drop = FALSE])
      }
    })
    names(slot(object, "layers")) <- names(layers)
  }

  meta <- slot(object, "meta.data")
  if (ncol(meta) == 0) {
    slot(object, "meta.data") <- .empty_node_meta(length(node_names))
  } else if (!is.null(.explicit_rownames(meta))) {
    slot(object, "meta.data") <- .align_meta_data(meta, node_names, call = call)
  } else {
    slot(object, "meta.data") <- .drop_row_ids(meta[idx, , drop = FALSE])
  }

  reductions <- slot(object, "reductions")
  if (length(reductions) > 0) {
    slot(object, "reductions") <- lapply(names(reductions), function(nm) {
      dr <- reductions[[nm]]
      emb <- slot(dr, "embeddings")
      if (!is.null(rownames(emb))) {
        .align_node_dimreduc(dr, node_names, arg = nm, call = call)
      } else {
        slot(dr, "embeddings") <- .drop_row_ids(emb[idx, , drop = FALSE])
        dr
      }
    })
    names(slot(object, "reductions")) <- names(reductions)
  }

  slot(object, "nodes") <- node_names
  object
}
