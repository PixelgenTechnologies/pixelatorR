#' @param verbose Print messages
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent tbl_graph
#'
#' @examples
#' library(pixelatorR)
#' library(tidygraph)
#'
#' pxl_file <- minimal_mpx_pxl_file()
#'
#' # Read edgelist
#' edgelist <- ReadMPX_arrow_edgelist(pxl_file)
#'
#' # Load graph from edge list and store in a CellGraph object
#' cg <- LoadCellGraphs(edgelist, cells = "RCVCMP0000217", data_type = "MPX")[[1]]
#' cg
#'
#' # Fetch tbl_graph from CellGraph object
#' g <- CellGraphData(cg, slot = "cellgraph")
#' g
#'
#' # Break graph by removing random edges
#' set.seed(132)
#' g <- g %E>%
#'   filter(from %in% sample(from, n() - 500))
#'
#' # Fetch largest component from a tbl_graph
#' g_largest <- KeepLargestComponent(g)
#' g_largest
#'
#' @export
#'
KeepLargestComponent.tbl_graph <- function(
  object,
  verbose = TRUE,
  ...
) {
  if (igraph::is_connected(object)) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("Graph is already connected")
    }
    return(object)
  } else {
    graph_type_attribute <- attr(object, "type")
    split_components <- object %>%
      to_components()
    largest_component <- which.max(sapply(split_components, length))
    object_largest <- split_components[[largest_component]]
    attr(object_largest, "type") <- graph_type_attribute
    if (verbose && check_global_verbosity()) {
      cli_alert_info("Removed {length(object) - length(object_largest)} out of {length(object)} nodes")
    }
  }
  return(object_largest)
}


#' @param verbose Print messages
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent CellGraph
#'
#' @examples
#' # Add new graph to CellGraph object
#' cg@cellgraph <- g_largest
#'
#' # Fetch largest component from a CellGraph
#' cg_largest <- KeepLargestComponent(cg)
#' cg_largest
#'
#' @export
#'
KeepLargestComponent.CellGraph <- function(
  object,
  verbose = TRUE,
  ...
) {
  # Fetch node ids before filtering
  node_ids <- object@cellgraph %N>% pull(name)

  filtered_graph <- KeepLargestComponent(slot(object, name = "cellgraph"), verbose = verbose, ...)
  slot(object, name = "cellgraph") <- filtered_graph

  # fetch node ids after filtering
  node_ids_filtered <- filtered_graph %N>% pull(name)

  # Filter counts
  counts <- slot(object, name = "counts")
  if (!is.null(counts)) {
    counts <- counts[match(node_ids_filtered, node_ids), ]
    slot(object, name = "counts") <- counts
  }

  # Filter layouts if available
  layouts <- slot(object, name = "layout")
  if (!is.null(layouts)) {
    layouts <- lapply(layouts, function(layout) {
      layout[match(node_ids_filtered, node_ids), ]
    })
    slot(object, name = "layout") <- layouts
  }

  return(object)
}


#' @param verbose Print messages
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent MPXAssay
#'
#' @export
#'
KeepLargestComponent.MPXAssay <- function(
  object,
  verbose = TRUE,
  ...
) {
  # Check cellgraphs
  cellgraphs <- slot(object, name = "cellgraphs")
  loaded_graphs <- !sapply(cellgraphs, is.null)

  if (sum(loaded_graphs) == 0) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("No 'cellgraphs' loaded in object. Returning object unmodified.")
    }
    return(object)
  }

  # Only keep loaded graphs
  cellgraphs_loaded <- cellgraphs[loaded_graphs]

  if (verbose && check_global_verbosity()) {
    cli_alert_info("Keeping largest component for {length(cellgraphs_loaded)} graphs")
  }

  cellgraphs_loaded <- lapply(cellgraphs_loaded, function(g) {
    KeepLargestComponent(g, verbose = verbose, ...)
  })
  slot(object, name = "cellgraphs")[names(cellgraphs_loaded)] <- cellgraphs_loaded

  return(object)
}

#' @param verbose Print messages
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent PNAAssay
#'
#' @export
#'
KeepLargestComponent.PNAAssay <- KeepLargestComponent.MPXAssay

#' @param verbose Print messages
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent PNAAssay5
#'
#' @export
#'
KeepLargestComponent.PNAAssay5 <- KeepLargestComponent.MPXAssay

#' @param assay Name of a \code{CellGraphAssay}, \code{CellGraphAssay5},
#' \code{PNAAssay} or\code{PNAAssay5}
#' assay stored on the \code{Seurat} object
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent Seurat
#'
#' @export
#'
KeepLargestComponent.Seurat <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  cg_assay <- object[[assay]]
  assert_pixel_assay(cg_assay)

  cg_assay <- KeepLargestComponent(cg_assay, verbose = verbose, ...)
  object[[assay]] <- cg_assay

  return(object)
}
