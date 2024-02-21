# Declarations used in package check
globalVariables(
  names = c('nodes', 'group'),
  package = 'pixelatorR',
  add = TRUE
)


#' @param verbose Print messages
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent tbl_graph
#'
#' @examples
#' library(pixelatorR)
#' # Set arrow data output directory to temp for tests
#' options(pixelatorR.arrow_outdir = tempdir())
#'
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
#'                         package = "pixelatorR")
#'
#' # Read edgelist
#' edgelist <- ReadMPX_arrow_edgelist(pxl_file, overwrite = TRUE)
#'
#' # Load graph from edge list and store in a CellGraph object
#' cg <- LoadCellGraphs(edgelist, cells = "RCVCMP0000000")[[1]]
#' cg
#'
#' # Fetch tbl_graph from CellGraph object
#' g <- CellGraphData(cg, slot = "cellgraph")
#' g
#'
#' # Fetch largest component from a tbl_graph
#' g_largest <- KeepLargestComponent(g)
#' g_largest
#'
#' @export
#'
KeepLargestComponent.tbl_graph <- function (
  object,
  verbose = TRUE,
  ...
) {
  if (igraph::is_connected(object)) {
    if (verbose && check_global_verbosity())
      cli_alert_info("Graph is already connected")
    return(object)
  } else {
    graph_type_attribute <- attr(object, "type")
    split_components <- object %>%
      to_components()
    largest_component <- which.max(sapply(split_components, length))
    object_largest <- split_components[[largest_component]]
    attr(object_largest, "type") <- graph_type_attribute
    if (verbose && check_global_verbosity())
      cli_alert_info("Removed {length(object) - length(object_largest)} out of {length(object)} nodes")
  }
  return(object_largest)
}


#' @param verbose Print messages
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent CellGraph
#'
#' @examples
#' # Fetch largest component from a CellGraph
#' cg_largest <- KeepLargestComponent(cg)
#' cg_largest
#'
#' @export
#'
KeepLargestComponent.CellGraph <- function (
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
#' @method KeepLargestComponent CellGraphAssay
#'
#' @export
#'
KeepLargestComponent.CellGraphAssay <- function (
  object,
  verbose = TRUE,
  ...
) {

  # Check cellgraphs
  cellgraphs <- slot(object, name = "cellgraphs")
  loaded_graphs <- !sapply(cellgraphs, is.null)

  if (sum(loaded_graphs) == 0) {
    if (verbose && check_global_verbosity())
      cli_alert_info("No 'cellgraphs' loaded in 'CellGraphAssay'. Returning unmodified 'CellGraphAssay'.")
    return(object)
  }

  # Only keep loaded graphs
  cellgraphs_loaded <- cellgraphs[loaded_graphs]

  if (verbose && check_global_verbosity())
    cli_alert_info("Keeping largest component for {length(cellgraphs_loaded)} graphs")

  cellgraphs_loaded <- lapply(cellgraphs_loaded, function(g) {
    KeepLargestComponent(g, verbose = verbose, ...)
  })
  slot(object, name = "cellgraphs")[names(cellgraphs_loaded)] <- cellgraphs_loaded

  return(object)
}


#' @param assay Name of a \code{CellGraphAssay} stored on the \code{Seurat} object
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent Seurat
#'
#' @export
#'
KeepLargestComponent.Seurat <- function (
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  cg_assay <- object[[assay]]
  if (!inherits(cg_assay, what = "CellGraphAssay")) {
    abort(glue("assay '{assay}' is not a 'CellGraphAssay'"))
  }

  cg_assay <- KeepLargestComponent(cg_assay, verbose = verbose, ...)
  object[[assay]] <- cg_assay

  return(object)
}
