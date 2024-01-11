# Declarations used in package check
globalVariables(
  names = c('nodes', 'group'),
  package = 'pixelatorR',
  add = TRUE
)


#' @param verbose Print messages
#'
#' @importFrom tidygraph `%N>%`
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent tbl_graph
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
    filtered_object <- object %N>%
      mutate(group = group_components()) %>%
      filter(group == 1)
    if (verbose && check_global_verbosity())
      cli_alert_info("Removed {length(object) - length(filtered_object)} out of {length(object)} nodes")
  }
  return(filtered_object)
}


#' @param verbose Print messages
#'
#' @rdname KeepLargestComponent
#' @method KeepLargestComponent CellGraph
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

  # Fix counts
  counts <- object@counts
  counts <- counts[match(node_ids_filtered, node_ids), ]
  slot(object, name = "counts") <- counts

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
