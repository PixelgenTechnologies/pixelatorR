#' Plot 2D graph layouts
#'
#' Plot 2D component graph layouts computed with \code{\link{ComputeLayout}} and
#' optionally color nodes by certain attributes. Edges can also be visualized
#' by setting \code{map_edges}; however, since component graphs tend to be very
#' large, this can take a long time to draw.
#'
#' @param object A \code{Seurat} object
#' @param cells A character vector with cell IDs
#' @param marker Name of a marker to colors nodes/edges by
#' @param assay Name of assay to pull data from
#' @param layout_method Select appropriate layout previously computed with
#' \code{\link{ComputeLayout}}
#' @param colors A character vector of colors to color marker counts by
#' @param map_nodes,map_edges Should nodes and/or edges be mapped?
#' @param node_size Size of nodes
#' @param edge_width Not yet implemented TODO
#' @param ... Not yet implemented
#'
#' @rdname Plot2DGraph
#'
#' @import rlang
#' @import glue
#' @import ggplot2
#' @import patchwork
#' @importFrom ggraph ggraph geom_node_point geom_edge_link
#' @importFrom tidygraph `%N>%`
#'
#' @return An object of class \code{patchwork}
#'
#' @examples
#' library(pixelatorR)
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
#'                         package = "pixelatorR")
#'
#' seur <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
#' seur <- LoadCellGraphs(seur, load_as = "Anode", cells = colnames(seur)[1:10])
#' seur[["mpxCells"]] <- KeepLargestComponent(seur[["mpxCells"]])
#' seur <- ComputeLayout(seur, layout_method = "pmds", dim = 2)
#'
#' Plot2DGraph(seur, cells = colnames(seur)[1], marker = "HLA-ABC")
#'
#' @export
#'
Plot2DGraph <- function (
  object,
  cells = NULL,
  marker = NULL,
  assay = NULL,
  layout_method = c("pmds", "wpmds", "fr", "kk", "drl"),
  colors = c("lightgrey", "mistyrose", "red", "darkred"),
  map_nodes = TRUE,
  map_edges = FALSE,
  node_size = 0.5,
  edge_width = 0.3,
  ...
) {

  # Validate input parameters
  stopifnot(
    "'colors' must be a character vector with at least 2 colors" =
      is.character(colors) &&
      (length(colors) > 1),
    "'map_nodes' must be one of TRUE or FALSE" =
      is.logical(map_nodes) &&
      (length(map_nodes) == 1),
    "'map_edges' must be one of TRUE or FALSE" =
      is.logical(map_edges) &&
      (length(map_edges) == 1),
    "'map_nodes' and 'map_edges' cannot be deactivated at the same time" =
      map_nodes ||
      map_edges,
    "'cells' must be a non-empty character vector with cell IDs" =
      is.character(cells) &&
      (length(cells) > 0)
  )

  if (!is.null(marker)) {
    stopifnot(
      "'marker' must be a character of length 1" =
        is.character(marker) &&
        (length(marker) == 1)
    )
  }

  # Check and select a layout method
  layout_method <- match.arg(layout_method, choices = c("pmds", "wpmds", "fr", "kk", "drl"))
  layout_method_ext <- switch (layout_method,
                           "fr" = "Fruchterman Reingold (fr)",
                           "kk" = "Kamada Kawai (kk)",
                           "drl" = "DrL graph layout generator (drl)",
                           "pmds" = "pivot MDS (pmds)"
  )

  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    stopifnot(
      "'assay' must be a character of length 1" =
        is.character(assay) &&
        (length(assay) == 1)
    )
  } else {
    assay <- DefaultAssay(object)
  }

  # Validate assay
  cg_assay <- object[[assay]]
  if (!inherits(cg_assay, what = "CellGraphAssay")) {
    abort(glue("Invalid assay type '{class(cg_assay)}'. Expected a 'CellGraphAssay'"))
  }

  # Fetch data
  data_list <- lapply(cells, function(cell_id) {

    # Fetch component graph
    component_graph <- CellGraphs(cg_assay)[[cell_id]]
    if (is.null(component_graph))
      abort(glue("Missing cellgraph for component '{cell_id}'"))

    # unpack values
    graph <- component_graph@cellgraph
    layout <- component_graph@layout[[layout_method]]
    if (length(graph) == 0)
      abort(glue("Missing cellgraph for component '{cell_id}'"))
    if (length(layout) == 0)
      abort(glue("Missing layout '{layout_method}' for component '{cell_id}'"))

    # Add node marker counts if needed
    if (!is.null(marker)) {
      if (marker != "node_type") {
          graph <- graph %N>%
            mutate(marker = component_graph@counts[, marker]) %>%
            arrange(marker)
      }
    }

    data <- list(graph = graph, layout = layout, type = attr(graph, "type"), layout_type = layout_method)
    return(data)
  }) %>% setNames(nm = cells)

  # Create plots
  plots <- lapply(names(data_list), function(nm) {
    # Visualize the graph with ggraph
    p <- data_list[[nm]]$graph %>%
      ggraph(layout = data_list[[nm]]$layout) +
      {
        # Add edges if map_edges=TRUE
        if (map_edges) {
          geom_edge_link()
        }
      } +
      {
        # Add nodes if map_nodes=TRUE
        if (map_nodes) {
          geom_node_point(mapping = aes(color = marker), size = node_size)
        }
      } +
      coord_fixed() +
      labs(title = glue("{nm}"), color = marker) +
      theme_void() +
      theme(plot.title = element_text(size = 10))
  })

  # Wrap plots
  p <- wrap_plots(plots)
  if (!is.null(marker)) {
    if (marker != "node_type") {
      p <- p & scale_color_gradientn(colours = colors)
    }
  }
  p <- p + plot_annotation(title = glue("Layout with {layout_method_ext}"),
                           theme = theme(plot.title = element_text(size = 14, face = "bold")))

  return(p)
}
