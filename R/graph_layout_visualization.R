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
#' @param log_scale Convert node counts to log-scale with \code{logp}
#' @param node_size Size of nodes
#' @param edge_width Set the width of the edges if \code{map_edges = TRUE}
#' @param show_Bnodes Should B nodes be included in the visualization?
#' This option is only applicable to bipartite graphs. Note that by removing
#' the B nodes, all edges are removed from the graph and hence, \code{map_edges}
#' will have no effect.
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
  cells,
  marker = NULL,
  assay = NULL,
  layout_method = c("pmds", "wpmds", "fr", "kk", "drl"),
  colors = c("lightgrey", "mistyrose", "red", "darkred"),
  map_nodes = TRUE,
  map_edges = FALSE,
  log_scale = TRUE,
  node_size = 0.5,
  edge_width = 0.3,
  show_Bnodes = FALSE,
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

    # Validate marker
    if (!is.null(marker)) {
      if (marker == "node_type") {
        stopifnot(
          "marker = 'node_type' can only be used for bipartite graphs" =
            attr(graph, "type") == "bipartite"
        )
      } else {
        if (!marker %in% colnames(component_graph@counts)) {
          abort(glue("'{marker}' is missing from node count matrix ",
                     "for component {cell_id}"))
        }
      }
    }

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
            {
              if (log_scale) {
                mutate(., marker = log1p(marker))
              } else {
                .
              }
            }
      }
    }

    # Remove B nodes if show_Bnodes=FALSE
    if ((attr(graph, "type") == "bipartite") && !show_Bnodes) {
      inds_keep <- (graph %>% pull(node_type)) == "A"
      graph <- graph %>%
        filter(node_type == "A")
      layout <- layout[inds_keep, ]
    }

    # Rearrange by marker
    if (!is.null(marker)) {
      if (marker != "node_type") {
        # Rearrange layout
        order <- order(graph %>% pull(marker))
        layout <- data.frame(layout, row.names = graph %>% pull(name))
        graph <- graph %>%
          arrange(marker)
        layout <- layout[order, ] %>% as_tibble()
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
          geom_edge_link(edge_width = edge_width)
        }
      } +
      {
        # Add nodes if map_nodes=TRUE
        if (map_nodes) {
          if (!is.null(marker)) {
            if (marker == "node_type") {
              geom_node_point(mapping = aes(color = node_type), size = node_size)
            } else {
              geom_node_point(mapping = aes(color = marker), size = node_size)
            }
          } else {
            geom_node_point(size = node_size)
          }
        }
      } +
      coord_fixed() +
      {
        if (!is.null(marker)) {
          if (marker != "node_type") {
            if (log_scale) {
              labs(title = glue("{nm}"), color = paste0(marker, "\n(log-scaled)"))
            } else {
              labs(title = glue("{nm}"), color = marker)
            }
          } else {
            labs(title = glue("{nm}"))
          }
        }
      } +
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


#' Plot 3D graph layouts
#'
#' Plot a 3D component graph layout computed with \code{\link{ComputeLayout}} and
#' color nodes by a marker. 
#' 
#' @param object  A \code{Seurat} object
#' @param cell_id  ID of component to visualize
#' @param marker Name of marker to color the nodes by
#' @param assay Name of assay to pull data from
#' @param layout_method  Select appropriate layout previously computed with
#' \code{\link{ComputeLayout}}
#' @param project Project the nodes onto a sphere. Default FALSE
#' @param aspectmode Set aspect ratio to one of "data", "auto" or "cube". 
#' If "cube", this scene's axes are drawn as a cube, regardless of the axes' ranges. 
#' If "data", this scene's axes are drawn in proportion with the axes' ranges. 
#' If "auto", this scene's axes are drawn using the results of "data" except when one axis is more than four times the size of the two others, where in that case the results of "cube" are used. 
#' 
#' Default "data"
#' 
#' @param color Color the nodes expressing a marker. Default "darkred"
#' @param log_scale  Convert node counts to log-scale with \code{logp}
#' @param node_size Size of nodes
#' @param show_Bnodes  Should B nodes be included in the visualization?
#' This option is only applicable to bipartite graphs.
#' @param ... Additional parameters
#' @param showgrid Show the grid lines. Default TRUE
#' 
#' @rdname Plot3DGraph
#' 
#' @import plotly
#' @import glue
#' @import ggplot2
#'
#' @return A interactive 3D plot of a component graph layout as a \code{plotly} object
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
#' seur <- ComputeLayout(seur, layout_method = "pmds", dim = 3)
#'
#' Plot3DGraph(seur, cells = colnames(seur)[1], marker = "HLA-ABC")
#'
#' @export
Plot3DGraph <- function (
    object,
    cell_id,
    marker = NULL,
    assay = NULL,
    layout_method = c("pmds", "wpmds", "fr", "kk", "drl"),
    project = FALSE,
    aspectmode = c("data", "auto"),
    color = "darkred",
    showgrid = TRUE,
    log_scale = TRUE,
    node_size = 2,
    show_Bnodes = FALSE,
    ...
) {
  
  # Validate input parameters
  stopifnot(
    "'color' must be a character vector with a single color" =
      is.character(color) &&
      (length(color) == 1),
    "'cell_id' must be a non-empty character vector with a single cell ID" =
      is.character(cell_id) &&
      (length(cell_id) == 1)
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
  
  # Check and select an aspectmode
  aspectmode <- match.arg(aspectmode)
  
  
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
  
  
  # Fetch component graph
  component_graph <- CellGraphs(cg_assay)[[cell_id]]
  
  if (is.null(component_graph)) abort(glue("Missing cellgraph for component '{cell_id}'"))
  
  # unpack values
  graph <- component_graph@cellgraph
  
  # Validate marker
  if (!is.null(marker)) {
    if (marker == "node_type") {
      stopifnot(
        "marker = 'node_type' can only be used for bipartite graphs" =
          attr(graph, "type") == "bipartite"
      )
    } else {
      if (!marker %in% colnames(component_graph@counts)) {
        abort(glue("'{marker}' is missing from node count matrix ",
                   "for component {cell_id}"))
      }
    }}
    
    layout <- component_graph@layout[[layout_method]]
    if (length(graph) == 0)
      abort(glue("Missing cellgraph for component '{cell_id}'"))
    if (length(layout) == 0)
      abort(glue("Missing layout '{layout_method}' for component '{cell_id}'"))
    if (length(layout) < 3)
      abort(glue("Too few dimensions for a 3D visualization of layout '{layout_method}' for component '{cell_id}'"))
    
    # Add node marker counts if needed
    if (!is.null(marker)) {
      if (marker != "node_type") {
        graph <- graph %N>%
          mutate(marker = component_graph@counts[, marker]) %>%
          {
            if (log_scale) {
              mutate(., marker = log1p(marker))
            } else {
              .
            }
          }
      }
    }
    
    # Remove B nodes if show_Bnodes=FALSE
    if ((attr(graph, "type") == "bipartite") && !show_Bnodes) {
      inds_keep <- (graph %>% pull(node_type)) == "A"
      graph <- graph %>%
        filter(node_type == "A")
      layout <- layout[inds_keep, ]
    }
    
    # Rearrange by marker
    if (!is.null(marker)) {
      if (marker != "node_type") {
        # Rearrange layout
        order <- order(graph %>% pull(marker))
        layout <- data.frame(layout, row.names = graph %>% pull(name))
        graph <- graph %>%
          arrange(marker)
        layout <- layout[order, ] %>% as_tibble()
      }
    }
    
    data_list <- list(graph = graph, layout = layout, type = attr(graph, "type"), layout_type = layout_method)
  
  # Create plot
    plot_data <- data_list$layout %>% mutate(marker = data_list$graph %>% pull(marker))
    colorscale <- list(c(0, 1), c("lightgrey", color))
    # Plot 3D graph using plotly
    if(!project){
    fig <- plotly::plot_ly(plot_data, 
                           x = ~x, 
                           y = ~y, 
                           z = ~z, 
                           marker = list(color = ~marker, 
                                         colorscale = colorscale, 
                                         showscale = TRUE, 
                                         size = node_size))
    fig <- fig %>% plotly::add_markers() %>% plotly::layout(scene = list(aspectmode=aspectmode,
                                                                         xaxis = list(visible = showgrid),
                                                                         yaxis = list(visible = showgrid),
                                                                         zaxis = list(visible = showgrid)),
                                                            annotations = list(x = 1,
                                                                               y = 0.98,
                                                                               text = marker,
                                                                               showarrow = FALSE))
    } else {
      plot_data <- 
        plot_data %>% 
        
        # Normalize 3D coordinates to a sphere
        mutate(norm_factor = 
                 select(., x, y, z) %>% 
                 apply(MARGIN = 1, function(x) {
                   as.matrix(x) %>% 
                     norm(type = "F")
                 }),
               x_norm = x / norm_factor,
               y_norm = y / norm_factor,
               z_norm = z / norm_factor) 
      fig <- plotly::plot_ly(plot_data, 
                             x = ~x_norm, 
                             y = ~y_norm, 
                             z = ~z_norm, 
                             marker = list(color = ~marker, 
                                           colorscale = colorscale, 
                                           showscale = TRUE, 
                                           size = node_size))
      fig <- fig %>% plotly::add_markers()  %>% plotly::layout(scene = list(xaxis = list(visible = showgrid),
                                                                            yaxis = list(visible = showgrid),
                                                                            zaxis = list(visible = showgrid)),annotations = list(x = 1,
                                                                                  y = 0.98,
                                                                                  text = marker,
                                                                                  showarrow = FALSE))
      warning("Projection to a sphere overrides the aspect ratio setting", call. = FALSE)
    }
    return(fig)
  }
  