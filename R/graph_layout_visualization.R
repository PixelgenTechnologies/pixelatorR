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
#' @param map_nodes,map_edges Should nodes and/or edges be mapped? Note that
#' component graphs can have >100k edges which can be very slow to draw.
#' @param log_scale Convert node counts to log-scale with \code{log1p}. This
#' parameter is ignored for PNA graphs.
#' @param node_size Size of nodes
#' @param edge_width Set the width of the edges if \code{map_edges = TRUE}
#' @param show_Bnodes Should B nodes be included in the visualization?
#' This option is only applicable to bipartite MPX graphs. Note that by removing
#' the B nodes, all edges are removed from the graph and hence, \code{map_edges}
#' will have no effect.
#' @param collect_scales Collect color scales so that their limits are the same.
#' This can be used to make sure that the colors are comparable across markers.
#' @param return_plot_list Instead of collecting the plots in a grid, return a
#' list of \code{ggplot} objects.
#' @param ... Not yet implemented
#'
#' @rdname Plot2DGraph
#'
#' @import rlang
#'
#' @return An object of class \code{patchwork}
#'
#' @seealso [Plot2DGraphM()]
#'
#' @examples
#' library(pixelatorR)
#'
#' # MPX
#' pxl_file <- minimal_mpx_pxl_file()
#' seur <- ReadMPX_Seurat(pxl_file)
#' seur <- LoadCellGraphs(seur, load_as = "Anode")
#' seur <- ComputeLayout(seur, layout_method = "pmds", dim = 2)
#' Plot2DGraph(seur, cells = colnames(seur)[1], layout_method = "pmds", marker = "CD3E")
#'
#' # PNA
#' pxl_file <- minimal_pna_pxl_file()
#' seur <- ReadPNA_Seurat(pxl_file)
#' seur <- LoadCellGraphs(seur, cells = colnames(seur)[1], add_layouts = TRUE)
#' Plot2DGraph(seur, cells = colnames(seur)[1], marker = "CD3e")
#'
#' @export
#'
Plot2DGraph <- function(
  object,
  cells,
  marker = NULL,
  assay = NULL,
  layout_method = c("wpmds_3d", "pmds_3d", "wpmds", "pmds"),
  colors = c("lightgrey", "mistyrose", "red", "darkred"),
  map_nodes = TRUE,
  map_edges = FALSE,
  log_scale = TRUE,
  node_size = 0.5,
  edge_width = 0.3,
  show_Bnodes = TRUE,
  collect_scales = FALSE,
  return_plot_list = FALSE,
  ...
) {
  assert_class(object, "Seurat")
  assert_vector(colors)
  assert_single_value(map_nodes, type = "bool")
  assert_single_value(map_edges, type = "bool")
  if ((map_nodes == FALSE) && (map_edges == FALSE)) {
    cli::cli_abort(
      c("x" = "One of {.var map_nodes} or {.var map_edges} must be TRUE")
    )
  }
  assert_vector(cells, type = "character", n = 1)
  assert_single_value(node_size, type = "numeric")
  assert_single_value(edge_width, type = "numeric")
  assert_single_value(marker, type = "string", allow_null = TRUE)

  # Check and select a layout method
  layout_method <- match.arg(layout_method, choices = c("wpmds_3d", "pmds_3d", "wpmds", "pmds"))
  layout_method_ext <-
    switch(layout_method,
      "pmds" = "pivot MDS (pmds)",
      "wpmds" = "weighted pivot MDS (wpmds)",
      "pmds_3d" = "pivot MDS (pmds)",
      "wpmds_3d" = "weighted pivot MDS (wpmds)"
    )

  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    assert_single_value(assay, type = "string")
  } else {
    assay <- DefaultAssay(object)
  }

  # Validate assay
  pixel_assay <- object[[assay]]
  assert_pixel_assay(pixel_assay)

  # Fetch data
  data_list <- lapply(cells, function(cell_id) {
    # Fetch component graph
    component_graph <- CellGraphs(pixel_assay)[[cell_id]]
    if (is.null(component_graph)) {
      cli::cli_abort(
        c(
          "i" = "All selected cells must have a component graph loaded.",
          " " = "Did you forget to run {.fn LoadCellGraphs}?",
          "x" = "Missing {.cls CellGraph} for component {.str {cell_id}}"
        )
      )
    }

    # unpack values
    graph <- component_graph@cellgraph

    # Validate marker
    if (!is.null(marker)) {
      if (marker == "node_type") {
        if (!attr(graph, "type") == "bipartite") {
          cli::cli_abort(
            c("x" = "{.var marker = 'node_type'} can only be used for bipartite graphs")
          )
        }
      } else {
        if (!marker %in% colnames(component_graph@counts)) {
          cli_alert_danger(
            glue(
              "'{marker}' is missing from node count matrix ",
              "for component {cell_id}"
            )
          )
          return(NULL)
        }
      }
    }

    layout <- component_graph@layout[[layout_method]]
    if (length(graph) == 0) {
      cli::cli_abort(
        c("x" = "Missing cellgraph for component '{.val {cell_id}}'")
      )
    }
    if (length(layout) == 0) {
      cli::cli_abort(
        c("x" = "Missing layout {.str {layout_method}} for component '{.val {cell_id}}'")
      )
    }
    layout <- layout %>% select(x, y)

    # Add node marker counts if needed
    if (!is.null(marker)) {
      if (marker != "node_type") {
        graph <- graph %N>%
          mutate(marker = component_graph@counts[, marker]) %>%
          {
            if (log_scale && !inherits(pixel_assay, c("PNAAssay", "PNAAssay5"))) {
              mutate(., marker = log1p(marker))
            } else {
              .
            }
          }
      }
    }

    # Remove B nodes if show_Bnodes=FALSE
    if ((attr(graph, "type") == "bipartite") && !show_Bnodes && !inherits(pixel_assay, c("PNAAssay", "PNAAssay5"))) {
      inds_keep <- (graph %N>% pull(node_type)) %in% c("A", "umi1")
      graph <- graph %N>%
        filter(node_type %in% c("A", "umi1"))
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
  }) %>%
    set_names(nm = cells)

  # Set limits
  if (!is.null(marker)) {
    if (collect_scales) {
      max_val <- sapply(data_list, function(x) {
        if (is.null(x)) {
          return(NULL)
        }
        x$graph %>% pull(marker)
      }) %>%
        unlist() %>%
        max()
      limits <- c(rep(list(c(0, max_val)), length(data_list))) %>% set_names(names(data_list))
    } else {
      limits <- lapply(data_list, function(x) {
        if (is.null(x)) {
          return(NULL)
        }
        c(0, max(x$graph %>% pull(marker)))
      })
    }
  }

  # Create plots
  plots <- lapply(names(data_list), function(nm) {
    if (is.null(data_list[[nm]])) {
      return(patchwork::plot_spacer() + theme(plot.background = element_rect(fill = NA, colour = NA)))
    }
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
            if (log_scale && !inherits(pixel_assay, c("PNAAssay", "PNAAssay5"))) {
              labs(title = glue("{nm}"), color = paste0(marker, "\n(log-scaled)"))
            } else {
              labs(title = glue("{nm}"), color = marker)
            }
          }
        } else {
          labs(title = glue("{nm}"))
        }
      } +
      theme_void()

    # Add color scale
    if (!is.null(marker)) {
      if (marker != "node_type") {
        p <- p + scale_color_gradientn(colours = colors, limits = limits[[nm]])
      }
    }
    return(p)
  })

  if (return_plot_list) {
    return(plots)
  }

  # Wrap plots
  p <- wrap_plots(plots)
  p <- p + plot_annotation(
    title = glue("Layout with {layout_method_ext}"),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

  return(p)
}


#' Plot multiple markers on multiple graphs
#'
#' In contrast to \code{\link{Plot2DGraph}}, which only draw 1 marker at the time,
#' this function makes it possible to arrange plots into a grid with markers in rows
#' and components in columns. The color scales are fixed for each marker so that their
#' limits are the same across all components.
#'
#' @param markers A character vector with marker names
#' @param titles A named character vector with optional titles. The
#' names of \code{titles} should match \code{cells}
#' @param titles_theme A \code{theme} used to style the titles
#' @param titles_size The size of the text in the plot titles
#' @param titles_col The color of the plot titles
#' @param ... Parameters passed to Plot2DGraph
#' @inheritParams Plot2DGraph
#'
#' @return A \code{patchwork} object
#'
#' @seealso [Plot2DGraph()]
#'
#' @examples
#' library(pixelatorR)
#'
#' # MPX
#' pxl_file <- minimal_mpx_pxl_file()
#' seur <- ReadMPX_Seurat(pxl_file)
#' seur <- LoadCellGraphs(seur, load_as = "Anode")
#' seur <- ComputeLayout(seur, layout_method = "pmds", dim = 2)
#' Plot2DGraphM(seur, cells = colnames(seur)[2:3], layout_method = "pmds", markers = c("CD20", "CD4"))
#'
#' # PNA
#' pxl_file <- minimal_pna_pxl_file()
#' seur <- ReadPNA_Seurat(pxl_file)
#' seur <- LoadCellGraphs(seur, cells = colnames(seur)[2:3], add_layouts = TRUE)
#' Plot2DGraphM(seur, cells = colnames(seur)[2:3], markers = c("CD20", "CD4"))
#'
#' @export
#'
Plot2DGraphM <- function(
  object,
  cells,
  markers,
  assay = NULL,
  layout_method = c("wpmds_3d", "pmds_3d", "wpmds", "pmds"),
  colors = c("lightgrey", "mistyrose", "red", "darkred"),
  map_nodes = TRUE,
  map_edges = FALSE,
  log_scale = TRUE,
  node_size = 0.5,
  edge_width = 0.3,
  show_Bnodes = TRUE,
  titles = NULL,
  titles_theme = NULL,
  titles_size = 10,
  titles_col = "black",
  ...
) {
  assert_vector(cells, type = "character", n = 1)
  assert_vector(markers, type = "character", n = 1)
  assert_single_value(titles_size, type = "numeric")
  assert_single_value(titles_col, type = "string")

  # Set layout options
  ncols <- length(cells)
  nrows <- length(markers)

  # Validate titles and titles_theme
  if (!is.null(titles)) {
    assert_vector(titles, type = "character", n = 1)
    assert_vectors_match(cells, names(titles))
    titles <- titles[cells]
    assert_class(titles_theme, "theme", allow_null = TRUE)
  } else {
    titles <- cells
  }

  # Create columnn titles for patchwork
  title_plots <- lapply(titles, function(label) {
    # Here we create a title for each column
    # titles_size and titles_col controls the size and color of the text
    p_title <- wrap_elements(
      textGrob(
        label = label,
        gp = gpar(fontsize = titles_size, col = titles_col)
      )
    ) +
      theme(plot.background = element_rect(fill = "lightgrey", colour = NA))

    # If a theme is provided for the titles, this can be added here
    if (!is.null(titles_theme)) {
      p_title <- p_title + titles_theme
    }
    return(p_title)
  })

  # Add an empty space for the last column which is reserved
  # for the color legends
  title_plots <- c(title_plots, list(plot_spacer()))

  # Create rows
  plots <- list()
  legends <- list()
  for (marker in markers) {
    # Create a list of plots for each component
    # and the selected marker. By setting collect_scales = TRUE
    # we make sure that the limits of the color scale are the same
    # across components.
    plots_rows <- Plot2DGraph(
      object,
      cells = cells,
      layout_method = layout_method,
      marker = marker,
      colors = colors,
      assay = assay,
      map_nodes = map_nodes,
      map_edges = map_edges,
      log_scale = log_scale,
      node_size = node_size,
      edge_width = edge_width,
      show_Bnodes = show_Bnodes,
      collect_scales = TRUE,
      return_plot_list = TRUE,
      ...
    )

    # Fetch legend and modify themes
    for (i in seq_along(plots_rows)) {
      p <- plots_rows[[i]] + theme(plot.title = element_blank())

      # Fetch color legend from first plot
      # and save it for later
      if (ncols == i) {
        legends[[marker]] <- .get_legend(p) %>% wrap_elements()
      }

      # Remove the color legend and plot title.
      # These will be added to the patchwork in the final step
      p <- p + theme(legend.position = "none", plot.title = element_blank())
      plots_rows[[i]] <- p
    }

    # Each element in plots stores a single row of plots
    plots[[marker]] <- plots_rows
  }

  # Add legends to plots. Each element stored length(cells)
  # plots and after adding the legend, each row will be
  # number of cells + 1
  plots <- lapply(seq_along(plots), function(i) {
    c(plots[[i]], list(legends[[i]]))
  })

  # Combine all plots into a list
  plots <- plots %>% Reduce(c, .)

  # Final plot
  p <- wrap_plots(c(title_plots, plots),
    ncol = ncols + 1, # Sets the number of columns to length(cells) + 1
    heights = c(1, rep(4, nrows)), # Adjust the heights to make the title row a bit smaller
    widths = c(rep(4, ncols), 2)
  ) # Adjust the widths to make the legend column a bit smaller
  return(p)
}

#' Extract legend from a ggplot
#'
#' Adapted from .get_legend in \code{ggpubr}
#'
#' @param p A \code{ggplot} object
#'
#' @noRd
#'
.get_legend <- function(
  p
) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) {
    leg <- tmp$grobs[[leg]]
  } else {
    leg <- NULL
  }
  leg
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
#' @param aspectmode Set aspect ratio to one of "data" or "cube".
#' If "cube", this scene's axes are drawn as a cube, regardless of the axes' ranges.
#' If "data", this scene's axes are drawn in proportion with the axes' ranges.
#'
#' Default "data"
#'
#' @param colors Color the nodes expressing a marker. Must be a character vector
#' with at least two color names.
#' @param log_scale  Convert node counts to log-scale with \code{logp}
#' @param node_size Size of nodes
#' @param show_Bnodes  Should B nodes be included in the visualization?
#' This option is only applicable to bipartite graphs.
#' @param ... Additional parameters passed to \code{plot_ly}
#' @param showgrid Show the grid lines. Default TRUE
#'
#' @rdname Plot3DGraph
#'
#' @return A interactive 3D plot of a component graph layout as a \code{plotly} object
#'
#' @examples
#' library(pixelatorR)
#'
#' # MPX
#' pxl_file <- minimal_mpx_pxl_file()
#' seur <- ReadMPX_Seurat(pxl_file)
#' seur <- LoadCellGraphs(seur, cells = colnames(seur)[5])
#' seur <- ComputeLayout(seur, layout_method = "wpmds", dim = 3, pivots = 50)
#' Plot3DGraph(seur, cell_id = colnames(seur)[5], marker = "CD50", layout_method = "wpmds_3d")
#'
#' # PNA
#' pxl_file <- minimal_pna_pxl_file()
#' seur <- ReadPNA_Seurat(pxl_file)
#' seur <- LoadCellGraphs(seur, cells = colnames(seur)[1], add_layouts = TRUE)
#' Plot3DGraph(seur, cell_id = colnames(seur)[1], marker = "CD16", layout_method = "wpmds_3d")
#'
#' @export
#'
Plot3DGraph <- function(
  object,
  cell_id,
  marker = NULL,
  assay = NULL,
  layout_method = c("wpmds_3d", "pmds_3d"),
  project = FALSE,
  aspectmode = c("data", "cube"),
  colors = c("lightgrey", "mistyrose", "red", "darkred"),
  showgrid = TRUE,
  log_scale = TRUE,
  node_size = 2,
  show_Bnodes = FALSE,
  ...
) {
  # Validate input parameters
  assert_class(object, "Seurat")
  assert_vector(colors, type = "character")
  assert_single_value(cell_id, type = "string")
  if (!(cell_id %in% colnames(object))) {
    cli::cli_abort(
      c(
        "i" = "Please provide a valid {.var cell_id}",
        "x" = "Cell ID {.var {cell_id}} not found in the {.cls Seurat} object"
      )
    )
  }

  assert_single_value(marker, type = "string", allow_null = TRUE)

  # Check and select an aspectmode
  aspectmode <- match.arg(aspectmode, choices = c("data", "cube"))

  # Use default assay if assay = NULL
  if (!is.null(assay)) {
    assert_single_value(assay, type = "string")
  } else {
    assay <- DefaultAssay(object)
  }

  # Validate assay
  pixel_assay <- object[[assay]]
  assert_pixel_assay(pixel_assay)

  # Fetch component graph
  component_graph <- CellGraphs(pixel_assay)[[cell_id]]
  if (is.null(component_graph)) abort(glue("Missing cellgraph for component '{cell_id}'"))

  # unpack values
  graph <- component_graph@cellgraph

  # Validate marker
  if (!is.null(marker)) {
    if (marker == "node_type") {
      if (!attr(graph, "type") == "bipartite") {
        cli::cli_abort(
          c("x" = "{.var marker = 'node_type'} can only be used for bipartite graphs")
        )
      }
    } else {
      if (!marker %in% colnames(component_graph@counts)) {
        cli::cli_abort(
          "marker {.val {marker}} is missing from node count matrix for component {.val {cell_id}}"
        )
      }
    }
  }


  if (!layout_method %in% names(component_graph@layout)) {
    cli::cli_abort(
      c("x" = "Missing layout {.str {layout_method}} for component '{.val {cell_id}}'")
    )
  }
  layout <- component_graph@layout[[layout_method]]

  if (length(graph) == 0) {
    cli::cli_abort(
      c("x" = "Missing cellgraph for component '{.val {cell_id}}'")
    )
  }
  if (ncol(layout) != 3) {
    cli::cli_abort(
      c(
        "x" = "Layout {.val layout_method} must have 3 columns but has {ncol(layout)} columns"
      )
    )
  }

  # Add node marker counts if needed
  if (!is.null(marker)) {
    if (marker != "node_type") {
      layout <- layout %>%
        mutate(marker = component_graph@counts[, marker, drop = TRUE]) %>%
        {
          if (log_scale && !inherits(pixel_assay, c("PNAAssay", "PNAAssay5"))) {
            mutate(., marker = log1p(marker))
          } else {
            .
          }
        }
    }
  }

  # Remove B nodes if show_Bnodes=FALSE
  if ((attr(graph, "type") == "bipartite")) {
    layout$node_type <- graph %N>% pull(node_type)
    if (!show_Bnodes && !inherits(pixel_assay, c("PNAAssay", "PNAAssay5"))) {
      layout <- layout %>%
        filter(node_type %in% c("A", "UMI1", "umi1"))
    }
  }

  # Project data to sphere if project=TRUE
  if (project) {
    # Normalize 3D coordinates to a sphere
    layout_sphere <- layout %>%
      select(x, y, z) %>%
      project_layout_coordinates_on_unit_sphere()
    layout[, c("x", "y", "z")] <- layout_sphere
  }

  # Plot 3D graph using plotly
  fig <- plotly::plot_ly(layout,
    x = ~x,
    y = ~y,
    z = ~z,
    marker = list(size = node_size),
    mode = "scatter3d",
    ...
  )
  if (!is.null(marker)) {
    if (marker == "node_type") {
      fig <- fig %>%
        plotly::add_markers(color = ~node_type, colors = c("blue", "orange"))
    } else {
      fig <- fig %>%
        plotly::add_markers(color = ~marker, colors = colors)
    }
  } else {
    fig <- fig %>%
      plotly::add_markers()
  }

  fig <- fig %>%
    plotly::layout(
      scene = list(
        aspectmode = aspectmode,
        xaxis = list(visible = showgrid),
        yaxis = list(visible = showgrid),
        zaxis = list(visible = showgrid)
      ),
      annotations = list(
        x = 1,
        y = 0.98,
        text = ifelse(!is.null(marker), marker, ""),
        showarrow = FALSE
      )
    )
  return(fig)
}
